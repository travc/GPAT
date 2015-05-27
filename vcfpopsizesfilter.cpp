#include "Variant.h"
#include "split.h"
#include "cdflib.h"
#include "pdflib.h"
#include "var.h"

#include <assert.h>
#include <string>
#include <iostream>
#include <math.h>  
#include <cmath>
#include <stdexcept>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <getopt.h>

using namespace std;
using namespace vcflib;

void printVersion(void){
    //cerr << "INFO: version 1.0.0 ; date: April 2014 ; author: Zev Kronenberg; email : zev.kronenberg@utah.edu " << endl;
    // @TCC TODO
    cerr << "version ???? " << endl;
}


void printSummary(char** argv){
    cerr << "description: "
         << "Remove lines from a vcf if any population has less than requested number of called genotypes " << endl
         << endl;
    cerr << "usage: " << argv[0] << " [options] [<vcf_file>]" << endl
         << "    -h, --help                 display this message and exit" << endl
         << "    -p, --pop <sample_list>    (can give multiple times) comma separated list of sample names" << endl
         << "    -P, --pop-file <filename>  (can give multiple times) file listing sample names, one per line, '#' comments allowed" << endl
         << "    -N, --min-size <int>       (can give multiple times) minimum number of non-null samples require in a pop to pass" << endl
         << "                               NOTE: must have a -N for each pop (-p or -P) given" << endl
         << "    -y, --fmt <format_tag>     (can give multiple times) genotype FORMAT field to require; defualt is \"GT\"" << endl
         << "    -r  --region <string>      a tabix compliant genomic range : seqid or seqid:start-end" << endl
         << "    -v, --verbose              (can give multiple times) increase (stderr) output verbosity" << endl
         << "    -V, --version              output version number (TODO)" << endl
         << endl
         << "  If <vcf_file> is omitted, reads from stdin" <<endl
         << endl
         << "example: "<< argv[0] <<" -p 'samp1,samp2,samp3' -N 2 -P pop2_samps.txt -N 5" << endl
         << endl;
    printVersion() ;
}


template <typename T>
ostream & printVector(ostream & out, const vector<T> & v, const string & delim=" "){
    for( size_t i=0; i<v.size(); ++i ){
        if( i != 0 ){
            out << delim;
        }
        out << v[i];
    }
    return out;
}


void loadSampleList(const string & in_string, const bool string_is_filename, const bool samples_by_name, const VariantCallFile & variantFile, 
        vector<string> & samp_names_out){
    vector<string> indviduals;
    if( !string_is_filename ){
        indviduals = split(in_string, ",");
    }else{
        // read individuals from file
        string line;
        ifstream fs(in_string.c_str());
        while( getline(fs, line) ){
            // strip '#' comments
            line = line.substr(0, line.find_first_of('#'));
            // trim trailing whitespace
            line.erase(find_if(line.rbegin(), line.rend(), not1(ptr_fun<int, int>(isspace))).base(), line.end());
            if( !line.empty() ){ // skip empty lines
                indviduals.push_back(line);
            }
        }
        fs.close();
    }
    int idx;
    for( vector<string>::const_iterator i = indviduals.begin(); i != indviduals.end(); i++){
        if( !samples_by_name ){ // list of indicies
            idx = atoi( (*i).c_str() );
            if( idx >= variantFile.sampleNames.size() ){
                cerr << "FATAL: request sample index "<< idx << " but only " << variantFile.sampleNames.size() << " samples in input" << endl;
                exit(1);
            }
            samp_names_out.push_back(variantFile.sampleNames[idx]);
        }else{ // list of sample names
            samp_names_out.push_back(*i);
            vector<string>::const_iterator tmp_it = find(variantFile.sampleNames.begin(), variantFile.sampleNames.end(), *i);
            if( tmp_it == variantFile.sampleNames.end() ){
                cerr << "FATAL: cannot find sample '" << *i << "' in input" << endl;
                exit(1);
            }
        }
    }
}

void getNonNullSamples(const vector<string> & sample_names, 
            const map<string, map<string, vector<string> > > & in_samples, 
            const vector<string> & required_fmts, // FORMAT fields to check exists and isn't empty or "." (or contains '.' in the case of GT)
            vector<map<string, vector<string> > > & out_samples,
            bool clear_out_flag=false){
    if( clear_out_flag ){
        out_samples.clear();
    }
    map<string, map<string, vector<string> > >::const_iterator in_it;
    map<string, vector<string> >::const_iterator fmt_it;
    for( vector<string>::const_iterator i = sample_names.begin(); i != sample_names.end(); ++i ){
        in_it = in_samples.find(*i);
        if( in_it == in_samples.end() ){
            goto sample_failed; // totally missing entry, skip
        }
        for( size_t req_fmt_i=0; req_fmt_i<required_fmts.size(); ++req_fmt_i ){
            const string & req_fmt = required_fmts[req_fmt_i];
            fmt_it = in_it->second.find(req_fmt);
            if( fmt_it == in_it->second.end() ){
                cerr << "FATAL: Cannot find " << req_fmt <<" FORMAT value for sample " << *i << endl;
                exit(1);
            }
            const string & val = fmt_it->second.front();
            if( req_fmt == "GT" ){
                if( val.find('.') != string::npos ){ // skip if any allele in genotype is uncalled
                    goto sample_failed;
                }
            }else{
                if( val.empty() || val == "." ){ // empty or "." value, so skip sample
                    goto sample_failed; 
                }
            }
        }
        // made it this far, so good to add
        out_samples.push_back(in_it->second);
        sample_failed: (void)0;  // target of skipping sample goto
    }
}


int main(int argc, char** argv) {

    VariantCallFile variantFile; // using vcflib; thanks to Erik Garrison 

    // args/options
    int verbose;
    string vcf_filename;
    string region;
    vector<string> required_fmts;
    vector<string> pop_args;
    vector<bool> pop_arg_is_file_flags;
    vector< vector<string> > pops;
    vector<int> min_pop_sizes;

    const struct option longopts[] = {
        {"version",            0, 0, 'V'},
        {"help",               0, 0, 'h'},
        {"verbose",            0, 0, 'v'},
        {"pop",                1, 0, 'p'},
        {"pop-file",           1, 0, 'P'},
        {"min-size",           1, 0, 'N'},
        {"fmt",                1, 0, 'y'},
        {"region",             0, 0, 'r'},
        {0,0,0,0}
    };

    int option;
    while(true){
        int option_index = 0;
        option = getopt_long(argc, argv, "hVvp:P:N:y:r", longopts, &option_index);
        if( option == -1 ){
            break;
        }

        switch( option ){

            case 'V':
                printVersion();
                return 0;
            case 'v':
                ++verbose;
                break;

            case 'y':     
                required_fmts.push_back(optarg);
                break;

            case 'p':
                pop_args.push_back(optarg);
                pop_arg_is_file_flags.push_back(false);
                break;
            case 'P':
                pop_args.push_back(optarg);
                pop_arg_is_file_flags.push_back(true);
                break;
            case 'N':
                min_pop_sizes.push_back(atoi(optarg));        
                break;
            case 'r':
                region = optarg;
                break;

            case 'h':
            default:
                printSummary(argv);
                return 0;
        }
    }
    // remaining arguments (maybe vcf_filename)
    while( optind < argc ){
        if( argc - optind > 1 ){ // only allow one argument
            printSummary(argv);
            cerr<<"ERROR: Can only give a single vcf filename"<<endl;
            exit(1);
        }
        vcf_filename = argv[optind++];
    }

    // some more default option values
    if( required_fmts.empty() ){
        required_fmts.push_back("GT");
    }
    

    // option/argmunet sanity checking


    cerr<<"required_fmts : ";
    printVector(cerr, required_fmts)<<endl;
    cerr<<"pop_args : ";
    printVector(cerr, pop_args)<<endl;

    // open variant file
    try{
        if( vcf_filename.empty() ){
            variantFile.open(cin);
        }else{
            variantFile.open(vcf_filename);
        }
    }catch( const std::out_of_range& oor ){ // variantFile.open() throws out_of_range exception when file does not exit
        cerr << "FATAL: failed to open variant file '" << vcf_filename << "'" << endl;
        exit(1);
    }
    if( !variantFile.is_open() ){
        cerr << "FATAL: failed to open variant file '" << vcf_filename << "'" << endl;
        exit(1);
    }

    // load the pops (lists of sample names per population)
    assert( pop_args.size() == pop_arg_is_file_flags.size() );
    for( size_t i=0; i < pop_args.size(); ++i ){
        pops.push_back(vector<string>());
        loadSampleList(pop_args[i], pop_arg_is_file_flags[i], true, variantFile, pops.back());
    }

    // must have a -N for each pop
    if( min_pop_sizes.size() != pops.size() ){
        printSummary(argv);
        cerr << "FATAL: must give min-size (-N) for each pop (-p or -P)" << endl;
        exit(1);
    }
    
    // informative output
    for( size_t i=0; i < pops.size(); ++i ){
        cerr << "POP " << i << ": ";
        printVector(cerr, pops[i])<<endl;
        cerr << "  num samples: " << pops[i].size() << endl;
        cerr << "  min size: " << min_pop_sizes[i] << endl;
    }

    // set region (if given as an option)
    if( !region.empty() ){
        if( !variantFile.setRegion(region) ){
            cerr <<"FATAL: unable to set region '" << region << "'" << endl;
            return 1;
        }
    }

    // add commandline comment to header
    ostringstream ss;
    ss << "##commandline=\"";
    ss << argv[0];
    for( int i=1; i<argc; ++i ){
        ss << " " << argv[i];
    }
    ss << "\"";
    variantFile.addHeaderLine(ss.str());
    // output header
    cout << variantFile.header << endl;

    // main loop
    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {

        // skip monomorphic?
        //if(var.alt.size() > 1){
        //   goto failed_variant_line;
        //}

        vector < map< string, vector<string> > > samples; 
        for( size_t i=0; i < pops.size(); ++i ){
            samples.clear();
            getNonNullSamples(pops[i], var.samples, required_fmts, samples);
            if( samples.size() < min_pop_sizes[i] ){
                goto failed_variant_line; // a population failed, skip output
            }
        }

        cout << var << endl;
        continue;

        failed_variant_line: void(0); // target of skipping variant goto
        if( verbose > 0 ){
            cerr<<var.sequenceName<<"\t"<<var.position<<"\tFailed"<<endl;;
        }
    }
    return 0;           
}
