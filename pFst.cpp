#include "Variant.h"
#include "split.h"
#include "cdflib.h"
#include "pdflib.h"
#include "var.h"

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
    cerr << "INFO: version ???? " << endl;
}


void printSummary(char** argv){
    cerr << "description: "
         << "A probabilistic approach for detecting differences in allele frequencies between two populations" << endl
         << endl;
    cerr << "usage: " << argv[0] << " [options] --target <sample_list> --background <sample_list> --type <type> [--file <vcf file>]" << endl
         << endl
         << "  <sample_list> is either:" << endl
         << "      (default behavior) a zero based comma separated list of samples indicies corresponding to VCF columns" << endl
         << "      (with -S option)   a comma separated list of sample names" << endl
         << "      (with -F option)   the name of a file listing sample indicies (one per line, # comments allowed)" << endl
         << "      (with -SF option)  the name of a file listing sample names (one per line, # comments allowed)" << endl
         << endl
         << "options:" << endl
         << "    -t, --target <sample_list>      see <sample_list>" << endl
         << "    -t, --background <sample_list>  see <sample_list>" << endl
         << "    -y  --type <string>     genotype likelihood format ; genotypes: GT (implies --counts), GP, GL or PL; pooled: PO" << endl
         << "    -f  --file <string>     a properly formatted VCF file (if omitted, uses stdin)" << endl
         << "    -d  --deltaaf <float>   skip sites where the difference in allele frequencies is less than deltaaf, default is zero" << endl
         << "    -r  --region <string>   a tabix compliant genomic range : seqid or seqid:start-end" << endl
         << "    -c  --counts            use genotype counts rather than genotype likelihoods to estimate parameters, default false" << endl
         << "    -S  --sample-names      all <sample_list> use sample names instead of indicies" << endl
         << "    -F  --sample-files      all <sample_list> are the names of files listing samples" << endl
         << endl
         << "example:  pFst --target 0,1,2,3,4,5,6,7 --background 11,12,13,16,17,19,22 --file my.vcf --deltaaf 0.1 --type PL" << endl
         << endl;
    printVersion() ;
}


double bound(double v){
    if(v <= 0.00001){
        return  0.00001;
    }
    if(v >= 0.99999){
        return 0.99999;
    }
    return v;
}

double logLbinomial(double x, double n, double p){
    double ans = lgamma(n+1)-lgamma(x+1)-lgamma(n-x+1) + x * log(p) + (n-x) * log(1-p);
    return ans;
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
            const string & type, // FORMAT field to check exists and isn't empty or "."; set to "" to not check
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
            continue; // totally missing entry, skip
        }
        if( !type.empty() ){ // check type field exists if it is set
            fmt_it = in_it->second.find(type);
            if( fmt_it == in_it->second.end() ){
                cerr << "FATAL: Cannot find " << type <<" FORMAT value for sample " << *i << endl;
                exit(1);
            }
            const string & tmp = fmt_it->second.front();
            if( tmp.empty() || tmp == "." ){
                continue; // empty or "." value, so skip sample
            }
        }
        // check GT is fully called
        fmt_it = in_it->second.find("GT");
        if( fmt_it == in_it->second.end() ){
            cerr << "FATAL: Cannot find GT FORMAT value for sample " << *i << endl;
            exit(1);
        }else{
            const string & gt = fmt_it->second.front();
            if( gt.find('.') != string::npos ){ // Skip if any allele in genotype is uncalled
                continue;
            }
        }
        // made it this far, so good to add
        out_samples.push_back(in_it->second);
    }
}


int main(int argc, char** argv) {
    // pooled or genotyped
    int pool = 0;

    // the filename
    string filename;

    // set region to scaffold
    string region; 

    // using vcflib; thanks to Erik Garrison 
    VariantCallFile variantFile;

    // deltaaf is the difference of allele frequency we bother to look at 
    string deltaaf ;
    double daf  = 0;

    // use counts instead of likelihoods switch
    int counts = 0;

    // type pooled GL PL
    string type = "NA";

    const struct option longopts[] = {
        {"version",            0, 0, 'v'},
        {"help",               0, 0, 'h'},
        {"counts",             0, 0, 'c'},
        {"file",               1, 0, 'f'},
        {"target",             1, 0, 't'},
        {"background",         1, 0, 'b'},
        {"sample-names",       1, 0, 'S'},
        {"sample-files",       1, 0, 'F'},
        {"deltaaf",            1, 0, 'd'},
        {"region",             1, 0, 'r'},
        {"type",               1, 0, 'y'},
        {0,0,0,0}
    };

    string target_optarg;
    string background_optarg;
    bool tg_by_sample_name = false;
    bool tg_in_files = false;

    int index;
    int iarg=0;

    while(true){
        iarg = getopt_long(argc, argv, "r:d:t:b:f:y:chvSF", longopts, &index);
        if( iarg == -1 ){
            break;
        }

        switch (iarg)
        {

            case 'v':
                printVersion();
                return 0;
            case 'y':     
                type = optarg;
                cerr << "INFO: genotype likelihoods set to: " << type << endl;
                if(type == "GT"){
                    cerr << "INFO: using counts flag as GT was specified" << endl;
                }
                break;
            case 'c':
                cerr << "INFO: using genotype counts rather than genotype likelihoods" << endl;
                counts = 1;
                break;
            case 't':
                target_optarg = optarg;
                break;
            case 'b':
                background_optarg = optarg;
                break;
            case 'S':
                tg_by_sample_name = true;
                break;
            case 'F':
                tg_in_files = true;
                break;
            case 'f':
                cerr << "INFO: File: " << optarg  <<  endl;
                filename = optarg;
                break;
            case 'd':
                cerr << "INFO: only scoring sites where the allele frequency difference is greater than: " << optarg << endl;
                deltaaf = optarg;
                daf = atof(deltaaf.c_str());        
                break;
            case 'r':
                cerr << "INFO: set seqid region to : " << optarg << endl;
                region = optarg; 
                break;

            case 'h':
            default:
                printSummary(argv);
                return 0;
        }
    }

    // open variant file
    try{
        if( filename.empty() ){
            variantFile.open(cin);
        }else{
            variantFile.open(filename);
        }
    }catch( const std::out_of_range& oor ){ // variantFile.open() throws out_of_range exception when file does not exit
        cerr << "FATAL: failed to open variant file '" << filename << "'" << endl;
        exit(1);
    }
    if( !variantFile.is_open() ){
        cerr << "FATAL: failed to open variant file '" << filename << "'" << endl;
        exit(1);
    }
    vector<string> samples = variantFile.sampleNames;
    int nsamples = samples.size();

    // parse the target and background lists
    vector<string> target_names;
    if( target_optarg.empty() ){
        cerr << "FATAL: Must provide targets" << endl;
    }else{
        loadSampleList(target_optarg, tg_in_files, tg_by_sample_name, variantFile, target_names);
    }
    vector<string> background_names;
    if( background_optarg.empty() ){
        cerr << "FATAL: Must provide background" << endl;
    }else{
        loadSampleList(background_optarg, tg_in_files, tg_by_sample_name, variantFile, background_names);
    }

    // make a list of samples which are in target but not in background (used to get total list)
    vector<string> t_not_in_bg_names;
    for( vector<string>::const_iterator i=target_names.begin(); i != target_names.end(); ++i ){
        if( find(background_names.begin(), background_names.end(), *i) == background_names.end() ){
            t_not_in_bg_names.push_back(*i);
        }
    }

    // verbose output of sample names
    cerr << "INFO: target argument: " << target_optarg << endl;
    cerr << "INFO: target samples (" << target_names.size() << "): ";
    printVector(cerr, target_names) << endl;
    cerr << "INFO: background argument: " << target_optarg << endl;
    cerr << "INFO: background samples (" << background_names.size() << "): ";
    printVector(cerr, background_names) << endl;
    cerr << "INFO: total samples used (" << background_names.size()+t_not_in_bg_names.size() << "): ";
    printVector(cerr, t_not_in_bg_names);
    printVector(cerr, background_names) << endl;

    if( !region.empty() ){
        if(! variantFile.setRegion(region)){
            cerr <<"FATAL: unable to set region" << endl;
            return 1;
        }
    }

    string okayGenotypeLiklihoods[] = {"GT","PL","PO","GL","GP"};
    if( type.empty() ){
        cerr << "FATAL: failed to specify genotype likelihood format : GT,PL,PO,GL,GP" << endl;
        printSummary(argv);
        return 1;
    }
    if(type == "GT"){ // --type=GT implies --counts
        counts = 1;
    }
    string * tmp_end = okayGenotypeLiklihoods+(sizeof(okayGenotypeLiklihoods)/sizeof(char*));
    if( find(okayGenotypeLiklihoods, tmp_end, type) == tmp_end ){
        cerr << "FATAL: genotype likelihood is incorrectly formatted, only use: GT,PL,PO,GL,GP" << endl;
        printSummary(argv);
        return 1;
    }    

    Variant var(variantFile);

    while (variantFile.getNextVariant(var)) {

        if(var.alt.size() > 1){
            continue;
        }

        vector < map< string, vector<string> > > target; 
        vector < map< string, vector<string> > > background; 
        vector < map< string, vector<string> > > total; 

        getNonNullSamples(target_names, var.samples, type, target);
        getNonNullSamples(background_names, var.samples, type, background);
        total = background;
        getNonNullSamples(t_not_in_bg_names, var.samples, type, total, false);


        zvar * populationTarget        ;
        zvar * populationBackground    ;
        zvar * populationTotal         ;

        if(type == "PO"){
            populationTarget     = new pooled();
            populationBackground = new pooled();
            populationTotal      = new pooled();
        }
        if(type == "PL"){
            populationTarget     = new pl();
            populationBackground = new pl();
            populationTotal      = new pl();
        }
        if(type == "GL"){
            populationTarget     = new gl();
            populationBackground = new gl();
            populationTotal      = new gl();    
        }
        if(type == "GP"){
            populationTarget     = new gp();
            populationBackground = new gp();
            populationTotal      = new gp();
        }
        if(type == "GT"){
            populationTarget     = new gt();
            populationBackground = new gt();
            populationTotal      = new gt();
        }   

        populationTotal->loadPop(total          , var.sequenceName, var.position);  
        populationTarget->loadPop(target        , var.sequenceName, var.position);
        populationBackground->loadPop(background, var.sequenceName, var.position);

        if(populationTarget->npop < 2 || populationBackground->npop < 2){
            delete populationTarget;
            delete populationBackground;
            delete populationTotal;
            continue;
        }

        populationTotal->estimatePosterior();   
        populationTarget->estimatePosterior();
        populationBackground->estimatePosterior();

        if(populationTarget->alpha == -1 || populationBackground->alpha == -1){
            delete populationTarget;
            delete populationBackground;
            delete populationTotal;
            continue;
        }

        if(counts == 1){

            populationTotal->alpha  = 0.001 + populationTotal->nref;
            populationTotal->beta   = 0.001 + populationTotal->nalt;

            populationTarget->alpha = 0.001 + populationTarget->nref;
            populationTarget->beta  = 0.001 + populationTarget->nalt;

            populationBackground->alpha = 0.001 + populationBackground->nref;
            populationBackground->beta  = 0.001 + populationBackground->nalt;
        }

        double populationTotalEstAF       = bound(populationTotal->beta      / (populationTotal->alpha      + populationTotal->beta)     );
        double populationTargetEstAF      = bound(populationTarget->beta     / (populationTarget->alpha     + populationTarget->beta)    );
        double populationBackgroundEstAF  = bound(populationBackground->beta / (populationBackground->alpha + populationBackground->beta));

        // cout << populationTotalEstAF << "\t" << populationTotal->af << endl;

        // x, n, p
        double null = logLbinomial(populationTarget->beta, (populationTarget->alpha + populationTarget->beta),  populationTotalEstAF) +
            logLbinomial(populationBackground->beta, (populationBackground->alpha + populationBackground->beta),  populationTotalEstAF) ;
        double alt  = logLbinomial(populationTarget->beta, (populationTarget->alpha + populationTarget->beta),  populationTargetEstAF) +
            logLbinomial(populationBackground->beta, (populationBackground->alpha + populationBackground->beta),  populationBackgroundEstAF) ;

        double l = 2 * (alt - null);

        if(l <= 0){
            delete populationTarget;
            delete populationBackground;
            delete populationTotal;
            continue;
        }

        int     which = 1;
        double  p ;
        double  q ;
        double  x  = l;
        double  df = 1;
        int     status;
        double  bound ;

        cdfchi(&which, &p, &q, &x, &df, &status, &bound );

        cout << var.sequenceName << "\t"  << var.position << "\t" << 1-p << endl ;

        delete populationTarget;
        delete populationBackground;
        delete populationTotal;

        populationTarget     = NULL;
        populationBackground = NULL;
        populationTotal      = NULL;

    }
    return 0;           
}
