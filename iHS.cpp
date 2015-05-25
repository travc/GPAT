#include "Variant.h"
#include "split.h"
#include "cdflib.h"
#include "pdflib.h"
#include "var.h"

#include <string>
#include <iostream>
#include <math.h>  
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <getopt.h>

using namespace std;
using namespace vcflib;

void printVersion(void){
	    cerr << "INFO: version 1.1.0 ; date: April 2014 ; author: Zev Kronenberg; email : zev.kronenberg@utah.edu " << endl;
	    exit(1);

}

void printHelp(void){
  cerr << endl << endl;
  cerr << "INFO: help" << endl;
  cerr << "INFO: description:" << endl;
  cerr << "     iHS calculates the integrated ratio of haplotype decay between the reference and non-reference allele. " << endl;
  cerr << "     This implementation of iHS integrates over a SNP index, NOT a physical distance or genetic distance.   " << endl << endl;

  cerr << "Output : 4 columns :                  "    << endl;
  cerr << "     1. seqid                         "    << endl;
  cerr << "     2. position                      "    << endl;
  cerr << "     3. target allele frequency       "    << endl;
  cerr << "     4. integrated EHH (alternative)  "    << endl;
  cerr << "     5. integrated EHH (reference)    "    << endl;
  cerr << "     6. iHS log(iEHHalt/iEHHref)      "    << endl  << endl;

  cerr << "INFO: iHS  --target 0,1,2,3,4,5,6,7 --file my.phased.vcf  --region chr1:1-1000 " << endl << endl;
 
  cerr << "INFO: required: t,target  -- argument: a zero base comma separated list of target individuals corrisponding to VCF columns " << endl;
  cerr << "INFO: required: f,file    -- argument: proper formatted and phased VCF.                                                    " << endl;
  cerr << "INFO: required: y,type    -- argument: genotype likelihood format: PL,GL,GP                                                " << endl;
  cerr << "INFO: optional: r,region  -- argument: a tabix compliant genomic range : \"seqid:start-end\" or \"seqid\"                  " << endl; 
  cerr << endl;
 
  printVersion();

  exit(1);
}

void clearHaplotypes(string haplotypes[][2], int ntarget){
  for(int i= 0; i < ntarget; i++){
    haplotypes[i][0].clear();
    haplotypes[i][1].clear();
  }
}

void loadIndices(map<int, int> & index, string set){
  
  vector<string>  indviduals = split(set, ",");
  vector<string>::iterator it = indviduals.begin();
  
  for(; it != indviduals.end(); it++){
    index[ atoi( (*it).c_str() ) ] = 1;
  }
}

void calc(string haplotypes[][2], int nhaps, vector<double> afs, vector<long int> pos, vector<int> & target, vector<int> & background, string seqid){

   //cerr << "about to calc" << endl;
   //cerr << "nhap:        " << nhaps << endl;

  for(int snp = 0; snp < haplotypes[0][0].length(); snp++){
    
    int breakflag = 1;
    int count     = 0;

    double ehhA = 1;
    double ehhR = 1;

    double iHSA = 0;
    double iHSR = 0;

    int start = snp;
    int end   = snp;
    int core  = snp; 

    while( breakflag ) {
     
      if(count > 0){
	start -= 1;
	end   += 1;
      }
      if(start == -1){
	break;
      }
      if(end == haplotypes[0][0].length() - 1){
	break;
      }
      count += 1;

      map<string , int> targetH;

      double sumrT = 0;
      double sumaT = 0;
      double nrefT = 0;
      double naltT = 0;

      for(int i = 0; i < nhaps; i++){
	targetH[ haplotypes[i][0].substr(start, (end - start)) ]++;
	targetH[ haplotypes[i][1].substr(start, (end - start)) ]++;
      }     
      for( map<string, int>::iterator th = targetH.begin(); th != targetH.end(); th++){        
	if( (*th).first.substr((end-start)/2, 1) == "1"){     
	  sumaT += r8_choose(th->second, 2);  
	  naltT += th->second;
	}
	else{
	  sumrT += r8_choose(th->second, 2);  
	  nrefT += th->second;
	}
      }
      cerr << pos[snp] << " " << naltT  << " " << nrefT << endl; 
      
      double ehhAC = sumaT / (r8_choose(naltT, 2));
      double ehhRC = sumrT / (r8_choose(nrefT, 2));

      if(count == 1){
	ehhA = ehhAC;
	ehhR = ehhRC;
	continue;
      }

      iHSA += (ehhA + ehhAC) / 2;
      iHSR += (ehhR + ehhRC) / 2;

      ehhA = ehhAC;
      ehhR = ehhRC;

      if(ehhA < 0.05 && ehhR < 0.05){
	breakflag = 0;
      }
    } 

    cout << seqid << "\t" << pos[snp] << "\t" << afs[snp] << "\t" << iHSA << "\t" << iHSR << "\t" << log(iHSA/iHSR) << endl;
  }   
}

double EHH(string haplotypes[][2], int nhaps){

  map<string , int> hapcounts;

  for(int i = 0; i < nhaps; i++){
    hapcounts[ haplotypes[i][0] ]++;
    hapcounts[ haplotypes[i][1] ]++;
  }

  double sum = 0;
  double nh  = 0;

  for( map<string, int>::iterator it = hapcounts.begin(); it != hapcounts.end(); it++){
    nh  += it->second; 
    sum += r8_choose(it->second, 2);
  }

  double max = (sum /  r8_choose(nh, 2));
  
  return max;

}

void loadPhased(string haplotypes[][2], genotype * pop, int ntarget){
  
  int indIndex = 0;

  for(vector<string>::iterator ind = pop->gts.begin(); ind != pop->gts.end(); ind++){
    string g = (*ind);
    vector< string > gs = split(g, "|");
    haplotypes[indIndex][0].append(gs[0]);
    haplotypes[indIndex][1].append(gs[1]);
    indIndex += 1;
  }
}

int main(int argc, char** argv) {

  // set the random seed for MCMC

  srand((unsigned)time(NULL));

  // the filename

  string filename = "NA";

  // set region to scaffold

  string region = "NA"; 

  // using vcflib; thanks to Erik Garrison 

  VariantCallFile variantFile;

  // zero based index for the target and background indivudals 
  
  map<int, int> it, ib;
  
  // deltaaf is the difference of allele frequency we bother to look at 

  // ancestral state is set to zero by default


  int counts = 0;
  
  // phased 

  int phased = 0;

  string type = "NA";

    const struct option longopts[] = 
      {
	{"version"   , 0, 0, 'v'},
	{"help"      , 0, 0, 'h'},
        {"file"      , 1, 0, 'f'},
	{"target"    , 1, 0, 't'},
	{"region"    , 1, 0, 'r'},
	{"type"      , 1, 0, 'y'},
	{0,0,0,0}
      };

    int findex;
    int iarg=0;

    while(iarg != -1)
      {
	iarg = getopt_long(argc, argv, "y:r:d:t:b:f:hv", longopts, &findex);
	
	switch (iarg)
	  {
	  case 'h':
	    printHelp();
	  case 'v':
	    printVersion();
	  case 'y':
	    type = optarg;
	    break;
	  case 't':
	    loadIndices(it, optarg);
	    cerr << "INFO: there are " << it.size() << " individuals in the target" << endl;
	    cerr << "INFO: target ids: " << optarg << endl;
	    break;
	  case 'f':
	    cerr << "INFO: file: " << optarg  <<  endl;
	    filename = optarg;
	    break;
	  case 'r':
            cerr << "INFO: set seqid region to : " << optarg << endl;
	    region = optarg; 
	    break;
	  default:
	    break;
	  }
      }

    map<string, int> okayGenotypeLikelihoods;
    okayGenotypeLikelihoods["PL"] = 1;
    okayGenotypeLikelihoods["GL"] = 1;
    okayGenotypeLikelihoods["GP"] = 1;
    okayGenotypeLikelihoods["GT"] = 1;
    

    if(type == "NA"){
      cerr << "FATAL: failed to specify genotype likelihood format : PL or GL" << endl;
      printHelp();
      return 1;
    }
    if(okayGenotypeLikelihoods.find(type) == okayGenotypeLikelihoods.end()){
      cerr << "FATAL: genotype likelihood is incorrectly formatted, only use: PL or GL" << endl;
      printHelp();
      return 1;
    }

    if(filename == "NA"){
      cerr << "FATAL: did not specify a file" << endl;
      printHelp();
      return(1);
    }

    if(it.size() < 2){
      cerr << "FATAL: target option is required -- or -- less than two individuals in target\n";
      printHelp();
      return(1);
    }

    variantFile.open(filename);
    
    if(region != "NA"){
      if(! variantFile.setRegion(region)){
	cerr <<"FATAL: unable to set region" << endl;
	return 1;
      }
    }

    
    if (!variantFile.is_open()) {
        return 1;
    }
    
    Variant var(variantFile);

    vector<int> target_h, background_h;


    int index   = 0; 
    int  indexi = 0;


    vector<string> samples = variantFile.sampleNames;
    int nsamples = samples.size();

    for(vector<string>::iterator samp = samples.begin(); samp != samples.end(); samp++){
      
      string sampleName = (*samp);
     
      if(it.find(index) != it.end() ){
	target_h.push_back(indexi);
	indexi++;
      }
      index++;
    }
    
    
    // cerr << "n in target : " << target_h.size() << endl;

    vector<long int> positions;
    
    vector<double> afs;

    string haplotypes [target_h.size()][2];    
    
    string currentSeqid = "NA";
    
    // cerr << "about to loop variants" << endl;

    while (variantFile.getNextVariant(var)) {

      if(!var.isPhased()){
	cerr << "FATAL: Found an unphased variant. All genotypes must be phased!" << endl;
	return(1);
      }

      if(var.alt.size() > 1){
	continue;
      }

      if(currentSeqid != var.sequenceName){
	if(haplotypes[0][0].length() > 10){
	  calc(haplotypes, target_h.size(), afs, positions, target_h, background_h, currentSeqid);
	}
	clearHaplotypes(haplotypes, target_h.size());
	positions.clear();
	currentSeqid = var.sequenceName;
	afs.clear();
      }


      vector < map< string, vector<string> > > target, background, total;
      
      int sindex = 0;
      
      for(int nsamp = 0; nsamp < nsamples; nsamp++){

	map<string, vector<string> > sample = var.samples[ samples[nsamp]];
	
	if(it.find(sindex) != it.end() ){
	  target.push_back(sample);
	}	
	sindex += 1;
      }
      
      genotype * populationTarget    ;

      
      if(type == "PL"){
	populationTarget     = new pl();
      }
      if(type == "GL"){
	populationTarget     = new gl();
      }
      if(type == "GP"){
	populationTarget     = new gp();
      }
      if(type == "GT"){
	populationTarget     = new gt();
      }

      populationTarget->loadPop(target, var.sequenceName, var.position);
      
      if(populationTarget->af > 0.95 || populationTarget->af < 0.05){
	delete populationTarget;
	continue;
      }
      positions.push_back(var.position);
      afs.push_back(populationTarget->af);
      loadPhased(haplotypes, populationTarget, populationTarget->gts.size()); 
    
      populationTarget = NULL;
      delete populationTarget;
    }
    
    calc(haplotypes, target_h.size(), afs, positions, target_h, background_h, currentSeqid);
    
    return 0;		    
}
