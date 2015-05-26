# GPAT
A version Zev Kronenberg's GPAT++ (Genotype-Phenotype Association Toolkit)  
See: https://github.com/jewmanchue/vcflib/wiki  
And the original repository: https://github.com/jewmanchue/vcflib  

I've separated out GPAT++ from ekg's original vcflib code.
It still depends in vcflib, but can be maintained independently and built using updated versions of vcflib.  
In fact, GPAT can be included as a submodule to vcflib, as my fork of vcflib does: https://github.com/travc/vcflib

## Install ##
- get ekg's vcflib https://github.com/ekg/vcflib  
`git clone --recursive https://github.com/ekg/vcflib.git`
- `cd vcflib`
- make vcflib  
`make`
- get my GPAT  
`git clone https://github.com/travc/GPAT.git`
- make GPAT
`cd GPAT`  
`make`

That should just work. 
**The executables will all be put into `vcflib/bin`.**

You could also put GPAT wherever you want.
Just set `VCFLIB_PATH` in the `Makefile` to point to the location of ekg's vcflib.
(You need a copy/clone of vcflib regardless.)
