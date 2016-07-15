# GLPhase

This is a cuda-enabled fork of
[SNPTools impute.cpp](http://sourceforge.net/p/snptools/code/ci/master/tree/). This
code should scale linearly with sample size up to a small multiple of
the number of CUDA cores (shaders) on the GPU being used. 

GLPhase also has an option for incorporating pre-existing haplotypes
into the phasing and imputation
process. [Release 1.4.13](https://github.com/wkretzsch/GLPhase/releases/tag/v1.4.13)
was used with this option 
to impute genotypes for the first release of the 
[Haplotype Reference Consortium](http://www.haplotype-reference-consortium.org/).

## Installation

### Dependencies

GLphase depends on [libgsl](https://www.gnu.org/software/gsl/) and [boost](http://www.boost.org/).

### Compilation

    # to compile all code (with all optimizations turned on)
    make

    # run the glphase executable to get a description of the
    # glphase command line arguments
    bin/glphase

    # run regression tests (turns off optimizations)
    make test

    # run regression tests + longer integration tests
    make disttest

    # compile without CUDA support
    # first clean the work dir
    make clean
    make NCUDA=1

    # compile without CUDA or OMP support (on MacOSX for example)
    make NCUDA=1 NOMP=1
   
## Running GLPhase (v1.4.13)


