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

GLphase depends on [libgsl](https://www.gnu.org/software/gsl/),
[boost](http://www.boost.org/), and [libz](http://www.zlib.net/).

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

## Converting a VCF to SNPTools `.bin` format

A perl script at `scripts/vcf2STBin.pl` can be used to convert a VCF
with PL format fields to a SNPTools conformant `.bin` file.  For
example, this command will convert a gzipped input VCF at
`input.vcf.gz` into a SNPTools `.bin` file at `input.bin`:

    scripts/vcf2STbin.pl input.vcf.gz

## Running GLPhase (v1.4.13)

### As a drop-in replacement for SNPTools/impute.cpp

GLPhase can be run as a CUDA-enabled drop-in replacement for
`SNPTools/impute.cpp`. Assuming a SNPTools style `.bin` file with
genotype likelihoods exists:

    bin/glphase input.bin

### Using pre-existing haplotypes

GLPhase can use pre-existing haplotypes to restrict the set of
possible haplotypes from which the MH sampler may choose surrogate
parent haplotypes. This approach is described in:

> The Haplotype Reference Consortium. A reference panel of 64,976
> haplotypes for genotype imputation. Nature Genetics (accepted) -- 
> [bioRxiv](http://biorxiv.org/content/early/2015/12/23/035170)


This command phases and imputes haplotypes on a SNPTools `.bin` file
using a genetic map and pre-existing haplotypes.  The output file is
a gzipped VCF file at `output_base_name.vcf.gz`.

```bash
glphase -B0 -i5 -m95 -q0 -Q1 -t2 -C100 -K200 \
    input.bin \
    -g genetic_map.txt \
    -h pre_existing_haplotypes.haps.gz \
    -s pre_existing_haplotypes.sample \
    -o output_base_name
```

The pre-existing haplotypes should be in
[WTCCCformat](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#hapsample), 
and a genetic map can be obtained from the [Impute2 website](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#reference).

## Ligating haplotypes

It is recommended to ligate haplotypes using
[hapfuse](https://bitbucket.org/wkretzsch/hapfuse/src).  Before
fusing, the output from GLPhase needs to be converted from gzipped VCF
to something htslib can read. Here an example using [bcftools](https://samtools.github.io/bcftools/bcftools.html):

    zcat output_base_name.vcf.gz | bcftools -Ob -o \
        output_base_name.bcf

