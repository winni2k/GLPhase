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

If you have [Singularity](https://www.sylabs.io/singularity/) version 3 installed, then you can run the glphase container located [here](https://cloud.sylabs.io/library/wkretzsch/default/glphase).

### Dependencies

GLPhase depends on [libgsl](https://www.gnu.org/software/gsl/),
[boost](http://www.boost.org/), and [libz](http://www.zlib.net/).

### Compilation

    # Clone this repository recursively
    git clone --recursive https://github.com/winni2k/GLPhase.git
    cd GLPhase

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
[WTCCC format](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#hapsample), 
and a genetic map can be obtained from the [Impute2 website](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#reference).

### Using a reference panel

GLPhase can use a reference panel of haplotypes to inform genotype
imputation of samples for which genotype likelihoods are available.
In contrast to pre-existing haplotypes, the haplotypes
in the reference panel do not need to be from the same samples that
are being imputed. In this mode, when surrogate parent haplotypes 
are being chosen for a sample, the haplotypes may come from the 
current estimate of sample haplotypes or the reference panel. `-k` 
can be specified to restrict the choice of surrogate parent haplotypes
to the reference panel in the first iteration of haplotype estimation.

```bash
glphase \
    input.bin \
    -g samples/hapGen/ex.map \
    -H samples/hapGen/ex.haps.gz \
    -L samples/hapGen/ex.leg \
    -k \
    -o output_base_name
```

The reference haplotypes and legend should be in
[Impute2 format](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#haplegsample), 
and a genetic map can be obtained from the [Impute2 website](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#reference).

## Ligating haplotypes

It is recommended to ligate haplotypes using
[hapfuse](https://bitbucket.org/wkretzsch/hapfuse/src).  Before
fusing, the output from GLPhase needs to be converted from gzipped VCF
to something htslib can read. Here an example using [bcftools](https://samtools.github.io/bcftools/bcftools.html):

    zcat output_base_name.vcf.gz | bcftools -Ob -o \
        output_base_name.bcf

