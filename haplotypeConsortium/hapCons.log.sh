# This log details work for the haplotype consortium

# GLs are located at
/data/mandarin/backed-up/haplotype_consortium/

# start with chr20 pilot data
# have made symlinks in 
/homes/kretzsch/feng/marchini/haplotypeConsortium/chr20_pilot/GLs/sanger

# the sanger data contains 428131 GLs from 7213 individuals.
# all of them are SNPS:
echo <<EOF
VCFtools - v0.1.11
(C) Adam Auton 2009

Parameters as interpreted:
        --gzvcf SC-7312.chr20.unionAC10.PLs.vcf.gz
        --remove-indels

Using zlib version: 1.2.7
Reading Index file.
File contains 428131 entries and 7213 individuals.
Applying Required Filters.
Filtering sites by allele type
After filtering, kept 7213 out of 7213 Individuals
After filtering, kept 428131 out of a possible 428131 Sites
Run Time = 220.00 seconds
@fenghuang sanger$ 
EOF

# starting work on phasing pipeline
cd /homes/kretzsch/feng/marchini/haplotypeConsortium/chr20_pilot/results/2013-09-17_run_insti_on_sanger_GLs

###########
# Wed Sep 18 10:21:39 BST 2013
echo <<EOF

1. assumption all missing data is homozygous reference allele
2. SHAPEIT2 with missing data
3. run on intersect

Measure imputation performance of haplotypes

4. Impute from phase 3 haplotypes as baseline
Impute from haplotype consortium

+ GL data

take a look at all sites
alignment issues: 

1 remove all sites that arennt in any reference panel

 tools:
 gtool
 impute2 for imputation

measure imputation accuracy:
   Brian Howie's script

EOF