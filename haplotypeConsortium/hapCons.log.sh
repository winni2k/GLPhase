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