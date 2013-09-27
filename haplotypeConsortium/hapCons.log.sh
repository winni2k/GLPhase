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

from GLs
BEAGLE -> posteriors + haplotypes
fix posteriors 
use haplotypes to initialize mcmc machinery

-> genotypes from Shapeit compare to CG data overlaps with 1kg


calibration
1. compare haplotypes in down-stream 
2. given threshold on posteriors, what amount of data is called and what is error rate.

CG data
two samples
- half sequenced in 1kg project, use for concordance measurements
- half not in 1kg, use for imputation experiments

EOF

###########
# Thu Sep 19 12:27:24 BST 2013
# got first set of bin files from chr20 on 7000 samples by running
ssh feng
cd /homes/kretzsch/feng/marchini/haplotypeConsortium/chr20_pilot/results/2013-09-17_run_insti_on_sanger_GLs/
./runall.pl -r -j 10

# then copied all bin files to well using git annex
# on well had to retouch prereqs

###########
# Thu Sep 19 13:51:43 BST 2013
# copied bin files to /well
rsync -navP ./regions/7213.427589/chr20/*.bin\
 dense:/well/marchini/winni/proj/marchini/haplotypeConsortium/chr20_pilot/results/2013-09-17_run_insti_on_sanger_GLs/regions/7213.427589/chr20/

# then ran command from cluster3
ssh cluster3
cd /well/marchini/winni/proj/marchini/haplotypeConsortium/chr20_pilot/results/2013-09-17_run_insti_on_sanger_GLs
./runall.pl -m /well/marchini/winni/proj/marchini -t -c 'localhost' 
./runall.pl -m /well/marchini/winni/proj/marchini -i -q 'short.qb' -P 'marchini.prjb'

# then copied concatenated vcf back to fenghuang for analysis

