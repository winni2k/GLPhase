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
rsync -navP ./regions/7213.427589.519/chr20/*.bin  dense:/well/marchini/winni/proj/marchini/haplotypeConsortium/results/chr20_pilot/2013-09-17_run_insti_on_sanger_GLs/regions/7213.427589.519/chr20/

# then ran command from cluster3
ssh cluster3
cd /well/marchini/winni/proj/marchini/haplotypeConsortium/results/chr20_pilot/2013-09-17_run_insti_on_sanger_GLs
./runall.pl -m /well/marchini/winni/proj/marchini -t -c 'localhost' 
./runall.pl -m /well/marchini/winni/proj/marchini -i -q 'short.qb' -P 'marchini.prjb'

# then copied concatenated vcf back to fenghuang for analysis
ssh feng
cd /homes/kretzsch/feng/marchini/haplotypeConsortium/chr20_pilot/results/2013-09-17_run_insti_on_sanger_GLs/
rsync -avP cluster3:/well/marchini/winni/proj/marchini/haplotypeConsortium/results/chr20_pilot/2013-09-17_run_insti_on_sanger_GLs/merged/7213.427589.519/chr20.concat.vcf.gz ./merged/7213.427589.519/

###########
# Sun Sep 29 20:24:17 BST 2013
# did the same for completely merged GLs

rsync -navP ./regions/14513.427589.519/chr20/*.bin  dense:/well/marchini/winni/proj/marchini/haplotypeConsortium/results/chr20_pilot/2013-09-30_run_insti_on_combined_GLs/regions/14513.427589.519/chr20/

cd /well/marchini/winni/proj/marchini/haplotypeConsortium/results/chr20_pilot/2013-09-30_run_insti_on_combined_GLs
./runall.pl -m /well/marchini/winni/proj/marchini -t -c 'localhost' 
./runall.pl -m /well/marchini/winni/proj/marchini -i -q 'short.qb' -P 'marchini.prjb'

# then copied concatenated vcf back to fenghuang for analysis
ssh feng
cd /homes/kretzsch/feng/marchini/haplotypeConsortium/results/chr20_pilot/2013-09-30_run_insti_on_combined_GLs/
rsync -navP cluster3:/well/marchini/winni/proj/marchini/haplotypeConsortium/results/chr20_pilot/2013-09-30_run_insti_on_combined_GLs/merged/14513.427589.519/chr20.concat.vcf.gz ./merged/14513.427589.519/

# touch up the missing files
./runall.pl -m /well/marchini/winni/proj/marchini -i -t

# run shapeit
./runall.pl -m /well/marchini/winni/proj/marchini -i -s

###########
# Wed Oct 02 14:29:35 BST 2013
# started implementing the down-stream imputation 
# downloaded Omni1MQuad-v1.0-H to the data dir
# see the README.sh in feng:/homes/kretzsch/feng/marchini/haplotypeConsortium/data/2013-10-02_omni_chip_data/README.sh


###########
# Thu Oct 03 15:25:11 BST 2013
# results from 14500 samples were still disappointing
# reverted to pre GoCD insti version 1.0.11 and am running with -C1000
ssh cluster3
cd /well/marchini/winni/proj/marchini/haplotypeConsortium/results/chr20_pilot/2013-09-30_run_insti_on_combined_GLs
rm gprobs/14513.427589.519/chr20/*.vcf.gz
./runall.pl -m /well/marchini/winni/proj/marchini -i -q 'short.qb' -P 'marchini.prjb'

ssh feng
cd /homes/kretzsch/feng/marchini/haplotypeConsortium/results/chr20_pilot/2013-09-30_run_insti_on_combined_GLs/
rsync -navP cluster3:/well/marchini/winni/proj/marchini/haplotypeConsortium/results/chr20_pilot/2013-09-30_run_insti_on_combined_GLs/merged/14513.427589.519/chr20.concat.vcf.gz ./merged/14513.427589.519/

###########
# Sun Oct 06 19:58:57 BST 2013
#kicked off another set of jobs with C5000
ssh cluster3
cd /well/marchini/winni/proj/marchini/haplotypeConsortium/results/chr20_pilot/2013-09-30_run_insti_on_combined_GLs
./runall.pl -m /well/marchini/winni/proj/marchini -i -q 'long.qb' -P 'marchini.prjb' -C 5000

###########
# Mon Oct 07 10:30:11 BST 2013
# need to split concat vcf into chunks to run shapeit
cd /well/marchini/winni/proj/marchini/haplotypeConsortium/results/chr20_pilot/2013-09-30_run_insti_on_combined_GLs
/well/marchini/winni/proj/marchini/converge/rare_snps/scripts/BED_maker.pl -c 20 -s chr20.biAllelic.sites  -L 12000 -O 1000 -i > lists/chr20.shapeit.overlap.bed
/well/marchini/winni/proj/marchini/converge/rare_snps/scripts/BED_deoverlapper.pl < lists/chr20.shapeit.overlap.bed > lists/chr20.shapeit.deOverlap.bed

# move concat files to new 'C1000' position
cd merged/14513.427589.519/
mv chr20.concat.vcf.gz chr20.concat.C1000.vcf.gz
mv chr20.concat.vcf.gz.tbi chr20.concat.C1000.vcf.gz.tbi

./runall.pl -m /well/marchini/winni/proj/marchini -i -s -C 1000 -q 'short.qb' -P 'marchini.prjb'

###########
# Mon Oct 07 15:09:56 BST 2013
# no more phasing
# need to extract the allele frequencies from the original GLs file
# is in README of data dir 

###########
# Mon Oct 07 17:24:33 BST 2013
# try running snptools with 2.5k on amd cluster
cd /well/marchini/winni/proj/marchini/haplotypeConsortium/results/chr20_pilot/2013-09-30_run_insti_on_combined_GLs
./runall.pl -m /well/marchini/winni/proj/marchini -i -C 2500 -q 'short.qa' -P 'marchini.prja'

###########
# Tue Oct 08 17:00:35 BST 2013
# killed all jobs
# ran insti with 2000 cycles on b nodes
./runall.pl -m /well/marchini/winni/proj/marchini -i -C 2000 -q 'short.qb' -P 'marchini.prjb'

# ran beagle on same data
# but on feng first
ssh feng
cd /homes/kretzsch/feng/marchini/haplotypeConsortium/results/chr20_pilot/2013-09-30_run_insti_on_combined_GLs
./runall.pl -m ~/feng/marchini -b -q 'short.qb' -P 'mott-flint.prjb'

rsync -navP ./regions/14513.427589.519/chr20/*.BEAGLE.PL.gz  cluster3:/well/marchini/winni/proj/marchini/haplotypeConsortium/results/chr20_pilot/2013-09-30_run_insti_on_combined_GLs/regions/14513.427589.519/chr20/

###########
# Tue Oct 08 21:36:37 BST 2013
# now run beagle on cluster3
ssh cluster3
cd /well/marchini/winni/proj/marchini/haplotypeConsortium/results/chr20_pilot/2013-09-30_run_insti_on_combined_GLs
./runall.pl -m /well/marchini/winni/proj/marchini -b -q 'short.qb' -P 'mott-flint.prjb' -t
./runall.pl -m /well/marchini/winni/proj/marchini -b -q 'short.qb' -P 'mott-flint.prjb' -i


###########
# Wed Oct 09 18:00:01 BST 2013
# timing instiv1-0-11
ssh feng
cd /homes/kretzsch/feng/marchini/haplotypeConsortium/results/chr20_pilot/2013-09-30_run_insti_on_combined_GLs
./runall.pl -m ~/feng/marchini -C 2000 -I -t
./runall.pl -m ~/feng/marchini -C 2000 -I -i

cd timing
for i in 500 1000 2000; do echo $i; time /homes/kretzsch/feng/marchini/insti/src/insti -C $i chr20_62764710_62870647.C$i.bin > $i.log 2>&1 & done


###########
# Fri Oct 18 18:55:12 BST 2013
# taking a peak at beagle data
ssh cluster3
cd /users/winni/winni_on_marchini/proj/marchini/haplotypeConsortium/results/chr20_pilot/2013-09-30_run_insti_on_combined_GLs
qsub -sync y -cwd -V -b yes -j y -o distributedmake.log -N gprobs_merger -P marchini.prjb -r no -q short.qb /well/marchini/winni/proj/marchini/converge/rare_snps/scripts/gprobs_merger.pl -b chr20.noOverlap.dose.fileNames.bed -o merged/14513.427589.519/chr20.concat.C5000.peek.dose.gz -s

ssh feng
cd /homes/kretzsch/feng/marchini/haplotypeConsortium/results/chr20_pilot/2013-09-30_run_insti_on_combined_GLs
rsync -navP cluster3:/users/winni/winni_on_marchini/proj/marchini/haplotypeConsortium/results/chr20_pilot/2013-09-30_run_insti_on_combined_GLs/merged/14513.427589.519/chr20.concat.C5000.peek.*.gz merged/14513.427589.519/
