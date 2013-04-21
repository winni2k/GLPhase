SNPTools	1.0
pileline for NGS SNP analysis
author	Yi Wang @ Fuli Yu' Group
Baylor College of Medicine Human Genome Sequencing Center.
All rights reserved.

Workflow summary:
pileup->varisite->bamodel->poprob->probin->impute->hapfuse
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bamodel
three-component binomial mixture modeling of BAMs on the SAME individual
author	Yi Wang @ Fuli Yu' Group @ BCM-HGSC
usage	bamodel [options] <out site.vcf in1.bam> [in2.bam in3.bam...]
	-a <FLOAT>	P(read=ref|geno=alt/alt) (0.010)
	-h <FLOAT>	P(read=ref|geno=ref/alt) (0.500)
	-r <FLOAT>	P(read=ref|geno=ref/ref) (0.995)
	-p <FLOAT>	precision of beta distribution prior (100)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bamodel
three-component binomial mixture modeling of BAMs on one SAME individual
author	Yi Wang @ Fuli Yu' Group @ BCM-HGSC
usage	bamodel [options] <out site.vcf in1.bam> [in2.bam in3.bam...]
	-a <FLOAT>	P(read=ref|geno=alt/alt) (0.010)
	-h <FLOAT>	P(read=ref|geno=ref/alt) (0.500)
	-r <FLOAT>	P(read=ref|geno=ref/ref) (0.995)
	-p <FLOAT>	precision of beta distribution prior (100)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

hapfuse
joint chunked haplotypes into chromosome wide haplotypes
author	Yi Wang @ Fuli Yu' Group @ BCM-HGSC
usage	hapfuse <out.vcf dir> [gender]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

impute
haplotype imputation by cFDSL distribution
author	Yi Wang @ Fuli Yu' Group @ BCM-HGSC
usage	impute [options] 1.bin 2.bin ...
	-d <density>	relative SNP density to Sanger sequencing (1)
	-b <burn>	burn-in generations (56)
	-l <file>	list of input files
	-m <mcmc>	sampling generations (200)
	-n <fold>	sample size*fold of nested MH sampler iteration (2)
	-t <thread>	number of threads (0=MAX)
	-v <vcf>	integrate known genotype in VCF format
	-c <conf>	confidence of known genotype (0.9998)
	-x <gender>	impute x chromosome data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pileup
extract read depth from BAM
author	Yi Wang @ Fuli Yu' Group @ BCM-HGSC
usage	pileup [1.bam 2.bam ...]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

poprob
merge and transpose individual likelihood files to a population likelihood file
author	Yi Wang @ Fuli Yu' Group @ BCM-HGSC
usage	poprob [options] <site.vcf raws.list out.prob>
	-b <INT>	buffer size in MB (1024)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

prob2vcf
convert binary .prob file to tabixed .vcf.gz file
author	Yi Wang @ Fuli Yu' Group @ BCM-HGSC
usage	prob2vcf <in.prob out.vcf.gz chr>
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

probin
generate chunked .bin files for each chromosome
author	Yi Wang @ Fuli Yu' Group @ BCM-HGSC
usage	probin [options] <in.prob chr>
	-b <INT>	number of SNPs in a bin (1024)
	-s <INT>	stride of binning (512)
	-f <STR>	folder (./)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

varisite
allele selection and site scoring of variants
author	Yi Wang @ Fuli Yu' Group @ BCM-HGSC
usage	varisite [options] <ebd.list chr chr.fa>
	-i	consider indel allele (false)
	-s <FLT>	significance cutoff (1.5)
	-v <site.vcf>	only scoring on give sites
	-g <bit_mask>	group files by fields specified with bit_mask (0)
	-l <FLT>	variance stabilizer (1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Installation
The software is based on samtools's API library: libbam.a
The software is based on tabix's API library: libtabix.a
The software is based on GNU Scientific library
User should make sure that they are installed before make SNPTools
Edit makefile if they are not installed in default directories
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
License
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the 'Software'), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
