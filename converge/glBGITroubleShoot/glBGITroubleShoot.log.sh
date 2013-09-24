###########
# Fri Sep 13 11:10:16 BST 2013
# trying to figure out if gls on the BGI side are ok

# pull out random site from chr20 in what I think is a good vcf
# got site from ~/rare_snps/lists/bgi_genotyping.chr20.sites.tsv
# 20:350490

# file exists and looks ok
vcftools  --gzvcf ~/rare_snps/shared/GLs/9300/chr20.prob.vcf.gz

# pull out site of interest 
tabix -h ~/rare_snps/shared/GLs/9300/chr20.prob.vcf.gz 20:350489-350490 > chr20.350490.9300.vcf

# ok files seem to be off
# what if I look at my poprob output?


