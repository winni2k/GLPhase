###########
# Wed Apr 15 14:27:25 BST 2015
# created VCF with GLs from .bin file

# first the body
cd glReader
gzip -cd ../20_011976121_012173018.bin.onlyThree.bin | ~/marchini/scripts/stbin2vcf.pl -c20 > 20_011976121_012173018.bin.onlyThree.vcf
bgzip -c 20_011976121_012173018.bin.onlyThree.vcf > 20_011976121_012173018.bin.onlyThree.vcf.gz
bcftools index 20_011976121_012173018.bin.onlyThree.vcf.gz
bcftools view 20_011976121_012173018.bin.onlyThree.vcf.gz -Ou -o 20_011976121_012173018.bin.onlyThree.bcf
bcftools view 20_011976121_012173018.bin.onlyThree.vcf.gz -Ob -o 20_011976121_012173018.bin.onlyThree.bcf.gz

