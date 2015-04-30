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


###########
# Tue Apr 21 16:30:43 BST 2015
# created simulated gls
cd hapGen
Rscript generateGLs.R

# convert expected haplotypes to vcf
(echo "sample population group sex"; zcat ex.bin |head -1 | tr '\t' '\n' |tail -n +4 | sed 's/$/ CEU EUR 1/') > ex.sample
perl -ane '$F[0]="20:$F[1]_$F[2]_$F[3]"; print join(" ", @F)."\n"' <ex.leg > ex.bt.leg
paste <(tail -n +2 ex.bt.leg) ex.haps | sed 's/\t/ /' | perl -ane 'print join(" ", ($F[0],@F))."\n"'| head -n 1000 > ex.wtccc.haps

(echo -e "ID_1 ID_2 missing\n0    0    0"; zcat ex.bin |head -1 | tr '\t' '\n' |tail -n +4 | perl -ne 'chomp; print "$_ $_ 0\n"') > ex.wtccc.sample

bcftools convert --hapsample2vcf ex.wtccc.haps,ex.wtccc.sample -Ov -oex.vcf

# need input gls in vcf format for beagle
zcat ex.gen.gz |perl -ane '$F[1] = "$F[0]:$F[2]_$F[3]_$F[4]"; print join(" ", @F)."\n"'| gzip -c > ex.forBT.gen.gz
bcftools convert --gensample2vcf ex.forBT.gen.gz,ex.wtccc.sample | bcftools +tag2tag -Oz -o ex.gls.vcf.gz -- -r --gp-to-gl

# and call genotypes with beagle
java -jar ~/opt/beagle4/b4.r1274.jar gl=ex.gls.vcf.gz out=ex.bgl

###########
# Wed Apr 29 14:30:30 BST 2015
# create pseudo multiallelic site
# edited simple.gls.v1.txt to have one quad allelic site
gzip -c simple.gls.v1.samePos.txt > simple.gls.v1.samePos.bin

###########
# Thu Apr 30 13:48:17 BST 2015
# added another allele that should cause insti to throw error
gzip -c simple.gls.v2.samePos.err.txt >simple.gls.v2.samePos.err.bin
