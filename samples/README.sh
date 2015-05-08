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

###########
# Fri May 01 11:17:37 BST 2015
# create a standard bin file with multiallelic sites
cp 20_011976121_012173018.bin.onlyThree.bin.txt  20_011976121_012173018.bin.onlyThree.withMultiall.bin.txt
# edit 20_011976121_012173018.bin.onlyThree.withMultiall.bin.txt 
gzip -c 20_011976121_012173018.bin.onlyThree.withMultiall.bin.txt  > 20_011976121_012173018.bin.onlyThree.withMultiall.bin

###########
# Fri May 01 18:31:55 BST 2015
# create simulated data with multiallelics
zcat ex.gls.vcf.gz | grep -v '^#' | awk 'BEGIN{OFS="\t"}{$4="C"; if(NR %2){$5="T"} else{$5 = "A"} $2=99999 + NR + NR%2; print}'  > tmp
(zcat ex.gls.vcf.gz | grep '^#'; cat tmp) | bgzip -c > ex.gls.multi.vcf.gz

# let's try calling genotypes with beagle
java -jar ~/opt/beagle4/beagle.r1399.jar gl=ex.gls.multi.vcf.gz out=ex.multi.bgl

###########
# Fri May 08 10:37:09 BST 2015
# need to create a multi version of ex.vcf.gz
zcat ex.vcf.gz | grep -v '^#' | awk 'BEGIN{OFS="\t"}{$4="C"; if(NR %2){$5="T"} else{$5 = "A"} $2=99999 + NR + NR%2; print}'  > tmp
(zcat ex.vcf.gz | grep '^#'; cat tmp) | bgzip -c > ex.multi.vcf.gz

# need to create multi version of ex.bin
zcat ex.bin | tail -n +2 | awk 'BEGIN{FS="\t"; OFS="\t"}{if(NR %2){$3="CT"} else{$3 = "CA"} $2=99999 + NR + NR%2; print}'  > tmp
(zcat ex.bin | head -1; cat tmp) | gzip -c > ex.multi.bin

# debugging
paste <(zcat ../../../samples/hapGen/ex.multi.vcf.gz | grep -v '^#' | cut -f1,2,4,5,10 )\
      <(zcat simulated_gls.gls.multi.vcf.gz.vcf.gz | grep -v '^#' | cut -f1,2,4,5,10 | perl -pne 's/([01]|[01]):\S+/$1/' )

###########
# Fri May 08 14:54:57 BST 2015
# let's see if the model simply performs poorly on closely spaced sites
# create bin file with all sites next to each other
zcat ex.gls.vcf.gz |  grep -v '^#' | awk 'BEGIN{FS="\t"; OFS="\t"}{if(NR %2){$3="CT"} else{$3 = "CA"} $2=99999 + NR; print}'  > tmp
(zcat ex.gls.vcf.gz | grep '^#'; cat tmp) | bgzip -c > ex.gls.almost_multi.vcf.gz

zcat ex.vcf.gz |  grep -v '^#' | awk 'BEGIN{FS="\t"; OFS="\t"}{if(NR %2){$3="CT"} else{$3 = "CA"} $2=99999 + NR; print}'  > tmp
(zcat ex.vcf.gz | grep '^#'; cat tmp) | bgzip -c > ex.almost_multi.vcf.gz

