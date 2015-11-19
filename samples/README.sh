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

###########
# Tue Nov 17 09:39:53 GMT 2015
# need to create a test file for multiallelics
mkdir samples/multi_gls
cd samples/multi_gls

# create known genotypes
bt view -r20 /gpfs1/well/marchini/winni/proj/marchini/haplotypeConsortium/data/2014-10-14_1kg_omni_sites/ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz -Ou | bt view -T ../../2014-10-19_concordance_analysis/hrc_pilot_ISEC_hrc_r1.positions -Ou | bt view -t20:5335724-5861377 -Ou | bt view -S../../2014-10-19_concordance_analysis/omni_ISEC_hrc_pilot_ISEC_hrc_r1.sample.list -Oz -o 20.5335724-5861377.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz

# create list of multiallelic GLs
bt view -r20:5335724-5861377 ../../../../data/2014-08-27_sanger_release1_GLs/merged-likelihoods-20140904-winniFilt/chr20.HRC.r1.AC5.32920_samples.likelihoods.winniFilter.bcf.gz -G  -H | cut -f2 | uniq -d | head -n 256 > multi.sites

# subset out multi sites (52 total)
bt view -Rmulti.sites ../../../../data/2014-08-27_sanger_release1_GLs/merged-likelihoods-20140904-winniFilt/chr20.HRC.r1.AC5.32920_samples.likelihoods.winniFilter.bcf.gz -Ob -o chr20.5335724-5861377.only_multi.HRC.r1.AC5.32920_samples.likelihoods.winniFilter.bcf.gz

# create site list of genotyped sites
bt view -H 20.5335724-5861377.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz | cut -f1,2 > genotyped.sites

# extract genotyped sites
bt view -r20:5335724-5861377 ../../../../data/2014-08-27_sanger_release1_GLs/merged-likelihoods-20140904-winniFilt/chr20.HRC.r1.AC5.32920_samples.likelihoods.winniFilter.bcf.gz -Ou | bt view -T genotyped.sites -Ob -o chr20.5335724-5861377.genotyped.HRC.r1.AC5.32920_samples.likelihoods.winniFilter.bcf.gz

### now fill up remaining sites
cat multi.sites genotyped.sites | sort | uniq | cut -f2 > exclude.pos

# sample another 468 sites
bt view -r20:5335724-5861377 ../../../../data/2014-08-27_sanger_release1_GLs/merged-likelihoods-20140904-winniFilt/chr20.HRC.r1.AC5.32920_samples.likelihoods.winniFilter.bcf.gz -G -H | cut -f1,2 | grepLarge.pl -v -fexclude.pos -c2 | shuf | head -n 479 > extra.sites

# extract extra site GLs
bt view -r20:5335724-5861377 ../../../../data/2014-08-27_sanger_release1_GLs/merged-likelihoods-20140904-winniFilt/chr20.HRC.r1.AC5.32920_samples.likelihoods.winniFilter.bcf.gz -Ou | bt view -Textra.sites -Ob -o chr20.5335724-5861377.extra.HRC.r1.AC5.32920_samples.likelihoods.winniFilter.bcf.gz
bt index chr20.5335724-5861377.extra.HRC.r1.AC5.32920_samples.likelihoods.winniFilter.bcf.gz


# create MSTMv3 sample intersect file
comm -12 <(bt query -l ../../2014-09-01_GL_hap_intersect/unionMAC5.winni_filt_subset.with_multi/20.union.filteredAC5.onlyPhased.NM_HOMMAJORv3.inGLSamples.winni_filt_subset.with_multi.bcf.gz) ../../2014-09-01_GL_hap_intersect/GL_haplotype_intersect.noCG9.TGPP3_ISEC.samples > MSTMv3.TGPP3_ISEC.samples

# concatenate files
bt concat -D -a \
   chr20.5335724-5861377.only_multi.HRC.r1.AC5.32920_samples.likelihoods.winniFilter.bcf.gz\
   chr20.5335724-5861377.extra.HRC.r1.AC5.32920_samples.likelihoods.winniFilter.bcf.gz\
   chr20.5335724-5861377.genotyped.HRC.r1.AC5.32920_samples.likelihoods.winniFilter.bcf.gz\
   -Ou |\
    bt view -S MSTMv3.TGPP3_ISEC.samples\
       -Ob -o chr20.5335724-5861377.1024_site_subset.HRC.r1.AC5.TGPP3_samples.likelihoods.winniFilter.bcf.gz

bt index chr20.5335724-5861377.1024_site_subset.HRC.r1.AC5.TGPP3_samples.likelihoods.winniFilter.bcf.gz 
bt view -G -H  chr20.5335724-5861377.1024_site_subset.HRC.r1.AC5.TGPP3_samples.likelihoods.winniFilter.bcf.gz | wc -l


# get pre-existing haps
bt isec -w1 -n=2 ../../2014-09-01_GL_hap_intersect/unionMAC5.winni_filt_subset.with_multi/20.union.filteredAC5.onlyPhased.NM_HOMMAJORv3.inGLSamples.winni_filt_subset.with_multi.bcf.gz chr20.5335724-5861377.1024_site_subset.HRC.r1.AC5.TGPP3_samples.likelihoods.winniFilter.bcf.gz -Ou | bt view -S MSTMv3.TGPP3_ISEC.samples -Ob -o 20.5335724-5861377.1024_site_subset.union.filteredAC5.onlyPhased.NM_HOMMAJORv3.inGLSamples.winni_filt_subset.with_multi.bcf.gz

# create GL site list
bt view chr20.5335724-5861377.1024_site_subset.HRC.r1.AC5.TGPP3_samples.likelihoods.winniFilter.bcf.gz -G -Oz -o chr20.5335724-5861377.1024_site_subset.HRC.r1.AC5.TGPP3_samples.likelihoods.winniFilter.sites.vcf.gz

# order haps by site list
bt view 20.5335724-5861377.1024_site_subset.union.filteredAC5.onlyPhased.NM_HOMMAJORv3.inGLSamples.winni_filt_subset.with_multi.bcf.gz | ../../scripts/reorder_multi.pl -r chr20.5335724-5861377.1024_site_subset.HRC.r1.AC5.TGPP3_samples.likelihoods.winniFilter.sites.vcf.gz  | bt view -Ob -o 20.5335724-5861377.1024_site_subset.union.filteredAC5.onlyPhased.NM_HOMMAJORv3.inGLSamples.winni_filt_subset.with_multi.ordered.bcf.gz

# generate tabhaps for older versions
PHBASE=20.5335724-5861377.1024_site_subset.union.filteredAC5.onlyPhased.NM_HOMMAJORv3.inGLSamples.winni_filt_subset.with_multi.ordered
bt convert --hapsample $PHBASE $PHBASE.bcf.gz

zcat $PHBASE.hap.gz | tr ' ' '\t' | bgzip -c> $PHBASE.tabhaps.gz

# tabix tabhaps
tabix -s 1 -b 3 -e 3 $PHBASE.tabhaps.gz 


###########
# Tue Nov 17 14:57:24 GMT 2015
# remove multiallelics
bt view chr20.5335724-5861377.1024_site_subset.HRC.r1.AC5.TGPP3_samples.likelihoods.winniFilter.bcf.gz -G -H |cut -f1,2 | sort |uniq -d > duplicate.sites

bt view -T^duplicate.sites chr20.5335724-5861377.1024_site_subset.HRC.r1.AC5.TGPP3_samples.likelihoods.winniFilter.bcf.gz -Ob -o chr20.5335724-5861377.1024_site_subset.HRC.r1.AC5.TGPP3_samples.likelihoods.winniFilter.no_multi.bcf.gz

###########
# Tue Nov 17 15:59:34 GMT 2015
# let's find real error
/data/fenghuang/not-backed-up/kretzsch/dev/insti/t/../bin/insti.1.4.13b.ncuda -g /data/fenghuang/not-backed-up/kretzsch/dev/insti/t/../samples/geneticMap/genetic_map_chr20_combined_b37.txt.gz -C100 -m 10 -B0 -i5 -h /data/fenghuang/not-backed-up/kretzsch/dev/insti/t/../samples/multi_gls/20.5335724-5861377.1024_site_subset.union.filteredAC5.onlyPhased.NM_HOMMAJORv3.inGLSamples.winni_filt_subset.with_multi.ordered.tabhaps.gz -s /data/fenghuang/not-backed-up/kretzsch/dev/insti/t/../samples/multi_gls/20.5335724-5861377.1024_site_subset.union.filteredAC5.onlyPhased.NM_HOMMAJORv3.inGLSamples.winni_filt_subset.with_multi.ordered.sample gls.bin

zcat gls.bin.vcf.gz | bgzip -c> tmp && mv -f tmp gls.bin.vcf.gz
bt index gls.bin.vcf.gz

OSAMP=../samples/multi_gls/20.5335724-5861377.chip.omni_broad_sanger_combined.20140818.snps.genotypes.samples
OGT=../samples/multi_gls/20.5335724-5861377.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz
bt stats $OGT gls.bin.vcf.gz -S $OSAMP > gls.bin.vcf.gz.stats
grep NRD gls.bin.vcf.gz.stats

# ok, it looks like the expected error is below 15 %

###########
# Wed Nov 18 09:26:19 GMT 2015
# Let's see if this error is just due to faulty bcf reading
bt view chr20.5335724-5861377.1024_site_subset.HRC.r1.AC5.TGPP3_samples.likelihoods.winniFilter.bcf.gz -Oz -o chr20.5335724-5861377.1024_site_subset.HRC.r1.AC5.TGPP3_samples.likelihoods.winniFilter.vcf.gz

vcf2STbin.pl chr20.5335724-5861377.1024_site_subset.HRC.r1.AC5.TGPP3_samples.likelihoods.winniFilter.vcf.gz

# yes it is!

###########
# Wed Nov 18 10:32:52 GMT 2015
# convert known bin to VCF so we can try it
zcat 20_011976121_012173018.bin.onlyThree.bin | perl ../scripts/bin2vcf.pl | bgzip -c > 20_011976121_012173018.bin.onlyThree.vcf.gz
bt index 20_011976121_012173018.bin.onlyThree.vcf.gz
bt view 20_011976121_012173018.bin.onlyThree.vcf.gz -Ob -o 20_011976121_012173018.bin.onlyThree.bcf.gz
bt index 20_011976121_012173018.bin.onlyThree.bcf.gz
