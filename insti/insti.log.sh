### log of work I've done on insti

#### copied over from rare_snps.log.sh

###########
# Mon Apr 29 16:55:34 BST 2013
# created a list of 5 random .bin files
cd dense:~/rare_snps/asian_ref_panel
ls probin/ |grep "^16" |shuf | head -n 5 > ../lists/asian_ref_panel.probin.rand5.bin.list

#### end of copy over

###########
# Wed May 01 13:00:39 BST 2013
# grab sample 0 and sample 4500 from log files
ssh fenghuang
cd ~/fenghuang/marchini/insti/results/2013-04-29
zcat 16_034437067_034765685.bin.log.gz | grep -P '^\d+\t(0|4500|9299)\t' | gzip -c > 16_034437067_034765685.bin.log.samples0_4500_9299.tsv.gz

###########
# Sat May 04 21:19:17 BST 2013
# here a list of things I ran yesterday and today
cd ~/fenghuang/marchini/insti/results/
./2013-05-02/analyze.pl -r -j 5 -i # this is the long reference run
./2013-05-04/analyze.fold1.pl -r -j 5 # half as many cycles (9303 == sample size)

###########
# Sat May 04 21:21:33 BST 2013
# jfsag: running at 0.5 fold
./2013-05-04/analyze.fold0.5.pl -r -j 5 # quarter as many cycles (4500)

###########
# Sat Jun 01 16:52:48 BST 2013
# created a dir in which to summarize all the previous runs into one
mkdir results/2013-06-01_create-int-plot

# started work on integrated plot
touch results/2013-06-01_create-int-plot/create-int-plot.pl

# renamed some dirs in results to easier understand what the experiments were


###########
# Tue Jun 11 14:45:17 BST 2013
## implementing reference panel for insti
# create test data for gls
perl -ane '$o = q//; $o .= join("\t", (@F[0..2]));if(++$i==1){$o .= "\t".join("\t", qw/samp1 samp2 samp3/)} else{$o .= "\t". join("\t",( join(" ",(0,1)), join(" ", (1,0)), join(" ",(0,0))))} print $o ."\n"'< 20_011976121_012173018.bin.onlyThree.gls | bgzip -c > 20_011976121_012173018.bin.onlyThree.bin

# now create test data for legend
perl -ane 'if(++$i==1){print "rsid position a0 a1\n"}else{print join(" ",("$F[0]:$F[1]",$F[1], split("",$F[2]))). "\n"}' 20_011976121_012173018.bin.onlyThree.gls > 20_011976121_012173018.bin.onlyThree.legend

# create haplotypes file
perl -e 'for (0..600){print join(" ",qw/0 1 0 1/)."\n"}; for (601..1023){ print join(" ", qw/1 0 0 1/)."\n"}' > 20_011976121_012173018.bin.onlyThree.haps

###########
# Thu Aug 15 16:43:46 BST 2013
# let's work on the kickstart option

# that worked quite nicely

###########
# Wed Aug 21 11:52:35 BST 2013
# let's do some memory profiling
ssh feng
mkdir -p /homes/kretzsch/feng/marchini/insti/results/2013-08-21_Massif_memcheck
cd /homes/kretzsch/feng/marchini/insti/results/2013-08-21_Massif_memcheck
mkdir gls
cp ../2013-08-16_impute_C100_asianRefPanel/gls/20_059875288_060024976.bin.plusThree.bin gls/

# run massif memory check on insti in kickstart mode
valgrind --tool=massif ../../src/insti -m 10 -b 10 -C 1 -k -L ../2013-07-18b_insti_amh_100_cycles_kickstart/refPanelDir/ALL.chr20.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.only_snps.asians.polymorphic.20_059875288_060024976.bin.minusThree.legend.gz -H ../2013-07-18b_insti_amh_100_cycles_kickstart/refPanelDir/ALL.chr20.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.only_snps.asians.polymorphic.20_059875288_060024976.bin.minusThree.hap.gz -e 20_059875288_060024976.bin.IMPUTE_1000C_ASIAN_REF_KS.log.gz gls/20_059875288_060024976.bin.plusThree.bin


###########
# Thu Dec 05 13:24:43 GMT 2013
# create some extra test data 

# adding ref haps with extra sites
cd ~/marchini/insti/samples
perl -e 'BEGIN{@G = qw/A T G C/; print "id position a0 a1\n"}; for (1..1000){$pos = int(rand(62000001)); $out = "20:$pos $pos"; $g1 = $G[int(rand(@G))]; $g2 = $g1; while($g1 eq $g2){$g2=$G[int(rand(@G))]}; $out .= " $g1 $g2\n"; print $out;}' > 20_0_62000000.bin.onlyThree.legend

perl -e 'for(1..1000){ @a=(); for(1..4){ push @a, int(rand(2))}; print join(" ", @a)."\n";}' > 20_0_62000000.bin.onlyThree.haps

# combining snps we want with ones we don't want
paste -d"\n" 20_011976121_012173018.bin.onlyThree.legend 20_0_62000000.bin.onlyThree.legend |tail -n +2 | grep -P '.' > 20_0_62000000.011976121_012173018.paste.onlyThree.legend

paste -d"\n" 20_011976121_012173018.bin.onlyThree.hap 20_0_62000000.bin.onlyThree.hap | grep -P '.' > 20_0_62000000.011976121_012173018.paste.onlyThree.hap

# create haps file as well
~/marchini/scripts/hapLeg2haps.pl -c 20 -h 20_011976121_012173018.bin.onlyThree.hap -l 20_011976121_012173018.bin.onlyThree.legend > 20_011976121_012173018.bin.onlyThree.haps
~/marchini/scripts/hapLeg2haps.pl -c 20 -h 20_0_62000000.011976121_012173018.paste.onlyThree.hap -l 20_0_62000000.011976121_012173018.paste.onlyThree.legend > 20_0_62000000.011976121_012173018.paste.onlyThree.haps

# create sample files
echo -e "ID_1 ID_2 missing\n0 0 0\nsamp1 samp1 0\nsamp2 samp2 0\nsamp3 samp3 0\n" > onlyThree.hapsSample.extraLine.sample
echo -e "sample population group sex\nsamp1 CEU EUR 1\nsamp2 CEU EUR 2\nsamp3 CEU EUR 1\n" > onlyThree.hapLegSamp.extraLine.sample

echo -e "ID_1 ID_2 missing\n0 0 0\nsamp1 samp1 0\nsamp2 samp2 0\nsamp3 samp3 0" > onlyThree.hapsSample.sample
echo -e "sample population group sex\nsamp1 CEU EUR 1\nsamp2 CEU EUR 2\nsamp3 CEU EUR 1" > onlyThree.hapLegSamp.sample

# create scaffold haps
perl -ne 'chomp; print $_." 0 0\n"' 20_011976121_012173018.bin.onlyThree.haps| gshuf |head -50 > 20_011976121_012173018.bin.onlyThree.scaffold50.haps
~/marchini/scripts/haps2hapLegend.pl -i 20_011976121_012173018.bin.onlyThree.scaffold50.haps -o 20_011976121_012173018.bin.onlyThree.scaffold50
gunzip 20_011976121_012173018.bin.onlyThree.scaffold50.hap.gz

