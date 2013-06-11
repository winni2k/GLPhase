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
cd ~/fenghuang/marchini/wimpute/results/2013-04-29
zcat 16_034437067_034765685.bin.log.gz | grep -P '^\d+\t(0|4500|9299)\t' | gzip -c > 16_034437067_034765685.bin.log.samples0_4500_9299.tsv.gz

###########
# Sat May 04 21:19:17 BST 2013
# here a list of things I ran yesterday and today
cd ~/fenghuang/marchini/wimpute/results/
./2013-05-02/analyze.pl -r -j 5 -i # this is the long reference run
./2013-05-04/analyze.fold1.pl -r -j 5 # half as many cycles (9303 == sample size)

###########
# Sat May 04 21:21:33 BST 2013
# jfsag: running at 0.5 fold
./2013-05-04/analyze.fold0.5.pl -r -j 5 # quarter as many cycles (4500)






