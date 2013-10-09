###########
# Fri Jul 19 10:33:25 BST 2013
# this log contains all work that I have done to backup BAM files from /data/BGI on sparse to external HDs

# on sparse
BGI=/Net/sparse/data/BGI
bgiBackup.pl\
  bamList=BGI/
  
# creating valid test bam files
java -jar opt/picard-tools-1.78/DownsampleSam.jar I=/data/BGI/01-10-2012/HKC11064_HUMzwrR/data/MD_CHW_AAS_10179.bam  O=MD_CHW_AAS_10179.bam PROBABILITY=0.00005

FILE=MD_CHW_AAS_10179.downsamp.bam
FILE=MD_CHW_AAS_10011.downsamp.bam
FILE=MD_CHW_AAS_10011.invalid.bam
java -jar ~/opt/ValidateSamFile.jar I=$FILE O=$FILE.validation MODE=SUMMARY VALIDATE_INDEX=TRUE && echo hi world

# see ~/dev/bamTrackLib for further work on the backup suite

# now run on sparse
ssh sparse

# register bams
bamTracker.pl validateSamJar=~/opt/ValidateSamFile.jar addBamList=/data/BGI/bamlist.curr

# now back them up to externalHD1
# first test 17 problem bams
bamTracker.pl validateSamJar=~/opt/ValidateSamFile.jar backupBamList=doubleSlash.bam.list backupTargetDir=cairparavel:/mnt/extHD1/BGI backupDeviceName=externalhd1 doNotValidate=1 doNotIndex=1

# ok, now do this on first 1000 of bamlist.curr
head -n 1000 /data/BGI/bamlist.curr > bamlist.curr.head1000
bamTracker.pl validateSamJar=~/opt/ValidateSamFile.jar backupBamList=bamlist.curr.head1000 backupTargetDir=cairparavel:/mnt/extHD1/BGI backupDeviceName=externalhd1 doNotValidate=1 doNotIndex=1

# now backup all files that exist
bamTracker.pl validateSamJar=~/opt/ValidateSamFile.jar backupBamList=bamlist.curr.exists backupTargetDir=cairparavel:/mnt/extHD1/BGI backupDeviceName=converge4 doNotValidate=1 doNotIndex=1

###########
# Mon Oct 07 14:16:52 BST 2013
# create two lists
wc bamlist.curr.exists 
shuf bamlist.curr.exists > bamlist.curr.exists.shuf
head -n 4920 bamlist.curr.exists.shuf >  bamlist.curr.exists.shuf.A
tail -n 4920 bamlist.curr.exists.shuf >  bamlist.curr.exists.shuf.B

# send the first list to converge7
HD=converge7
bamTracker.pl validateSamJar=~/opt/ValidateSamFile.jar backupBamList=bamlist.curr.exists.shuf.A backupTargetDir=cairparavel:/mnt/$HD/BGI backupDeviceName=$HD doNotValidate=1 doNotIndex=1

HD=converge8
LIST=bamlist.curr.exists.shuf.B
bamTracker.pl validateSamJar=~/opt/ValidateSamFile.jar backupBamList=$LIST backupTargetDir=cairparavel:/mnt/$HD/BGI backupDeviceName=$HD doNotValidate=1 doNotIndex=1