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
bamTracker.pl validateSamJar=~/opt/ValidateSamFile.jar backupBamList=doubleSlash.bam.list backupTargetDir=cairparavel:/mnt/externalhd1/BGI backupDeviceName=externalhd1