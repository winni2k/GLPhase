# first agenda item: report bug to broad
# this fails
/usr/local/bin/java -Xmx128G -jar ~/src/GenomeAnalysisTK-2.7-2-g6bda569/GenomeAnalysisTK.jar -T UnifiedGenotyper -I /data/itch/winni/proj/marchini/converge/ir_bqsr_evaluation/data.dense.wg/recalibrated_bams/MD_CHW_AAC_2699.recalibrated.bam -R /data/1kg/reference_v37d5/hs37d5.fa -o ./test.vcf

# create minimal bam
ssh dense
cd /home/winni/marchini/converge/variantCalling
samtools view -bh /data/itch/winni/proj/marchini/converge/ir_bqsr_evaluation/data.dense.wg/recalibrated_bams/MD_CHW_AAC_2699.recalibrated.bam 1:1-100000 > MD_CHW_AAC_2699.recalibrated.bam

# try minimal version
/usr/local/bin/java -Xmx128G -jar ~/src/GenomeAnalysisTK-2.7-2-g6bda569/GenomeAnalysisTK.jar -T UnifiedGenotyper -I ./MD_CHW_AAC_2699.recalibrated.bam -R /data/1kg/reference_v37d5/hs37d5.fa -o ./test.vcf

###########
# Tue Sep 03 18:54:27 BST 2013
### new bug, if index file is older gatk still tries accessing, even if there is no newer index available
/usr/local/bin/java -Xmx128G -jar /data/itch/caina/software/GenomeAnalysisTK-2.6-5-gba531bd/GenomeAnalysisTK.jar -T UnifiedGenotyper -I /data/itch/winni/proj/marchini/converge/ir_bqsr_evaluation/data.dense.wg/recalibrated_bams/MD_CHW_AAO_3974.recalibrated.bam -R /data/1kg/reference_v37d5/hs37d5.fa -o test.vcf

# some bai files are simply outdated
# regenerating
cd ~/ir_bqsr_evaluation/data.dense.wg/recalibrated_bams/

for i in $(find . -name "*.bam"); do j=$(basename $i .bam); test $j.bam.bai -ot $i && echo $PWD/$i >> ~/marchini/converge/variantCalling/bamsWithOldBai.bam.list; done

cd /home/winni/marchini/converge/variantCalling
parallel -a bamsWithOldBai.bam.list -j 2 ddIndex {}

###########
# Wed Sep 04 12:13:54 BST 2013
# new bug something wrong with this file:
/data/itch/winni/proj/marchini/converge/ir_bqsr_evaluation/data.dense.wg/recalibrated_bams/MD_CHW_AAB_2226.recalibrated.bam

/usr/local/bin/java -Xmx128G -jar /data/itch/caina/software/GenomeAnalysisTK-2.6-5-gba531bd/GenomeAnalysisTK.jar -T UnifiedGenotyper -I /data/itch/winni/proj/marchini/converge/ir_bqsr_evaluation/data.dense.wg/recalibrated_bams/MD_CHW_AAB_2226.recalibrated.bam -R /data/1kg/reference_v37d5/hs37d5.fa -o test.vcf