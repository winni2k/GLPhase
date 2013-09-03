# first agenda item: report bug to broad
# this fails
/usr/local/bin/java -Xmx128G -jar ~/src/GenomeAnalysisTK-2.7-2-g6bda569/GenomeAnalysisTK.jar -T UnifiedGenotyper -I /data/itch/winni/proj/marchini/converge/ir_bqsr_evaluation/data.dense.wg/recalibrated_bams/MD_CHW_AAC_2699.recalibrated.bam -R /data/1kg/reference_v37d5/hs37d5.fa -o ./test.vcf

# create minimal bam
ssh dense
cd /home/winni/marchini/converge/variantCalling
samtools view -bh /data/itch/winni/proj/marchini/converge/ir_bqsr_evaluation/data.dense.wg/recalibrated_bams/MD_CHW_AAC_2699.recalibrated.bam 1:1-100000 > MD_CHW_AAC_2699.recalibrated.bam

# try minimal version
/usr/local/bin/java -Xmx128G -jar ~/src/GenomeAnalysisTK-2.7-2-g6bda569/GenomeAnalysisTK.jar -T UnifiedGenotyper -I /data/itch/winni/proj/marchini/converge/ir_bqsr_evaluation/data.dense.wg/recalibrated_bams/MD_CHW_AAC_2699.recalibrated.bam -R /data/1kg/reference_v37d5/hs37d5.fa -o ./test.vcf