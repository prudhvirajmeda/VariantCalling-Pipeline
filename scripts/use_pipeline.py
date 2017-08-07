import sys
import os.path
import os
import subprocess

if len(sys.argv) < 3:
    print 'error'
    sys.exit()

inFile1 = sys.argv[1]
inFile2 = sys.argv[2]



outFile1 = inFile1+'.trimmed.paired'
outFile2 = inFile1+'.trimmed.unpaired'
outFile3 = inFile2+'.trimmed.paired'
outFile4 = inFile2+'.trimmed.unpaired'

# trimming the sequences
trim_cmd ='java -jar /home/parikhhi/.pyenv/versions/2.7.3/bin/trimmomatic-0.33.jar PE -phred64'+" "+inFile1+" "+inFile2+" "+outFile1+" "+ outFile2+" "+" "+outFile3+" "+ outFile4+" "+ 'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33'  
subprocess.call(trim_cmd,shell=True)

# trimming metrics
trim_metrtics_cmd='python'+" "+'~/count_pretrim.py'+" "+inFile1+" "+inFile2+" "+outFile1+" "+outFile3+" "+'trim_metrics.txt'
subprocess.call(trim_metrtics_cmd,shell=True)
#print trim_metrtics_cmd

#bwa alignment
samplename =' FCC1BRRACXX'
#outFile='aligned.sam'
bwa_cmd ='/usr/global/blp/bin/bwa mem -t 20 -M -R'+'\'@RG\\tID:'+samplename+'\\tLB:'+samplename+'\\tSM:'+samplename+'\\tPL:ILLUMINA\''+" "+'/gpfs_fs/data1/refdb/genomes/Homo_sapiens/UCSC/bwa_index/hg19.fa'+" "+ outFile1+" "+outFile3+" "+\
'>'+" "+'aligned.sam'
#print cmd
subprocess.call(bwa_cmd,shell=True);

#convertin sam to bam 
sam_2_bam_cmd = 'java -jar /usr/global/blp/picard-tools-1.95/SamFormatConverter.jar'+" "+"I="+'aligned.sam'+" "+"O="+'aligned.bam'
subprocess.call(sam_2_bam_cmd,shell=True);
#print sam_2_bam_cmd

#alignment metrics
#align_metrics_cmd='java -jar /usr/global/blp/picard-tools-1.95/CollectAlignmentSummaryMetrics.jar'+" "+"R=/gpfs_fs/data1/refdb/genomes/Homo_sapiens/UCSC/bwa_index/hg19.fa"+" "+"I="+'aligned.bam'+" "+"O=alignment_metrics.txt"
#subprocess.call(align_metrics_cmd,shell=True)
#print align_metrics_cmd


#sort co-ordinates
sort_cmd = 'java -jar /usr/global/blp/picard-tools-1.95/SortSam.jar'+" "+"I=aligned.bam"+" "+'O=sorted.bam'+" "+'SO=coordinate'
subprocess.call(sort_cmd,shell=True);
#print sort_cmd

#alignment metrics                                                                                                                                                     
align_metrics_cmd='java -jar /usr/global/blp/picard-tools-1.95/CollectAlignmentSummaryMetrics.jar'+" "+"R=/gpfs_fs/data1/refdb/genomes/Homo_sapiens/UCSC/bwa_index/hg19.fa"+" "+"I="+'sorted.bam'+" "+"O=alignment_metrics.txt"
subprocess.call(align_metrics_cmd,shell=True)                                                                                                                         
#print align_metrics_cmd




#re-order with the reference file
re_order_cmd = 'java -jar /usr/global/blp/picard-tools-1.95/ReorderSam.jar'+" "+"I=sorted.bam"+" "+'O=sorted_aligned_ref.bam'+" "+"R=/gpfs_fs/data1/refdb/genomes/Homo_sapiens/UCSC/reordered/hg19.reorder.fa"
subprocess.call(re_order_cmd,shell=True);
#print re_order_cmd

#mark duplicates
mark_cmd = 'java -jar /usr/global/blp/picard-tools-1.95/MarkDuplicates.jar'+" "+'I=sorted_aligned_ref.bam'+" "+"O=marked_duplicates_SorAligRef.bam"+" "+'METRICS_FILE=marked_dup_metrics_SortAlignRef.txt'
subprocess.call(mark_cmd,shell=True); 
#print mark_cmd

#bam index                                                                                                                                                                          
index_cmd = 'java -jar /usr/global/blp/picard-tools-1.95/BuildBamIndex.jar'+" "+'I=marked_duplicates_SorAligRef.bam'
subprocess.call(index_cmd,shell=True);
#print index_cmd  

#realigner target creater
realign_tar_cmd = '/usr/global/blp/jdk1.7.0_45/bin/java'+" "+'-jar'+" "+"/usr/global/blp/GenomeAnalysisTK-3.1.1/GenomeAnalysisTK.jar"+" "+'-T'+" "+"RealignerTargetCreator"+" "+"-R"+" "+'/gpfs_fs/data1/refdb/genomes/Homo_sapiens/UCSC/reordered/hg19.reorder.fa'+" "+'-known'+" "+'/gpfs_fs/bnfo620/exome_data/dbsnp_138.hg19.vcf'+" "+'-I'+" "+'marked_duplicates_SorAligRef.bam'+" "+'-o'+" "+'realignertargetcreator.intervals'
subprocess.call(realign_tar_cmd,shell=True);
#print realign_tar_cmd

#left indel
left_indel_cmd = '/usr/global/blp/jdk1.7.0_45/bin/java'+" "+'-jar'+" "+"/usr/global/blp/GenomeAnalysisTK-3.1.1/GenomeAnalysisTK.jar"+" "+"-R"+" "+'/gpfs_fs/data1/refdb/genomes/Homo_sapiens/UCSC/reordered/hg19.reorder.fa'+" "+'-T'+" "+'LeftAlignIndels'+" "+'-I'+" "+'marked_duplicates_SorAligRef.bam'+" "+'-o'+" "+'output_with_leftaligned_indels.bam'
subprocess.call(left_indel_cmd,shell=True);
#print left_indel_cmd

#indel realigner
indel_realigner_cmd='/usr/global/blp/jdk1.7.0_45/bin/java'+" "+'-jar'+" "+"/usr/global/blp/GenomeAnalysisTK-3.1.1/GenomeAnalysisTK.jar"+" "+'-T'+" "+'IndelRealigner'+" "+"-R"+" "+'/gpfs_fs/data1/refdb/genomes/Homo_sapiens/UCSC/reordered/hg19.reorder.fa'+" "+'-targetIntervals'+" "+'realignertargetcreator.intervals'+" "+'-known'+" "+'/gpfs_fs/bnfo620/exome_data/dbsnp_138.hg19.vcf'+" "+'-I'+" "+'output_with_leftaligned_indels.bam'+" "+'-o'+" "+'indel_reliagned.bam'
subprocess.call(indel_realigner_cmd,shell=True);
#print indel_realigner_cmd


#base recalibration
base_call_cmd = '/usr/global/blp/jdk1.7.0_45/bin/java'+" "+'-jar'+" "+"/usr/global/blp/GenomeAnalysisTK-3.1.1/GenomeAnalysisTK.jar"\
+" "+'-T'+" "+'BaseRecalibrator'+" "+"-R"+" "+'/gpfs_fs/data1/refdb/genomes/Homo_sapiens/UCSC/reordered/hg19.reorder.fa'+" "+'-I'+" "+'indel_reliagned.bam'+" "+'-knownSites'+" "+'/gpfs_fs/bnfo620/exome_data/dbsnp_138.hg19.vcf'+" "+'-knownSites'+" "+'/gpfs_fs/bnfo620/exome_data/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf'+" "+'-o'+" "+'recal_data.table'
subprocess.call(base_call_cmd,shell=True);
#print base_call_cmd

#print reads
                                                                       
print_reads_cmd = '/usr/global/blp/jdk1.7.0_45/bin/java'+" "+'-jar'+" "+"/usr/global/blp/GenomeAnalysisTK-3.1.1/GenomeAnalysisTK.jar"+" "+'-T'+" "+'PrintReads'+" "+"-R"+" "+'/gpfs_fs/data1/refdb/genomes/Homo_sapiens/UCSC/reordered/hg19.reorder.fa'+" "+'-I'+" "+'indel_reliagned.bam'+" "+'-BQSR'+" "+'recal_data.table'+" "+'-o'+" "+'recal_reads.bam'
subprocess.call(print_reads_cmd,shell=True);
#print print_reads_cmd

#vcf calling

vcf_cmd=  '/usr/global/blp/jdk1.7.0_45/bin/java'+" "+'-jar'+" "+"/usr/global/blp/GenomeAnalysisTK-3.1.1/GenomeAnalysisTK.jar"+" "+'-T'+" "+'HaplotypeCaller'+" "+"-R"+" "+'/gpfs_fs/data1/refdb/genomes/Homo_sapiens/UCSC/reordered/hg19.reorder.fa'+" "+'-I'+" "+'recal_reads.bam'+" "+'--genotyping_mode'+" "+'DISCOVERY'+" "+'-stand_emit_conf'+" "+'10'+" "+'-stand_call_conf'+" "+'30'+" "+'-o'+" "+'raw_variants.vcf'
subprocess.call(vcf_cmd,shell=True);
#print vcf_cmd

#vcu filetering
#saperate snp's form raw vcf file

sap_snp_cmd = '/usr/global/blp/jdk1.7.0_45/bin/java'+" "+'-jar'+" "+"/usr/global/blp/GenomeAnalysisTK-3.1.1/GenomeAnalysisTK.jar"+" "+'-T'+" "+'SelectVariants'+" "+"-R"+" "+'/gpfs_fs/data1/refdb/genomes/Homo_sapiens/UCSC/reordered/hg19.reorder.fa'+" "+'-V'+" "+'raw_variants.vcf'+" "+'-selectType'+" "+'SNP'+" "+'-o'+" "+'raw_snps.vcf'
subprocess.call(sap_snp_cmd,shell=True);
#print sap_snp_cmd

#apply fileter to snp call set

filter_snp_cmd = '/usr/global/blp/jdk1.7.0_45/bin/java'+" "+'-jar'+" "+"/usr/global/blp/GenomeAnalysisTK-3.1.1/GenomeAnalysisTK.jar"+" "+'-T'+" "+'VariantFiltration'+" "+"-R"+" "+'/gpfs_fs/data1/refdb/genomes/Homo_sapiens/UCSC/reordered/hg19.reorder.fa'+" "+'-V'+" "+'raw_snps.vcf'+" "+'--filterExpression'+" "+'"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"'+" "+'--filterName'+" "+'"my_snp_filter"'+" "+'-o'+" "+'filtered_snps.vcf'
subprocess.call(filter_snp_cmd,shell=True);
#print filter_snp_cmd 

#saparate indels from raw files

sap_indel_cmd ='/usr/global/blp/jdk1.7.0_45/bin/java'+" "+'-jar'+" "+"/usr/global/blp/GenomeAnalysisTK-3.1.1/GenomeAnalysisTK.jar"+" "+'-T'+" "+'SelectVariants'+" "+"-R"+" "+'/gpfs_fs/data1/refdb/genomes/Homo_sapiens/UCSC/reordered/hg19.reorder.fa'+" "+'-V'+" "+'raw_variants.vcf'+" "+'-selectType'+" "+'INDEL'+" "+'-o'+" "+'raw_indels.vcf'
subprocess.call(sap_indel_cmd,shell=True);
#print sap_indel_cmd

#apply filters to indel call set

filter_indel_cmd='/usr/global/blp/jdk1.7.0_45/bin/java'+" "+'-jar'+" "+"/usr/global/blp/GenomeAnalysisTK-3.1.1/GenomeAnalysisTK.jar"+" "+'-T'+" "+'VariantFiltration'+" "+"-R"+" "+'/gpfs_fs/data1/refdb/genomes/Homo_sapiens/UCSC/reordered/hg19.reorder.fa'+" "+'-V'+" "+'raw_indels.vcf'+" "+'--filterExpression'+" "+'"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"'+" "+'--filterName'+" "+'"my_indel_filter"'+" "+'-o'+" "+'filtered_indels.vcf'
subprocess.call(filter_indel_cmd,shell=True);
#print filter_indel_cmd


# run parser
run_parser_cmd='python'+" "+'~/test_parser.py'
subprocess.call(run_parser_cmd,shell=True);

#reference command to the above
#RealignerTargetCreator
#rtg_cmd= '
#/usr/global/blp/jdk1.7.0_45/bin/java -jar /usr/global/blp/GenomeAnalysisTK-3.1.1/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /gpfs_fs/data1/refdb/genomes/Homo_sapiens/UCSC/reordered/hg19.reorder.fa -known /gpfs_fs/bnfo620/exome_data/dbsnp_138.hg19.vcf -I sorted_aligned_ref.bam -o realignertargetcreator.intervals

#indel to left
#left_




#/usr/global/blp/jdk1.7.0_45/bin/java -jar /usr/global/blp/GenomeAnalysisTK-3.1.1/GenomeAnalysisTK.jar -T IndelRealigner -R /gpfs_fs/data1/refdb/genomes/Homo_sapiens/UCSC/reordered/hg19.reorder.fa -targetIntervals realignertargetcreator.intervals -known /gpfs_fs/bnfo620/exome_data/dbsnp_138.hg19.vcf -I sorted_aligned_ref.bam -o sorted_aligned_ref_indelrealigner.bam
