#import sys

#if len(sys.argv)<2:
#    print 'error'
#    sys.exit()

f1 = open('trim_metrics.txt')
f2 = open('alignment_metrics.txt')
f3 = open('marked_dup_metrics_SortAlignRef.txt')
f4 = open('raw_variants.vcf')
f5 = open('raw_snps.vcf')
f6 = open('raw_indels.vcf')
f7 = open('filtered_snps.vcf')
f8 = open('filtered_indels.vcf')
ofh = open('stat.txt','w')

capt = []
for line in f1:
    line = line.strip()
    capt.append(line)

total_reads = capt[0]
trimmed_reads = capt[1]

align_cap1=[]
align_cap2=[]
align_cap3=[]
for line in f2:
    line = line.strip()
    if line.startswith('#'):
        f2.next()
    elif line.startswith(" "):
        f2.next()
    
    else:
        line = line.strip()
        line = line.split('\t')
        if len(line[0])==0:
#            print len(line)
            f2.next()
#            print 'hi'
        else:
            row2 = line[1]       #total reads
            row3 =line[5]        #aligned reads
            row4= line[6]
            align_cap1.append(row2)
            align_cap2.append(row3)
            align_cap3.append(row4)
#print align_cap1
#print align_cap2
#print align_cap3
reads_bwa = align_cap1[len(align_cap1)-1]
aligned_bwa = align_cap2[len(align_cap2)-1]
perc_align=str(float(align_cap3[len(align_cap3)-1])*100) 

for line in f3:
    if line.startswith('LIBRARY'):
        line = f3.next()
        line = line.split("\t");
        
        lib = line[0]
    
        PERCENT_DUPLICATION=line[7]
        ESTIMATED_LIBRARY_SIZE=line[8].split('\n')
        ESTIMATED_LIBRARY_SIZE=ESTIMATED_LIBRARY_SIZE[0]
#print ESTIMATED_LIBRARY_SIZE

raw_var_count = 0
for line in f4:
    if line.startswith('chr'):
        raw_var_count +=1 
#print raw_var_count

raw_snp_count = 0
for line in f5:
    if line.startswith('chr'):
        raw_snp_count +=1
#print raw_snp_count

raw_indel_count = 0
for line in f6:
    if line.startswith('chr'):
        raw_indel_count +=1
#print raw_indel_count

fileter_snp_count = 0
for line in f7:
    if line.startswith('chr'):
        line = line.split('\t')
        #print line
        #print line[6]
        if line[6] == 'PASS':
            fileter_snp_count += 1
#   fileter_snp_count += 1
#print fileter_snp_count

fileter_indel_count = 0
for line in f8:
    if line.startswith('chr'):
        line = line.split('\t')
        #print line                                                                                                                                                                  
        #print line[6]                                                                                                                                                               
        if line[6] == 'PASS':
            fileter_indel_count += 1
#   fileter_snp_count += 1                                                                                                                                                           
#print fileter_indel_count



print >> ofh ,'Pre_Trim',"\t",'Post_Trim',"\t",'Total_BWA_Reads',"\t",'Aligned_BWA_Reads',"\t",'%Aligned',"\t",'%Duplication',"\t",'LIBRARY_SIZE',"\t",'No_raw_var',"\t",'No_snps',"\t",'no_indel',"\t",'no_filtd_snps',"\t",'no_filtd_indels'
print >>ofh , total_reads+"\t\t"+trimmed_reads+"\t\t"+reads_bwa+"\t\t\t"+aligned_bwa+"\t\t\t"+perc_align+"\t\t"+PERCENT_DUPLICATION+"\t"+ESTIMATED_LIBRARY_SIZE+"\t"+str(raw_var_count)+"\t"+str(raw_snp_count)+"\t"+str(raw_indel_count)+"\t"+str(fileter_snp_count)+"\t"+str(fileter_indel_count)+"\t"+lib
