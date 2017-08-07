#!/bin/bash


####qsub -pe "*" 8 shell_script_trial.sh
## command to run the shell script...  the pe is parallel environment, the * signifies any environment.  
##  I was getting error "bowtie 2:  command not found because even though I had added the path for bowtie, the grid had no idea.  Sooooo, when I added
##  the full path below: voila! 

##  There are also different ways to submit a job.  Straightfoward in the .sh script like /usr/global/blp/bowtie2-2.1.0/bowtie2 -x sk36_index -U inputfile -S outputfile
##  Or if running multiple jobs with same script you can do it as below in the .sh script with $1 and $2 being accepted at commandline with qsub.  The 
##  file names must be in separate quotes like the next line.
## qsub -pe "*" 8 shell_script_trial.sh "PingPool_CAGATC_R1.fastq" "pur_959"



######################################

#$ -cwd
#$ -e ./
#$ -o ./
#$ -S /bin/bash
#$ -N Nihar_qsub_example
######################################
FILE1=$1
FILE2=$2 
SAMPLE=$3

/home/bradleysp/bin/fastqstats $FILE1 > ${FILE1}.stats.txt
/home/bradleysp/bin/fastqstats $FILE2 > ${FILE2}.stats.txt
 
 
