#!/bin/bash

#Abort on any error,
set -e

echo Running on $HOSTNAME
echo Current PATH is $PATH
#source ~/.bash_profile
#export PATH=$HOME/bin:$PATH
#echo Revised PATH is $PATH

#(and thus quit the script right away)
#which blastx

cd /p/path_to/Mcerasi_wild_samples


#x = """Mc_PR_Cherry_1_Read1.fastq
#Mc_PR_Cherry_1_Read2.fastq
#Mc_PR_Cherry_2_Read1.fastq
#Mc_PR_Cherry_2_Read1.fastq
#Mc_PR_Cherry_2_Read2.fastq
#Mc_PR_Cherry_3_Read1.fastq
#Mc_PR_Cherry_3_Read2.fastq
#Mc_PR_cress_1_Read1.fastq
#Mc_PR_cress_1_Read2.fastq
#Mc_PR_cress_2_Read1.fastq
#Mc_PR_cress_2_Read2.fastq
#Mc_PR_cress_3_Read1.fastq
#Mc_PR_cress_3_Read2.fastq
#Mc_PR_Galium_1_Read1.fastq
#Mc_PR_Galium_1_Read2.fastq
#Mc_PR_Galium_2_Read1.fastq
#Mc_PR_Galium_2_Read2.fastq
#Mc_PR_Galium_3_Read1.fastq
#Mc_PR_Galium_3_Read2.fastq
#""".split()


#for i in x:
#	print """\njava -jar /p/path_to/Downloads/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 8 -phred33 %s_Read1.fastq %s_Read2.fastq %s_R1_paired.fq %s_R1_unpaired.fq %s_R2_paired.fq %s_R2_unpaired.fq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 HEADCROP:9 TRAILING:3 SLIDINGWINDOW:4:22 MINLEN:36 \n\n#fastqc %s_R1_paired.fq %s_R2_paired.fq """ %(i[:-12],i[:-12],i[:-12],i[:-12],i[:-12],i[:-12],i[:-12],i[:-12])


java -jar /p/path_to/Downloads/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 8 -phred33 Mc_PR_Galium_2_R1_paired.fq.gz Mc_PR_Galium_2_R2_paired.fq.gz Mc_PR_Galium_2_R1_Q22paired.fq Mc_PR_Galium_2_R1_unpaired.fq Mc_PR_Galium_2_R2_Q22paired.fq Mc_PR_Galium_2_R2_unpaired.fq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 HEADCROP:9 TRAILING:3 SLIDINGWINDOW:4:22 MINLEN:36 

#fastqc Mc_PR_Galium_2_R1_paired.fq Mc_PR_Galium_2_R2_paired.fq 

java -jar /p/path_to/Downloads/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 8 -phred33 Mc_PR_cress_3_R1_paired.fq.gz Mc_PR_cress_3_R2_paired.fq.gz Mc_PR_cress_3_R1_Q22paired.fq Mc_PR_cress_3_R1_unpaired.fq Mc_PR_cress_3_R2_Q22paired.fq Mc_PR_cress_3_R2_unpaired.fq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 HEADCROP:9 TRAILING:3 SLIDINGWINDOW:4:22 MINLEN:36 

#fastqc Mc_PR_cress_3_R1_paired.fq Mc_PR_cress_3_R2_paired.fq 

java -jar /p/path_to/Downloads/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 8 -phred33 Mc_PR_Cherry_3_R1_paired.fq.gz Mc_PR_Cherry_3_R2_paired.fq.gz Mc_PR_Cherry_3_R1_Q22paired.fq Mc_PR_Cherry_3_R1_unpaired.fq Mc_PR_Cherry_3_R2_Q22paired.fq Mc_PR_Cherry_3_R2_unpaired.fq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 HEADCROP:9 TRAILING:3 SLIDINGWINDOW:4:22 MINLEN:36 

#fastqc Mc_PR_Cherry_3_R1_paired.fq Mc_PR_Cherry_3_R2_paired.fq 

java -jar /p/path_to/Downloads/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 8 -phred33 Mc_PR_Galium_3_R1_paired.fq.gz Mc_PR_Galium_3_R2_paired.fq.gz Mc_PR_Galium_3_R1_Q22paired.fq Mc_PR_Galium_3_R1_unpaired.fq Mc_PR_Galium_3_R2_Q22paired.fq Mc_PR_Galium_3_R2_unpaired.fq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 HEADCROP:9 TRAILING:3 SLIDINGWINDOW:4:22 MINLEN:36 

#fastqc Mc_PR_Galium_3_R1_paired.fq Mc_PR_Galium_3_R2_paired.fq 

java -jar /p/path_to/Downloads/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 8 -phred33 Mc_PR_cress_2_R1_paired.fq.gz Mc_PR_cress_2_R2_paired.fq.gz Mc_PR_cress_2_R1_Q22paired.fq Mc_PR_cress_2_R1_unpaired.fq Mc_PR_cress_2_R2_Q22paired.fq Mc_PR_cress_2_R2_unpaired.fq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 HEADCROP:9 TRAILING:3 SLIDINGWINDOW:4:22 MINLEN:36 

#fastqc Mc_PR_cress_2_R1_paired.fq Mc_PR_cress_2_R2_paired.fq 

java -jar /p/path_to/Downloads/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 8 -phred33 Mc_PR_Cherry_1_R1_paired.fq.gz Mc_PR_Cherry_1_R2_paired.fq.gz Mc_PR_Cherry_1_R1_Q22paired.fq Mc_PR_Cherry_1_R1_unpaired.fq Mc_PR_Cherry_1_R2_Q22paired.fq Mc_PR_Cherry_1_R2_unpaired.fq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 HEADCROP:9 TRAILING:3 SLIDINGWINDOW:4:22 MINLEN:36 

#fastqc Mc_PR_Cherry_1_R1_paired.fq Mc_PR_Cherry_1_R2_paired.fq 

java -jar /p/path_to/Downloads/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 8 -phred33 Mc_PR_Cherry_2_R1_paired.fq.gz Mc_PR_Cherry_2_R2_paired.fq.gz Mc_PR_Cherry_2_R1_Q22paired.fq Mc_PR_Cherry_2_R1_unpaired.fq Mc_PR_Cherry_2_R2_Q22paired.fq Mc_PR_Cherry_2_R2_unpaired.fq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 HEADCROP:9 TRAILING:3 SLIDINGWINDOW:4:22 MINLEN:36 

#fastqc Mc_PR_Cherry_2_R1_paired.fq Mc_PR_Cherry_2_R2_paired.fq 

java -jar /p/path_to/Downloads/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 8 -phred33 Mc_PR_cress_1_R1_paired.fq.gz Mc_PR_cress_1_R2_paired.fq.gz Mc_PR_cress_1_R1_Q22paired.fq Mc_PR_cress_1_R1_unpaired.fq Mc_PR_cress_1_R2_Q22paired.fq Mc_PR_cress_1_R2_unpaired.fq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 HEADCROP:9 TRAILING:3 SLIDINGWINDOW:4:22 MINLEN:36 

#fastqc Mc_PR_cress_1_R1_paired.fq Mc_PR_cress_1_R2_paired.fq 

java -jar /p/path_to/Downloads/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 8 -phred33 Mc_PR_Galium_1_R1_paired.fq.gz Mc_PR_Galium_1_R2_paired.fq.gz Mc_PR_Galium_1_R1_Q22paired.fq Mc_PR_Galium_1_R1_unpaired.fq Mc_PR_Galium_1_R2_Q22paired.fq Mc_PR_Galium_1_R2_unpaired.fq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 HEADCROP:9 TRAILING:3 SLIDINGWINDOW:4:22 MINLEN:36 

#fastqc Mc_PR_Galium_1_R1_paired.fq Mc_PR_Galium_1_R2_paired.fq 