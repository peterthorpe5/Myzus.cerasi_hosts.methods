#!/bin/bash
#$ -cwd

cd /home/PATH/M.cerasi


#!/bin/bash
#$ -cwd

cd /home/PATH/M.cerasi


mkdir star_indicies
STAR --runMode genomeGenerate --runThreadN 6 --limitGenomeGenerateRAM 74554136874 --genomeDir /home/PATH/M.cerasi/star_indicies --genomeFastaFiles /home/PATH/M.cerasi/Mc_alt.fasta

java -jar /home/PATH/Downloads/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 6 -phred33 Mc_PR_Galium_2_Read1.fastq.gz Mc_PR_Galium_2_Read2.fastq.gz Mc_PR_Galium_2_R1_Q22paired.fq.gz Mc_PR_Galium_2_R1_unpaired.fq.gz Mc_PR_Galium_2_R2_Q22paired.fq.gz Mc_PR_Galium_2_R2_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 HEADCROP:9 TRAILING:3 SLIDINGWINDOW:4:22 MINLEN:36 
cd /home/PATH/M.cerasi


STAR --genomeDir star_indicies/  --runThreadN 6 --sjdbGTFfile GTOOLS_exon.gtf --quantMode TranscriptomeSAM --outSAMmapqUnique 255 --outFilterMismatchNmax 10 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --readFilesIn /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_Galium_2_R1_Q22paired.fq.gz /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_Galium_2_R2_Q22paired.fq.gz --outFileNamePrefix Mc_PR_Galium_2

STAR --genomeDir star_indicies/  --runThreadN 6 --sjdbGTFfile GTOOLS_exon.gtf --quantMode TranscriptomeSAM --outSAMmapqUnique 255 --outFilterMismatchNmax 10 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --readFilesIn /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_Cherry_1_R1_Q22paired.fq.gz /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_Cherry_1_R2_Q22paired.fq.gz --outFileNamePrefix Mc_PR_Cherry_1

STAR --genomeDir star_indicies/  --runThreadN 6 --sjdbGTFfile GTOOLS_exon.gtf --quantMode TranscriptomeSAM --outSAMmapqUnique 255 --outFilterMismatchNmax 10 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --readFilesIn /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_Cherry_3_R1_Q22paired.fq.gz /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_Cherry_3_R2_Q22paired.fq.gz --outFileNamePrefix Mc_PR_Cherry_3

STAR --genomeDir star_indicies/  --runThreadN 6 --sjdbGTFfile GTOOLS_exon.gtf --quantMode TranscriptomeSAM --outSAMmapqUnique 255 --outFilterMismatchNmax 10 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --readFilesIn /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_cress_2_R1_Q22paired.fq.gz /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_cress_2_R2_Q22paired.fq.gz --outFileNamePrefix Mc_PR_cress_2

STAR --genomeDir star_indicies/  --runThreadN 6 --sjdbGTFfile GTOOLS_exon.gtf --quantMode TranscriptomeSAM --outSAMmapqUnique 255 --outFilterMismatchNmax 10 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --readFilesIn /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_Galium_1_R1_Q22paired.fq.gz /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_Galium_1_R2_Q22paired.fq.gz --outFileNamePrefix Mc_PR_Galium_1

STAR --genomeDir star_indicies/  --runThreadN 6 --sjdbGTFfile GTOOLS_exon.gtf --quantMode TranscriptomeSAM --outSAMmapqUnique 255 --outFilterMismatchNmax 10 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --readFilesIn /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_Galium_3_R1_Q22paired.fq.gz /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_Galium_3_R2_Q22paired.fq.gz --outFileNamePrefix Mc_PR_Galium_3

STAR --genomeDir star_indicies/  --runThreadN 6 --sjdbGTFfile GTOOLS_exon.gtf --quantMode TranscriptomeSAM --outSAMmapqUnique 255 --outFilterMismatchNmax 10 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --readFilesIn /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_Cherry_1_R1_Q22paired.fq.gz /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_Cherry_1_R2_Q22paired.fq.gz --outFileNamePrefix Mc_PR_Cherry_1

STAR --genomeDir star_indicies/  --runThreadN 6 --sjdbGTFfile GTOOLS_exon.gtf --quantMode TranscriptomeSAM --outSAMmapqUnique 255 --outFilterMismatchNmax 10 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --readFilesIn /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_Cherry_3_R1_Q22paired.fq.gz /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_Cherry_3_R2_Q22paired.fq.gz --outFileNamePrefix Mc_PR_Cherry_3

STAR --genomeDir star_indicies/  --runThreadN 6 --sjdbGTFfile GTOOLS_exon.gtf --quantMode TranscriptomeSAM --outSAMmapqUnique 255 --outFilterMismatchNmax 10 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --readFilesIn /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_cress_2_R1_Q22paired.fq.gz /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_cress_2_R2_Q22paired.fq.gz --outFileNamePrefix Mc_PR_cress_2

STAR --genomeDir star_indicies/  --runThreadN 6 --sjdbGTFfile GTOOLS_exon.gtf --quantMode TranscriptomeSAM --outSAMmapqUnique 255 --outFilterMismatchNmax 10 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --readFilesIn /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_Galium_1_R1_Q22paired.fq.gz /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_Galium_1_R2_Q22paired.fq.gz --outFileNamePrefix Mc_PR_Galium_1

STAR --genomeDir star_indicies/  --runThreadN 6 --sjdbGTFfile GTOOLS_exon.gtf --quantMode TranscriptomeSAM --outSAMmapqUnique 255 --outFilterMismatchNmax 10 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --readFilesIn /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_Cherry_2_R1_Q22paired.fq.gz /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_Cherry_2_R2_Q22paired.fq.gz --outFileNamePrefix Mc_PR_Cherry_2

STAR --genomeDir star_indicies/  --runThreadN 6 --sjdbGTFfile GTOOLS_exon.gtf --quantMode TranscriptomeSAM --outSAMmapqUnique 255 --outFilterMismatchNmax 10 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --readFilesIn /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_cress_1_R1_Q22paired.fq.gz /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_cress_1_R2_Q22paired.fq.gz --outFileNamePrefix Mc_PR_cress_1

STAR --genomeDir star_indicies/  --runThreadN 6 --sjdbGTFfile GTOOLS_exon.gtf --quantMode TranscriptomeSAM --outSAMmapqUnique 255 --outFilterMismatchNmax 10 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --readFilesIn /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_cress_3_R1_Q22paired.fq.gz /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_cress_3_R2_Q22paired.fq.gz --outFileNamePrefix Mc_PR_cress_3

STAR --genomeDir star_indicies/  --runThreadN 6 --sjdbGTFfile GTOOLS_exon.gtf --quantMode TranscriptomeSAM --outSAMmapqUnique 255 --outFilterMismatchNmax 10 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --readFilesIn /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_Galium_2_R1_Q22paired.fq.gz /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_Galium_2_R2_Q22paired.fq.gz --outFileNamePrefix Mc_PR_Galium_2

STAR --genomeDir star_indicies/  --runThreadN 6 --sjdbGTFfile GTOOLS_exon.gtf --quantMode TranscriptomeSAM --outSAMmapqUnique 255 --outFilterMismatchNmax 10 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --readFilesIn /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_Cherry_2_R1_Q22paired.fq.gz /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_Cherry_2_R2_Q22paired.fq.gz --outFileNamePrefix Mc_PR_Cherry_2

STAR --genomeDir star_indicies/  --runThreadN 6 --sjdbGTFfile GTOOLS_exon.gtf --quantMode TranscriptomeSAM --outSAMmapqUnique 255 --outFilterMismatchNmax 10 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --readFilesIn /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_cress_1_R1_Q22paired.fq.gz /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_cress_1_R2_Q22paired.fq.gz --outFileNamePrefix Mc_PR_cress_1

STAR --genomeDir star_indicies/  --runThreadN 6 --sjdbGTFfile GTOOLS_exon.gtf --quantMode TranscriptomeSAM --outSAMmapqUnique 255 --outFilterMismatchNmax 10 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --readFilesIn /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_cress_3_R1_Q22paired.fq.gz /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_cress_3_R2_Q22paired.fq.gz --outFileNamePrefix Mc_PR_cress_3

STAR --genomeDir star_indicies/  --runThreadN 6 --sjdbGTFfile GTOOLS_exon.gtf --quantMode TranscriptomeSAM --outSAMmapqUnique 255 --outFilterMismatchNmax 10 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --readFilesIn /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_Galium_3_R1_Q22paired.fq.gz /mnt/shared/projects/aphids/201406_Mcerasi_RNAseq_cherry_cress_gallium/Mc_PR_Galium_3_R2_Q22paired.fq.gz --outFileNamePrefix Mc_PR_Galium_3

STAR --genomeDir star_indicies/  --runThreadN 6 --sjdbGTFfile GTOOLS_exon.gtf --quantMode TranscriptomeSAM --outSAMmapqUnique 255 --outFilterMismatchNmax 10 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --readFilesIn Mc_all_R1.fq.gz Mc_all_R2.fq.gz --outFileNamePrefix Mc_all

#
#

samtools sort -@ 6 Mc_PR_Galium_2Aligned.sortedByCoord.out.bam Mc_PR_Galium_2

samtools sort -@ 6 Mc_PR_Cherry_1Aligned.sortedByCoord.out.bam Mc_PR_Cherry_1
#
samtools sort -@ 6 Mc_PR_Cherry_2Aligned.sortedByCoord.out.bam Mc_PR_Cherry_2
#
samtools sort -@ 6 Mc_PR_Cherry_3Aligned.sortedByCoord.out.bam Mc_PR_Cherry_3
#
samtools sort -@ 6 Mc_PR_cress_1Aligned.sortedByCoord.out.bam Mc_PR_cress_1
#
samtools sort -@ 6 Mc_PR_cress_2Aligned.sortedByCoord.out.bam Mc_PR_cress_2
#
samtools sort -@ 6 Mc_PR_cress_3Aligned.sortedByCoord.out.bam Mc_PR_cress_3
#
samtools sort -@ 6 Mc_PR_Galium_1Aligned.sortedByCoord.out.bam Mc_PR_Galium_1
#
samtools sort -@ 6 Mc_PR_Galium_2Aligned.sortedByCoord.out.bam Mc_PR_Galium_2
#
samtools sort -@ 6 Mc_PR_Galium_3Aligned.sortedByCoord.out.bam Mc_PR_Galium_3

#
#
samtools index Mc_PR_Cherry_1.bam

samtools index Mc_PR_cress_1.bam

samtools index Mc_PR_cress_3.bam

samtools index Mc_PR_Galium_2.bam

samtools index Mc_PR_Cherry_2.bam

samtools index Mc_PR_cress_2.bam

samtools index Mc_PR_Galium_1.bam

samtools index Mc_PR_Galium_3.bam


#Mc_PR_Cherry_1.bam Mc_PR_Cherry_2.bam Mc_PR_Cherry_3.bam Mc_PR_Galium_1.bam Mc_PR_Galium_2.bam Mc_PR_Galium_3.bam Mc_PR_cress_1.bam Mc_PR_cress_2.bam Mc_PR_cress_3.bam

samtools merge -@ 6 -f -r -h header.txt merged.bam Mc_PR_Cherry_1.bam Mc_PR_Cherry_2.bam Mc_PR_Cherry_3.bam Mc_PR_Galium_1.bam Mc_PR_Galium_2.bam Mc_PR_Galium_3.bam Mc_PR_cress_1.bam Mc_PR_cress_2.bam Mc_PR_cress_3.bam

samtools index merged.bam


# call the SNP/ varients then calculate the density
/mnt/shared/synology/PATH/apps/conda/bin/freebayes -f Mc_alt.fasta --no-population-priors --min-alternate-count 5 --min-alternate-fraction 0.9 -v M.cerasi_Mc_PR_Cherry_1.vcf --ploidy 2 Mc_PR_Cherry_1.bam
vcftools --vcf M.cerasi_Mc_PR_Cherry_1.vcf --SNPdensity 10000 --out M.cerasi_Mc_PR_Cherry_1.SNP_snpdensity

/mnt/shared/synology/PATH/apps/conda/bin/freebayes -f Mc_alt.fasta --no-population-priors --min-alternate-count 5 --min-alternate-fraction 0.9 -v M.cerasi_Mc_PR_Cherry_3.vcf --ploidy 2 Mc_PR_Cherry_3.bam
vcftools --vcf M.cerasi_Mc_PR_Cherry_3.vcf --SNPdensity 10000 --out M.cerasi_Mc_PR_Cherry_3.SNP_snpdensity

/mnt/shared/synology/PATH/apps/conda/bin/freebayes -f Mc_alt.fasta --no-population-priors --min-alternate-count 5 --min-alternate-fraction 0.9 -v M.cerasi_Mc_PR_cress_2.vcf --ploidy 2 Mc_PR_cress_2.bam
vcftools --vcf M.cerasi_Mc_PR_cress_2.vcf --SNPdensity 10000 --out M.cerasi_Mc_PR_cress_2.SNP_snpdensity

/mnt/shared/synology/PATH/apps/conda/bin/freebayes -f Mc_alt.fasta --no-population-priors --min-alternate-count 5 --min-alternate-fraction 0.9 -v M.cerasi_Mc_PR_Galium_1.vcf --ploidy 2 Mc_PR_Galium_1.bam
vcftools --vcf M.cerasi_Mc_PR_Galium_1.vcf --SNPdensity 10000 --out M.cerasi_Mc_PR_Galium_1.SNP_snpdensity

/mnt/shared/synology/PATH/apps/conda/bin/freebayes -f Mc_alt.fasta --no-population-priors --min-alternate-count 5 --min-alternate-fraction 0.9 -v M.cerasi_Mc_PR_Galium_3.vcf --ploidy 2 Mc_PR_Galium_3.bam
vcftools --vcf M.cerasi_Mc_PR_Galium_3.vcf --SNPdensity 10000 --out M.cerasi_Mc_PR_Galium_3.SNP_snpdensity

/mnt/shared/synology/PATH/apps/conda/bin/freebayes -f Mc_alt.fasta --no-population-priors --min-alternate-count 5 --min-alternate-fraction 0.9 -v M.cerasi_Mc_PR_Cherry_2.vcf --ploidy 2 Mc_PR_Cherry_2.bam
vcftools --vcf M.cerasi_Mc_PR_Cherry_2.vcf --SNPdensity 10000 --out M.cerasi_Mc_PR_Cherry_2.SNP_snpdensity

/mnt/shared/synology/PATH/apps/conda/bin/freebayes -f Mc_alt.fasta --no-population-priors --min-alternate-count 5 --min-alternate-fraction 0.9 -v M.cerasi_Mc_PR_cress_1.vcf --ploidy 2 Mc_PR_cress_1.bam
vcftools --vcf M.cerasi_Mc_PR_cress_1.vcf --SNPdensity 10000 --out M.cerasi_Mc_PR_cress_1.SNP_snpdensity

/mnt/shared/synology/PATH/apps/conda/bin/freebayes -f Mc_alt.fasta --no-population-priors --min-alternate-count 5 --min-alternate-fraction 0.9 -v M.cerasi_Mc_PR_cress_3.vcf --ploidy 2 Mc_PR_cress_3.bam
vcftools --vcf M.cerasi_Mc_PR_cress_3.vcf --SNPdensity 10000 --out M.cerasi_Mc_PR_cress_3.SNP_snpdensity

/mnt/shared/synology/PATH/apps/conda/bin/freebayes -f Mc_alt.fasta --no-population-priors --min-alternate-count 5 --min-alternate-fraction 0.9 -v M.cerasi_Mc_PR_Galium_2.vcf --ploidy 2 Mc_PR_Galium_2.bam
vcftools --vcf M.cerasi_Mc_PR_Galium_2.vcf --SNPdensity 10000 --out M.cerasi_Mc_PR_Galium_2.SNP_snpdensity


/mnt/shared/synology/PATH/apps/conda/bin/freebayes -f Mc_alt.fasta --no-population-priors --min-alternate-count 5 --min-alternate-fraction 0.9 -v M.cerasi_cherry_galium_cress.vcf --ploidy 2 merged.bam

# caluclate the stats on the vcf file

/storage/home/users/pjt6/rtg-tools-3.10.1/rtg vcfstats --allele-lengths M.cerasi_cherry_galium_cress.vcf > Aphid_vcf_SNP_stats.txt

# claulate the PI per aphid/host population-priors
# the three rep on cherry
vcftools --vcf M.cerasi_cherry_galium_cress.vcf --window-pi 10000 --out M.cerasi_cherry.10000.sitepi --remove-indv  Mc_PR_cress_1  --remove-indv  Mc_PR_cress_2   --remove-indv  Mc_PR_cress_3  --remove-indv  Mc_PR_Galium_1  --remove-indv  Mc_PR_Galium_2 --remove-indv  Mc_PR_Galium_3

# the three reps on cress
vcftools --vcf M.cerasi_cherry_galium_cress.vcf --window-pi 10000 --out M.cerasi_cress.10000.sitepi --remove-indv  Mc_PR_Cherry_1  --remove-indv  Mc_PR_Cherry_2  --remove-indv  Mc_PR_Cherry_3  --remove-indv  Mc_PR_Galium_1  --remove-indv  Mc_PR_Galium_2 --remove-indv  Mc_PR_Galium_3

# the three reps on cleavers/ galium (sticky weed - what ever it is called)
vcftools --vcf M.cerasi_cherry_galium_cress.vcf --window-pi 10000 --out M.cerasi_cleavers.10000.sitepi  --remove-indv  Mc_PR_Cherry_1  --remove-indv  Mc_PR_Cherry_2  --remove-indv  Mc_PR_Cherry_3  --remove-indv  Mc_PR_cress_1  --remove-indv  Mc_PR_cress_2   --remove-indv  Mc_PR_cress_3  

# The resulting files were then parsed and analysed with 
python compare_genetic_diversity.py

# to investigate expected versus observer heterozygosity

vcftools --vcf M.cerasi_cherry_galium_cress.vcf --het --out M.cerasi_cherry_galium_cress.het

