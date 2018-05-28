#$ -l hostname="n13*"
#$ -cwd
cd /p/path_to/M.cerasi_genome_data/final_assembly_v1

mkdir star_indicies

# index the genome

STAR --genomeDir star_indicies/  --runThreadN 4 --outFilterMultimapNmax 3 --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 2 --outFilterMismatchNoverLmax 0.1 --readFilesCommand zcat --readFilesIn /p/path_to/M.cerasi_genome_data/RNAseq_data/Mc_R1.fq.gz /p/path_to/M.cerasi_genome_data/RNAseq_data/Mc_R2.fq.gz --outFileNamePrefix Mc_strict_RNAseq_mapped
 

# host samples


STAR --genomeDir star_indicies/  --runThreadN 4 --outFilterMultimapNmax 5 --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7 --readFilesCommand zcat --outFileNamePrefix /p/path_to/M.cerasi_genome_data/final_assembly_v1/Mc_Gallium1 --readFilesIn /p/path_to/Mcerasi_wild_samples/Mc_PR_Galium_1_R1_Q22paired.fq.gz /p/path_to/Mcerasi_wild_samples/Mc_PR_Galium_1_R2_Q22paired.fq.gz

STAR --genomeDir star_indicies/  --runThreadN 4 --outFilterMultimapNmax 5 --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7 --readFilesCommand zcat --outFileNamePrefix /p/path_to/M.cerasi_genome_data/final_assembly_v1/Mc_Gallium2 --readFilesIn /p/path_to/Mcerasi_wild_samples/Mc_PR_Galium_2_R1_Q22paired.fq.gz /p/path_to/Mcerasi_wild_samples/Mc_PR_Galium_2_R2_Q22paired.fq.gz

STAR --genomeDir star_indicies/  --runThreadN 4 --outFilterMultimapNmax 5 --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7 --readFilesCommand zcat --outFileNamePrefix /p/path_to/M.cerasi_genome_data/final_assembly_v1/Mc_Gallium3 --readFilesIn /p/path_to/Mcerasi_wild_samples/Mc_PR_Galium_3_R1_Q22paired.fq.gz /p/path_to/Mcerasi_wild_samples/Mc_PR_Galium_3_R2_Q22paired.fq.gz

STAR --genomeDir star_indicies/  --runThreadN 4 --outFilterMultimapNmax 5 --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7 --readFilesCommand zcat --outFileNamePrefix /p/path_to/M.cerasi_genome_data/final_assembly_v1/Mc_cress1 --readFilesIn /p/path_to/Mcerasi_wild_samples/Mc_PR_cress_1_R1_Q22paired.fq.gz /p/path_to/Mcerasi_wild_samples/Mc_PR_cress_1_R2_Q22paired.fq.gz

STAR --genomeDir star_indicies/  --runThreadN 6 --outFilterMultimapNmax 5 --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7 --readFilesCommand zcat --outFileNamePrefix /p/path_to/M.cerasi_genome_data/final_assembly_v1/Mc_cress2 --readFilesIn /p/path_to/Mcerasi_wild_samples/Mc_PR_cress_2_R1_Q22paired.fq.gz /p/path_to/Mcerasi_wild_samples/Mc_PR_cress_2_R2_Q22paired.fq.gz

STAR --genomeDir star_indicies/  --runThreadN 6 --outFilterMultimapNmax 5 --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7 --readFilesCommand zcat --outFileNamePrefix /p/path_to/M.cerasi_genome_data/final_assembly_v1/Mc_cress3 --readFilesIn /p/path_to/Mcerasi_wild_samples/Mc_PR_cress_3_R1_Q22paired.fq.gz /p/path_to/Mcerasi_wild_samples/Mc_PR_cress_3_R2_Q22paired.fq.gz

STAR --genomeDir star_indicies/  --runThreadN 6 --outFilterMultimapNmax 5 --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7 --readFilesCommand zcat --outFileNamePrefix /p/path_to/M.cerasi_genome_data/final_assembly_v1/Mc_cherry2 --readFilesIn /p/path_to/Mcerasi_wild_samples/Mc_PR_Cherry_2_R1_Q22paired.fq.gz /p/path_to/Mcerasi_wild_samples/Mc_PR_Cherry_2_R2_Q22paired.fq.gz

STAR --genomeDir star_indicies/  --runThreadN 6 --outFilterMultimapNmax 5 --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7 --readFilesCommand zcat --outFileNamePrefix /p/path_to/M.cerasi_genome_data/final_assembly_v1/Mc_cherry3 --readFilesIn /p/path_to/Mcerasi_wild_samples/Mc_PR_Cherry_3_R1_Q22paired.fq.gz /p/path_to/Mcerasi_wild_samples/Mc_PR_Cherry_3_R2_Q22paired.fq.gz

STAR --genomeDir star_indicies/  --runThreadN 6 --outFilterMultimapNmax 5 --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7 --readFilesCommand zcat --outFileNamePrefix /p/path_to/M.cerasi_genome_data/final_assembly_v1/Mc_cherry1 --readFilesIn /p/path_to/Mcerasi_wild_samples/Mc_PR_Cherry_1_R1_Q22paired.fq.gz /p/path_to/Mcerasi_wild_samples/Mc_PR_Cherry_1_R2_Q22paired.fq.gz

