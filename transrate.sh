#!/bin/bash
#$ -cwd

cd /p/path_to/M.cerasi_genome_data/final_assembly_v1


# transrate command
transrate --assembly Myzus_cerasi_genome_assembly.v1.fasta --left /p/path_to/M.cerasi_genome_data/RNAseq_data/Mc_R1_paired.fq.gz --right /p/path_to/M.cerasi_genome_data/RNAseq_data/Mc_R2_paired.fq.gz --threads 4

#bam read to identify chimeras. The probabilty that it is one gene - not more than one fused. 
#~/transrate-tools/src/bam-read Rp_sorted.bam chimera_indetification_outfile_95.csv 0.95

#~/transrate-tools/src/bam-read Rp_sorted.bam chimera_indetification_outfile_99.csv 0.99