#!/bin/bash
#$ -cwd

cd /home/pt40963/scratch/M.cerasi

#rnaseq quast

python ~/Downloads/rnaQUAST-1.2.0/rnaQUAST.py -t 16 --blat -c Trinity-GG002_trans_v2.fasta -r Mc_alt.fasta -gtf GTOOLS.gtf -1 Mc_all_R1.fq.gz -2 Mc_all_R2.fq.gz -o Mc_wild_rnaQUAST

