#!/bin/bash

#Abort on any error,
set -e

echo Running on $HOSTNAME
echo Current PATH is $PATH
#source ~/.bash_profile

#export PATH=$HOME/bin:$PATH
#echo Revised PATH is $PATH

#This will give an error if bowtie is not on the $PATH
#(and thus quit the script right away)
#which bowtie

#cd /p/path_to/Mcerasi_wild_samples/
#gunzip *.fq.gz
cd /p/path_to/Mcerasi_wild_samples/RNAseq_final_GG_v1
echo About to run Trinity!
ulimit -s unlimited
wait
# normalise by kmer coverage
#/p/path_to/Downloads/trinityrnaseq_r20131110/util/normalize_by_kmer_coverage.pl --seqType fq --JM 230G --KMER_SIZE 31 --max_cov 50 --left /p/path_to/Mcerasi_wild_samples/Mc_R1.fq --right /p/path_to/Mcerasi_wild_samples/Mc_R2.fq --output Mc_normalised --pairs_together --PARALLEL_STATS
echo Trinity-norm is done

# trinity assembly
#/p/path_to/Downloads/trinityrnaseq_r20140717/Trinity --seqType fq --JM 6G --left /p/path_to/Mcerasi_wild_samples/Mc_R1.fq.normalized_K31_C50_pctSD200.fq --right /p/path_to/Mcerasi_wild_samples/Mc_R2.fq.normalized_K31_C50_pctSD200.fq --min_contig_length 175 --genome /p/path_to/M.cerasi_genome_data/final_assembly_v1/Myzus_cerasi_genome_assembly.v1.fasta --genome_guided_max_intron 12000 --genome_guided_sort_buffer 5G --full_cleanup --genome_guided_CPU 2 --GMAP_CPU 14 --CPU 7 --min_kmer_cov 2 --SS_lib_type FR --output Mcerasi_GG_assembly


echo trinity is finished...!!!!!!!!!!!!!!!!!!
cd /p/path_to/Mcerasi_wild_samples/RNAseq_final_GG_v1/
#stats
#/p/path_to/Downloads/trinityrnaseq_r20140717/util/TrinityStats.pl Trinity-GG.fasta >trinity_assembly_stats.txt

#reformat fasta names for transrate
#python ~/misc_python/NGS/reformat_fasta_gene_names_for_transrate.py Trinity-GG.fasta Trinity-GG002.fasta

# transrate
#transrate --assembly Trinity-GG002.fasta --left /p/path_to/Mcerasi_wild_samples/Mc_R1.fq --right /p/path_to/Mcerasi_wild_samples/Mc_R2.fq --threads 4

#transrate --assembly Trinity-GG002.fasta --left /p/path_to/Mcerasi_wild_samples/Mc_R1.fq --right /p/path_to/Mcerasi_wild_samples/Mc_R2.fq --threads 4

#transrate --assembly Trinity-GG002_trans_v1.fasta --left /p/path_to/Mcerasi_wild_samples/Mc_R1.fq --right /p/path_to/Mcerasi_wild_samples/Mc_R2.fq --threads 8

#/p/path_to/Downloads/trinityrnaseq_r20140717/util/TrinityStats.pl Trinity-GG002_trans_v1.fasta > Trinity-GG002_trans_v1.fasta.stats.txt


# reformat names for trinity mapping:

#python ~/misc_python/NGS/reformat_fasta_gene_names_for_RSEM_after_transrate.py Trinity-GG002_trans_v1.fasta test.fasta

#/p/path_to/Downloads/trinityrnaseq_r20140717/util/TrinityStats.pl test.fasta > test.fasta.stats.fasta

#/p/path_to/Downloads/trinityrnaseq_r20140717/util/TrinityStats.pl Trinity-GG002.fasta >Trinity-GG002.fasta.stats.txt

#/p/path_to/Downloads/trinityrnaseq_r20140717/util/TrinityStats.pl Trinity-GG002_trans_v1.fasta > Trinity-GG002_trans_v1.fasta.stats.txt

#/p/path_to/Downloads/trinityrnaseq_r20140717/util/TrinityStats.pl Trinity-GG002_trans_v2.fasta > Trinity-GG002_trans_v2.fasta.stats.txt


#python ~/misc_python/NGS/reformat_fasta_gene_names_for_RSEM_after_transrate.py Trinity-GG002_trans_v2.fasta final_Mc_wild_GG.fasta

#/p/path_to/Downloads/trinityrnaseq_r20140717/util/TrinityStats.pl final_Mc_wild_GG.fasta > final_Mc_wild_GG.fasta.stats.txt

cd /p/path_to/Mcerasi_wild_samples/RNAseq_final_GG_v1/Mcerasi_GG_assembly 


# RSEM - get expression
/p/path_to/Downloads/trinityrnaseq_r20140717/util/align_and_estimate_abundance.pl --prep_reference --thread_count 8  --transcripts final_Mc_wild_GG.fasta --trinity_mode --est_method RSEM --aln_method bowtie --seqType fq --left /p/path_to/Mcerasi_wild_samples/Mc_PR_Cherry_1_R1_paired.fq --right /p/path_to/Mcerasi_wild_samples/Mc_PR_Cherry_1_R2_paired.fq --output_prefix Mc_cherry1

wait
/p/path_to/Downloads/trinityrnaseq_r20140717/util/align_and_estimate_abundance.pl --thread_count 6  --transcripts final_Mc_wild_GG.fasta --trinity_mode --est_method RSEM --aln_method bowtie --seqType fq --left /p/path_to/Mcerasi_wild_samples/Mc_PR_Cherry_2_R1_paired.fq --right /p/path_to/Mcerasi_wild_samples/Mc_PR_Cherry_2_R2_paired.fq --output_prefix Mc_cherry2

/p/path_to/Downloads/trinityrnaseq_r20140717/util/align_and_estimate_abundance.pl --thread_count 8  --transcripts final_Mc_wild_GG.fasta --trinity_mode --est_method RSEM --aln_method bowtie --seqType fq --left /p/path_to/Mcerasi_wild_samples/Mc_PR_Cherry_3_R1_paired.fq --right /p/path_to/Mcerasi_wild_samples/Mc_PR_Cherry_3_R2_paired.fq --output_prefix Mc_cherry3

/p/path_to/Downloads/trinityrnaseq_r20140717/util/align_and_estimate_abundance.pl --thread_count 8  --transcripts final_Mc_wild_GG.fasta --trinity_mode --est_method RSEM --aln_method bowtie --seqType fq --left /p/path_to/Mcerasi_wild_samples/Mc_PR_cress_1_R1_paired.fq --right /p/path_to/Mcerasi_wild_samples/Mc_PR_cress_1_R2_paired.fq --output_prefix Mc_cress1
                                                                                                                                                                                                                                                     
/p/path_to/Downloads/trinityrnaseq_r20140717/util/align_and_estimate_abundance.pl --thread_count 8  --transcripts final_Mc_wild_GG.fasta --trinity_mode --est_method RSEM --aln_method bowtie --seqType fq --left /p/path_to/Mcerasi_wild_samples/Mc_PR_cress_2_R1_paired.fq --right /p/path_to/Mcerasi_wild_samples/Mc_PR_cress_2_R2_paired.fq --output_prefix Mc_cress2
                                                                                                                                                                                                                       
#echo ...........................................................done3 in 2nd bacth
/p/path_to/Downloads/trinityrnaseq_r20140717/util/align_and_estimate_abundance.pl --thread_count 8  --transcripts final_Mc_wild_GG.fasta --trinity_mode --est_method RSEM --aln_method bowtie --seqType fq --left /p/path_to/Mcerasi_wild_samples/Mc_PR_cress_3_R1_paired.fq --right /p/path_to/Mcerasi_wild_samples/Mc_PR_cress_3_R2_paired.fq --output_prefix Mc_cress3

/p/path_to/Downloads/trinityrnaseq_r20140717/util/align_and_estimate_abundance.pl --thread_count 8  --transcripts final_Mc_wild_GG.fasta --trinity_mode --est_method RSEM --aln_method bowtie --seqType fq --left /p/path_to/Mcerasi_wild_samples/Mc_PR_Galium_1_R1_paired.fq --right /p/path_to/Mcerasi_wild_samples/Mc_PR_Galium_1_R2_paired.fq --output_prefix Mc_Gallium1
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
/p/path_to/Downloads/trinityrnaseq_r20140717/util/align_and_estimate_abundance.pl --thread_count 8  --transcripts final_Mc_wild_GG.fasta --trinity_mode --est_method RSEM --aln_method bowtie --seqType fq --left /p/path_to/Mcerasi_wild_samples/Mc_PR_Galium_2_R1_paired.fq --right /p/path_to/Mcerasi_wild_samples/Mc_PR_Galium_2_R2_paired.fq --output_prefix Mc_Gallium2

/p/path_to/Downloads/trinityrnaseq_r20140717/util/align_and_estimate_abundance.pl --thread_count 8  --transcripts final_Mc_wild_GG.fasta --trinity_mode --est_method RSEM --aln_method bowtie --seqType fq --left /p/path_to/Mcerasi_wild_samples/Mc_PR_Galium_3_R1_paired.fq --right /p/path_to/Mcerasi_wild_samples/Mc_PR_Galium_3_R2_paired.fq --output_prefix Mc_Gallium3

echo ................. all mapping is done .......................................
/p/path_to/Downloads/trinityrnaseq_r20131110/util/RSEM_util/merge_RSEM_frag_counts_single_table.pl Mc_cherry1.isoforms.results Mc_cherry2.isoforms.results Mc_cherry3.isoforms.results Mc_Gallium1.isoforms.results Mc_Gallium2.isoforms.results Mc_Gallium3.isoforms.results Mc_cress1.isoforms.results Mc_cress2.isoforms.results Mc_cress3.isoforms.results > Mc_wild.transcripts.counts.matrix
/p/path_to/Downloads/trinityrnaseq_r20131110/util/RSEM_util/merge_RSEM_frag_counts_single_table.pl Mc_cherry1.genes.results Mc_cherry2.genes.results Mc_cherry3.genes.results Mc_Gallium1.genes.results Mc_Gallium2.genes.results Mc_Gallium3.genes.results Mc_cress1.genes.results Mc_cress2.genes.results Mc_cress3.genes.results > Mc_wild.genes.counts.matrix
echo .............................................nearly everything done


#................. all mapping is done .......................................
/p/path_to/Downloads/trinityrnaseq_r20131110/util/RSEM_util/merge_RSEM_frag_counts_single_table.pl Mc_Gallium1.isoforms.results Mc_Gallium2.isoforms.results Mc_Gallium3.isoforms.results Mc_cress1.isoforms.results Mc_cress2.isoforms.results Mc_cress3.isoforms.results > Mc_wild_cress_gallium.transcripts.counts.matrix

/p/path_to/Downloads/trinityrnaseq_r20131110/util/RSEM_util/merge_RSEM_frag_counts_single_table.pl Mc_Gallium1.genes.results Mc_Gallium2.genes.results Mc_Gallium3.genes.results Mc_cress1.genes.results Mc_cress2.genes.results Mc_cress3.genes.results > Mc_wild_cress_galium.genes.counts.matrix


######################################################################################################################################################
# CDS prediction and annotation:

#!/bin/bash
#$ -l hostname="n13*"
#Abort on any error,
#set -e
echo Running on $HOSTNAME
echo Current PATH is $PATH
#########################################################################################################
##### shell script to annotate mpO using pfamAB domains	######################################
######################################################################################################

cd /p/path_to/Mcerasi_wild_samples/RNAseq_final_GG_v1/Mcerasi_GG_assembly/ 


#python ~/misc_python/NGS/reformat_fasta_gene_names_for_final_transcriptome.py final_Mc_wild_GG.fasta Mc_wild Mc_wild_RNAseq.fasta


# 				transdecoder picking the cds, guided by swiss prot and pfam domains

/p/path_to/Downloads/TransDecoder-2.0.1/TransDecoder.LongOrfs -t /p/path_to/Mcerasi_wild_samples/RNAseq_final_GG_v1/Mcerasi_GG_assembly/Mc_wild_RNAseq.fasta




#blastp -query /p/path_to/nematode/longi/Mc_wild_RNAseq_ass.fa.transdecoder_dir/longest_orfs.pep  -db /p/path_to/Downloads/TransDecoder-2.0.1/uniprot_sprot.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 8 > blastp.outfmt6
# use diamond instead - faster
diamond blastp -k 1 -p 8 --sensitive -v -q /p/path_to/Mcerasi_wild_samples/RNAseq_final_GG_v1/Mcerasi_GG_assembly/Mc_wild_RNAseq.fasta.transdecoder_dir/longest_orfs.pep -d /p/path_to/Downloads/blast_databases/uniref90.dmnd -a longestORF.blastp.uniref.DIAMOND.diamond
diamond view -a longestORF.blastp.uniref.DIAMOND.*.daa -f tab -o longestORF.blastp.uniref.DIAMOND.tab




hmmscan --cpu 8 --domtblout pfam.domtblout /p/path_to/Downloads/TransDecoder_r20131117/pfam/Pfam-AB.hmm.bin /p/path_to/Mcerasi_wild_samples/RNAseq_final_GG_v1/Mcerasi_GG_assembly/Mc_wild_RNAseq.fasta.transdecoder_dir/longest_orfs.pep

/p/path_to/Downloads/TransDecoder-2.0.1/TransDecoder.Predict -t /p/path_to/Mcerasi_wild_samples/RNAseq_final_GG_v1/Mcerasi_GG_assembly/Mc_wild_RNAseq.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits longestORF.blastp.uniref.DIAMOND.tab


#######################
# predicted cds vs uniref
diamond blastx -k 1 -p 8 --sensitive -v -q /p/path_to/Mcerasi_wild_samples/RNAseq_final_GG_v1/Mcerasi_GG_assembly/Mc_wild_RNAseq.fasta -d /p/path_to/Downloads/blast_databases/uniref90.dmnd -a Mc_wild.blastx.uniref.DIAMOND001.diamond
diamond view -a Mc_wild.blastx.uniref.DIAMOND001.*.daa -f tab -o Mc_wild.blastx.uniref.DIAMOND001.tab


diamond blastp -k 1 -p 8 --sensitive -v -q /p/path_to/Mcerasi_wild_samples/RNAseq_final_GG_v1/Mcerasi_GG_assembly/Mc_wild_RNAseq.fasta.transdecoder.pep -d /p/path_to/Downloads/blast_databases/uniref90.dmnd -a Mc_wild.blastp.uniref.DIAMOND001.diamond
diamond view -a Mc_wild.blastp.uniref.DIAMOND001.*.daa -f tab -o Mc_wild.blastp.uniref.DIAMOND001.tab

# Using a taxominoic database -old, but still useful. 
#diamond blastx -p 8 --sensitive -v -q /p/path_to/Mcerasi_wild_samples/RNAseq_final_GG_v1/Mcerasi_GG_assembly/Mc_wild_RNAseq.fasta -d /p/path_to/Downloads/uni_ref/uniref90.0.79.dmnd -a Mc_wild.blastx.uniref.tx_id.DIAMOND002.diamond
#diamond blastp -p 8 --sensitive -v -q /p/path_to/Mcerasi_wild_samples/RNAseq_final_GG_v1/Mcerasi_GG_assembly/Mc_wild_RNAseq.fasta.transdecoder.pep -d /p/path_to/Downloads/uni_ref/uniref90.0.79.dmnd -a Mc_wild.blastp.tx_id.uniref.DIAMOND002.diamond

# Vs NR
diamond blastp -p 8 --sensitive -v -q /p/path_to/Mcerasi_wild_samples/RNAseq_final_GG_v1/Mcerasi_GG_assembly/Mc_wild_RNAseq.fasta.transdecoder.pep -d ~/Downloads/blast_databases/nr.dmnd -a p.Mc_wild_vs_nr.diamond
diamond view -a p.Mc_wild_vs_nr*.daa -f tab -o Mc_wild.blastp.nr.tab
python ~/misc_python/diamond_blast_to_kingdom/Diamond_blast_to_taxid_add_kingdom_add_species_description.py -i Mc_wild.blastp.nr.tab -p /p/path_to/Downloads/blast_databases -o Mc_wild.blastp_nr_tax.tab


wait
mkdir for_trinnotate

########################################################################################################################################################
#				for trinnotate
#	pfam domains
hmmscan --cpu 8 --domtblout TrinotatePFAM.out /p/path_to/Downloads/TransDecoder_r20131117/pfam/Pfam-AB.hmm.bin /p/path_to/Mcerasi_wild_samples/RNAseq_final_GG_v1/Mcerasi_GG_assembly/Mc_wild_RNAseq.fasta.transdecoder.pep > pfam_pep.log


#		swissprot
#blastx -query /p/path_to/Mcerasi_wild_samples/RNAseq_final_GG_v1/Mcerasi_GG_assembly/Mc_wild_RNAseq.fasta-db /p/path_to/Downloads/Trinotate-2.0.2/uniprot_sprot.fasta -num_threads 8 -max_target_seqs 1 -outfmt 6 > blastx_trin.outfmt6

diamond blastx -k 1 -p 8 --sensitive -v -q /p/path_to/Mcerasi_wild_samples/RNAseq_final_GG_v1/Mcerasi_GG_assembly/Mc_wild_RNAseq.fasta -d /p/path_to/Downloads/blast_databases/uniprot.dmnd -a Mc_wild.blastx.unipot.DIAMOND002.diamond
diamond view -a Mc_wild.blastx.unipot.DIAMOND002.*.daa -f tab -o ./for_trinnotate/blastx_trin_uniprot.outfmt6


#blastp -query /p/path_to/Mcerasi_wild_samples/RNAseq_final_GG_v1/Mcerasi_GG_assembly/Mc_wild_RNAseq.fasta.transdecoder.pep -db /p/path_to/Downloads/Trinotate-2.0.2/uniprot_sprot.fasta -num_threads 8 -max_target_seqs 1 -outfmt 6 > blastp_trin.outfmt6
diamond blastp -k 1 -p 8 --sensitive -v -q /p/path_to/Mcerasi_wild_samples/RNAseq_final_GG_v1/Mcerasi_GG_assembly/Mc_wild_RNAseq.fasta -d /p/path_to/Downloads/blast_databases/uniprot.dmnd -a Mc_wild.blastp.unipot.DIAMOND002.diamond
diamond view -a Mc_wild.blastp.unipot.DIAMOND002.*.daa -f tab -o ./for_trinnotateblastp_trin_uniprot.outfmt6


#		uniref  (uniref90.dmnd)

#blastx -query /p/path_to/Mcerasi_wild_samples/RNAseq_final_GG_v1/Mcerasi_GG_assembly/Mc_wild_RNAseq.fasta -db /p/path_to/Downloads/Trinotate-2.0.2/uniprot_uniref90.trinotate.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 > uniref90.blastx.outfmt6
diamond blastx -k 1 -p 8 --sensitive -v -q /p/path_to/Mcerasi_wild_samples/RNAseq_final_GG_v1/Mcerasi_GG_assembly/Mc_wild_RNAseq.fasta -d /p/path_to/Downloads/blast_databases/uniref90.dmnd -a Mc_wild.blastx.uniref90.dmnd.DIAMOND003.diamond
diamond view -a Mc_wild.blastx.uniref90.dmnd.DIAMOND003.*.daa -f tab -o ./for_trinnotate/uniref90.blastx.outfmt6


#blastp -query /p/path_to/Mcerasi_wild_samples/RNAseq_final_GG_v1/Mcerasi_GG_assembly/Mc_wild_RNAseq.fasta.transdecoder.pep -db /p/path_to/Downloads/Trinotate-2.0.2/uniprot_uniref90.trinotate.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 > uniref90.blastp.outfmt6
diamond blastx -k 1 -p 8 --sensitive -v -q /p/path_to/Mcerasi_wild_samples/RNAseq_final_GG_v1/Mcerasi_GG_assembly/Mc_wild_RNAseq.fasta -d /p/path_to/Downloads/blast_databases/uniref90.dmnd -a Mc_wild.blastx.uniref90.dmnd.DIAMOND003.diamond
diamond view -a Mc_wild.blastx.uniref90.dmnd.DIAMOND003.*.daa -f tab -o ./for_trinnotate/uniref90.blastx.outfmt6

# 		generate the trans to genes map
/p/path_to/Downloads/trinityrnaseq_r20140717/util/support_scripts/get_Trinity_gene_to_trans_map.pl p/path_to/Mcerasi_wild_samples/RNAseq_final_GG_v1/Mcerasi_GG_assembly/Mc_wild_RNAseq.fasta > ./for_trinnotate/Trinity.fasta.gene_trans_map



#
#		Running signalP to predict signal peptides
signalp -f short -n ./for_trinnotate/signalp.out /p/path_to/Mcerasi_wild_samples/RNAseq_final_GG_v1/Mcerasi_GG_assembly/Mc_wild_RNAseq.fasta.transdecoder.pep 
#wait

#		Running tmHMM to predict transmembrane regions
tmhmm --short < /p/path_to/Mcerasi_wild_samples/RNAseq_final_GG_v1/Mcerasi_GG_assembly/Mc_wild_RNAseq.fasta.transdecoder.pep > ./for_trinnotate/tmhmm.out
#wait
#		Running RNAMMER to identify rRNA transcripts


###########################################################################################################################################################################
# BIOLINUX only

#/p/Downloads/Trinotate_r20131110/util/rnammer_support/RnammerTranscriptome.pl --transcriptome Mc_wild_RNAseq_ass.fa --path_to_rnammer /p/Downloads/rnammer

#wait




######################################################################################################################################################################################

#		TRINNOTATE

#wget "ftp://ftp.broadinstitute.org/pub/Trinity/Trinotate_v2.0_RESOURCES/Trinotate.sprot_uniref90.20150131.boilerplate.sqlite.gz" -O Trinotate.sqlite.gz



#gunzip Trinotate.sqlite.gz



# load protein hits
#Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp_trin.outfmt6

# load transcript hits
#Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx_trin.outfmt6

#Optional: load Uniref90 blast hits (requires the more comprehensive boilerplate database, as described above):

# load protein hits
#Trinotate Trinotate.sqlite LOAD_trembl_blastp uniref90.blastp.outfmt6

# load transcript hits
#Trinotate Trinotate.sqlite LOAD_trembl_blastx uniref90.blastx.outfmt6



#Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out



#Trinotate Trinotate.sqlite LOAD_tmhmm tmhmm.out



#Trinotate Trinotate.sqlite LOAD_signalp signalp.out


#Trinotate Trinotate.sqlite report -E 1-10 > Mc_wild_trinotate_annotation_report.xls


#########################################################################################################
##### shell script to annotate  using pfamAB domains swiss prot and uniref90	######################################
######################################################################################################

#wait
