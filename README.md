Methods used to investigate transcriptional responses in Myzus cerasi on different hosts
========================================================================================

The raw RNAseq reads for M. cerasi are available at study accession: PRJEB24338
Myzus cerasi genome and annotation was downloaded from http://bipaa.genouest.org/is/aphidbase/ and DOI:10.5281/zenodo.1252934.

1) Quality trim the reads and fastqc
``./read_trimming_Q22.sh

2) Once the reads were "cleaned", the final reads were converted mapped to the genome using:
``./Myzus.cerasi_map_to_genome_STAR.sh``

3) De novo, genome guided RNAseq assembly and post annotation was performed using:
``Myzus.cerasi_genome_guided_trinity_assembly.sh``

4) Transrate was used to reduce the RNAseq assembly to well supported contigs:
``transrate.sh``

5) Extra assembly QC:
``Mc_RNAQUAST.sh``

Blast to assembly against NR
``Mc_RNAseq_vs_nr_xml.sh``

6) Exon counts was perfomed using:
``Mc_exon_counts.sh``

7) Differential exon expression was perfomed using:
``differential_exon_R_commands.sh``

8) Map RNAseq to genome
``Myzus.cerasi_map_to_genome_STAR.sh``

9) count were generated:
``gene_counts.sh``

10) DE was perfomed using:
``Diff_Expression_trinity_trinity2.4_gene_models_order.sh``

GO enrichment was perfomed in Blast2Go

11) gene duplciation classes were obtained and counted using:
``duplication_quantify.py``

12) random gene duplcation types were obtained and quantified using:
``get_random_genes_duplication_types.py``

13) SNP/ variant analyses
``SNPs_Myzus_cerasi.sh``

14) compare genetic diversity from 13
``compare_genetic_diversity.py``
