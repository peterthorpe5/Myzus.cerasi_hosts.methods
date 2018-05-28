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
