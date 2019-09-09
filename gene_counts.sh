#mc wild

cd /home/pt40963/M.cerasi_genome_data/final_assembly_v1/wild_samples_DE


bedtools multicov -bams ./Mc_cherry1/sorted.bam ./Mc_cherry2/sorted.bam ./Mc_cherry3/sorted.bam ./Mc_cress1/sorted.bam ./Mc_cress2/sorted.bam ./Mc_cress3/sorted.bam ./Mc_Gallium1/sorted.bam ./Mc_Gallium2/sorted.bam ./Mc_Gallium3/sorted.bam -bed /home/pt40963/Aphid_gene_models/upstream/Mc/format_for_py_script.gff >  ./Mc_cherry123_cress_123_gallium_123.counts.gff
bedtools multicov -bams ./Mc_cherry1/sorted.bam ./Mc_cherry2/sorted.bam ./Mc_cherry3/sorted.bam ./Mc_cress1/sorted.bam ./Mc_cress2/sorted.bam ./Mc_cress3/sorted.bam ./Mc_Gallium1/sorted.bam ./Mc_Gallium2/sorted.bam ./Mc_Gallium3/sorted.bam -bed /home/pt40963/Aphid_gene_models/upstream/Mc/format_for_py_script.gff >  ./Mc_cherry123_cress_123_gallium_123.counts.gff
cat Mc_cherry123_cress_123_gallium_123.counts.gff | cut -f5,6,7,8,9,10,11,12,13,14 > Mc_cherry_cress_galium123.genes.counts.matrix
cat Mc_cherry123_cress_123_gallium_123.counts.gff | cut -f5,9,10,11,12,13,14 > Mc_cress_galium123.genes.counts.matrix