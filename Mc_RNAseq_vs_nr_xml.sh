#$ -l hostname="n13*"

cd /home/pt40963/scratch/M.cerasi/split_fasta_files

export BLASTDB=/mnt/shared/cluster/blast/ncbi/extracted

filenames=*.fasta


for f in ${filenames}
do
	cmd="blastx -db nr -query ${f} -evalue 0.00001 -seg no -num_threads 16 -outfmt 5 -out ${f}.xml" 
	echo ${cmd}
	eval ${cmd}
	wait
done



