#!/bin/bash
# usage: GENOME=<source_genome> OUTGROUPS=<file_containing_list_of_genomes> QUERY=<file_containing_curated_repeats> THREADS=<number of threads to use> bash HT_stage_2.sh 
# OUTGROUPS files should have species names in column 1 and genome names in column b

export GENOME
export OUTGROUPS

mkdir -p out/other_alignments

# unzip fastas and make genome databases
parallel --jobs 4 -a ${OUTGROUPS} --colsep="\t" 'gunzip < genomes/{2}_genomic.fna.gz > genomes/{2}.fasta'
parallel --jobs 4 -a ${OUTGROUPS} --colsep="\t" 'makeblastdb -in genomes/{2}.fasta -dbtype nucl'

# search genomes for repeats
parallel --jobs 8 -a ${OUTGROUPS} --colsep="\t" 'echo blastn -query out/final_HTT_${GENOME}_candiates.fasta -db genomes/{2}.fasta -outfmt "6 std qlen slen" -task dc-megablast -out data/${GENOME}_candiates_in_{1}.out'

# remove database files
rm genomes/*fasta*n*

# Rscript to identify candidates
Rscript HT_curator.R --genome ${GENOME} --outgroups ${OUTGROUPS} --threads ${THREADS}
