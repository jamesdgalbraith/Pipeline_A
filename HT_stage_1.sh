#!/bin/bash
# usage: GENOME=<source_genome> OUTGROUP=<outgroup_genome> QUERY=<file_containing_repeats> THREADS=<number of threads to use> bash HT_stage_1.sh

export GENOME
export OUTGROUP

mkdir -p data out/aligned seq genomes/
makeblastdb -in genomes/{GENOME} -dbtype nucl
makeblastdb -in genomes/{OUTGROUP} -dbtype nucl

# search genomes (initial sweep)
### cluster query
vsearch --cluster_size ${QUERY} --id 0.8 --centroids data/${QUERY}.centroids

### split query into pieces for faster search
Rscript splitter.R -f data/${QUERY}.centroids -p ${THREADS} -t DNA
ls split/${QUERY}.centroids_seq_* | sed 's/.*\///' > split_seq_list.txt

### search genome in parallel. do for source and outgroup(s)
parallel --bar --jobs ${THREADS} -a split_seq_list.txt 'blastn -query split/{} -db genomes/${GENOME} -out data/{}_${GENOME}.out -outfmt "6 std qlen slen" -task dc-megablast'
parallel --bar --jobs ${THREADS} -a split_seq_list.txt 'blastn -query split/{} -db genomes/${OUTGROUP} -out data/{}_${OUTGROUP}.out -outfmt "6 std qlen slen" -task dc-megablast'

### compile blast out
cat data/${QUERY}.centroids_seq_*${GENOME}.out > data/${QUERY}.centroids_${GENOME}.out
cat data/${QUERY}.centroids_seq_*${OUTGROUP}.out > data/${QUERY}.centroids_${OUTGROUP}.out

### remove raw split blast and working files
rm data/${QUERY}.centroids_seq_*_${GENOME}.out data/${QUERY}.centroids_seq_*_${OUTGROUP}.out split_seq_list.txt split/*

### make file of PATH for R to use
echo $PATH > path.txt

# Rscript to identify candidates and curate
Rscript HT_filter_1.R --query ${QUERY} --genome ${GENOME} --outgroup ${OUTGROUP} --threads ${THREADS}

# perform manual curation