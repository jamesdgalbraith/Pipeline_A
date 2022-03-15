#!/bin/bash

mkdir -p data out/aligned

# search genomes (initial sweep)
### cluster query
vsearch --cluster_size ${QUERY} --id 0.8 --centroids ${QUERY}.centroids

### split query into pieces for faster search
Rscript splitter.R -f ${QUERY}.centroids -p ${THREADS} -t DNA
ls split/${QUERY}.centroids_seq_* | sed 's/.*\///' > split_seq_list.txt

### search genome in parallel. do for source and outgroup(s)
parallel --bar --jobs ${THREADS} -a split_seq_list.txt 'blastn -query split/{} -db seq/${GENOME} -out data/{}_${GENOME}.out -outfmt "6 std qlen slen"'
parallel --bar --jobs ${THREADS} -a split_seq_list.txt 'blastn -query split/{} -db seq/${OUTGROUP} -out data/{}_${OUTGROUP}.out -outfmt "6 std qlen slen"'

### compile blast out
cat data/${QUERY}.centroids_seq_*${GENOME}.out > data/${QUERY}.centroids_${GENOME}.out
cat data/${QUERY}.centroids_seq_*${OUTGROUP}.out > data/${QUERY}.centroids_${OUTGROUP}.out

### remove raw split blast and working files
rm data/${QUERY}.centroids_seq_*_${GENOME}.out data/${QUERY}.centroids_seq_*_${OUTGROUP}.out split_seq_list.txt split/*

### make file of PATH for R to use
echo $PATH > path.txt

# Rscript to identify candidates and curate
Rscript HT_filter_1.R --query ${QUERY} --genome ${GENOME} -outgroup ${OUTGROUP} --threads ${THREADS}

# perform manual curation