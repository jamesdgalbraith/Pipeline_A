#!/bin/bash

mkdir -p data out/aligned

# search genomes (initial sweep)
### cluster query
vsearch --cluster_size ${QUERY} --id 0.8 --centroids ${QUERY}.centroids

### split query into pieces for faster search
Rscript splitter.R -f ${QUERY}.centroids -p ${THREADS} -t DNA
ls split/${QUERY}.centroids_seq_* > split_seq_list.txt
sed 's/.*\///' split_seq_list.txt -i

### search genome in parallel. do for source and outgroup(s)
parallel --bar --jobs ${THREADS} -a split_seq_list.txt 'blastn -query split/{} -db seq/${GENOME} -out data/{}_${GENOME}.out -outfmt "6 std qlen slen"'
parallel --bar --jobs ${THREADS} -a split_seq_list.txt 'blastn -query split/{} -db seq/${OUTGROUP} -out data/{}_${OUTGROUP}.out -outfmt "6 std qlen slen"'

### compile blast out
cat data/${QUERY}.centroids_seq_*${GENOME}.out > data/${QUERY}.centroids_${GENOME}.out
cat data/${QUERY}.centroids_seq_*${OUTGROUP}.out > data/${QUERY}.centroids_${OUTGROUP}.out

### remove raw split blast
rm data/latCol_rm.fa.centroids_seq_*

# Rscript to identify candidates and curate
Rscript HT_filter_1.R -source ${GENOME} -outgroup ${OUTGROUP}

# perform manual curation