#!/bin/bash

GENOME="latCor_2.0.fasta"
export GENOME
OUTGROUP="TS10Xv2-PRI.fasta"
export OUTGROUP

# search genomes (initial sweep)
### split query into pieces for faster search
Rscript splitter.R -f curated_${QUERY} -p ${THREADS} -t DNA
ls split/curated_${QUERY}_seq_* | sed 's/.*\///' > split_seq_list.txt

### search genome in parallel. do for source and outgroup(s)
parallel --bar --jobs ${THREADS} -a split_seq_list.txt 'blastn -query split/{} -db seq/${GENOME} -out data/{}_${GENOME}.out -outfmt "6 std qlen slen"'
parallel --bar --jobs ${THREADS} -a split_seq_list.txt 'blastn -query split/{} -db seq/${OUTGROUP} -out data/{}_${OUTGROUP}.out -outfmt "6 std qlen slen"'

### compile blast out
cat data/curated_${QUERY}_seq_*${GENOME}.out > data/curated_${QUERY}_${GENOME}.out
cat data/curated_${QUERY}_seq_*${OUTGROUP}.out > data/curated_${QUERY}_${OUTGROUP}.out

### remove raw split blast
rm data/curated_*_seq_*.out

# Rscript to identify candidates
Rscript HT_filter_2.R -genome ${GENOME} -outgroup ${OUTGROUP} -query curated_latCol_rm.fasta