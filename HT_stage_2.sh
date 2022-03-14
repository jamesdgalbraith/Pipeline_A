#!/bin/bash

# search genomes (initial sweep)
### split query into pieces for faster search
Rscript ~/projectDrive_2/DDE_Pipeline_2/splitter.R -f curated_${QUERY} -p ${THREADS} -t DNA
ls split/curated_${QUERY}_seq_* > split_seq_list.txt
sed 's/.*\///' split_seq_list.txt -i

### make genome databases
makeblastdb -in seq/${GENOME} -dbtype nucl
makeblastdb -in seq/${OUTGROUP} -dbtype nucl

### search genome in parallel. do for source and outgroup(s)
parallel --bar --jobs ${THREADS} -a split_seq_list.txt 'blastn -query split/{} -db seq/${GENOME} -out data/{}_${GENOME}.out -outfmt "6 std qlen slen"'
parallel --bar --jobs ${THREADS} -a split_seq_list.txt 'blastn -query split/{} -db seq/${OUTGROUP} -out data/{}_${OUTGROUP}.out -outfmt "6 std qlen slen"'

### compile blast out
cat data/curated_${QUERY}_seq_*${GENOME}.out > data/curated_${QUERY}_${GENOME}.out
cat data/curated_${QUERY}_seq_*${OUTGROUP}.out > data/curated_${QUERY}_${OUTGROUP}.out

### remove raw split blast
rm data/curated_*_seq_*.out

# Rscript to identify candidates
Rscript HT_filter_2.R -source ${GENOME} -outgroup ${OUTGROUP}