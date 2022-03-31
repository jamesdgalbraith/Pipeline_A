#!/bin/bash
# usage: SPECIES=<name_of_source_species> GENOME=<source_genome> OUTGROUP=<outgroup_genome> QUERY=<file_containing_curated_repeats> THREADS=<number of threads to use> bash HT_stage_2.sh 

export GENOME
export OUTGROUP

# search genomes (initial sweep)
### split query into pieces for faster search
Rscript splitter.R -f ${QUERY} -p ${THREADS} -t DNA
ls split/${QUERY}_seq_* | sed 's/.*\///' > split_seq_list.txt

### search genome in parallel. do for source and outgroup(s)
parallel --bar --jobs ${THREADS} -a split_seq_list.txt 'blastn -query split/{} -db genomes/${GENOME} -out data/{}_${GENOME}.out -outfmt "6 std qlen slen" -task dc-megablast'
parallel --bar --jobs ${THREADS} -a split_seq_list.txt 'blastn -query split/{} -db genomes/${OUTGROUP} -out data/{}_${OUTGROUP}.out -outfmt "6 std qlen slen" -task dc-megablast'

### compile blast out
cat data/${QUERY}_seq_*${GENOME}.out > data/${QUERY}_${GENOME}.out
cat data/${QUERY}_seq_*${OUTGROUP}.out > data/${QUERY}_${OUTGROUP}.out

### remove raw split blast
rm data/${QUERY}_seq_*.out split/${QUERY}_seq_*

# Rscript to identify candidates
Rscript HT_filter_2.R --genome ${GENOME} --outgroup ${OUTGROUP} --query ${QUERY} --species ${SPECIES}