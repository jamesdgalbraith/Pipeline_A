#!/bin/bash
# usage: GENOME=<source_genome> OUTGROUPS=<file_containing_list_of_genomes> QUERY=<file_containing_curated_repeats> THREADS=<number of threads to use> bash HT_stage_2.sh 
# OUTGROUPS files should have species names in column 1 and genome names in column b

export GENOME

mkdir out/other_alignments

### make database of and search each genome, remove database
while read a b; do

  makeblastdb -in seq/$b -dbtype nucl -out seq/$a
  blastn -query out/final_HTT_${GENOME}_candiates.fasta -db seq/${a} -outfmt "6 std qlen slen" -task dc-megablast | \
    awk '{if ($4/$13 > 0.1)}' > data/${GENOME}_candiates_in_${a}.out
  rm seq/$b.n*

done<${OUTGROUPS}

# Rscript to identify candidates
Rscript HT_curator.R --genome ${GENOME} --outgroups ${OUTGROUPS} --threads ${THREADS}