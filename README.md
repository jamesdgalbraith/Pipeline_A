## Pipeline A

This simple pipeline is designed to identify potential horizontally transferred sequences in newly sequenced genomes based on their absence from the genome of a closely related species.

Requirements: MAFFT, NCBI BLAST, R with tidyverse, plyranges, BSgenome and optranges packages installed

All input files should be FASTA files. Raw repeats can be from any *ab initio* software package.

Usage:

1\) Run HT_stage_1.sh

`GENOME=<source_genome> OUTGROUP=<outgroup_genome> QUERY=<file_containing_raw_repeats> THREADS=<number of threads to use> bash HT_stage_1.sh`

This script performs an initial search for sequences which have 2 or more copies in the query genome and are absent from the outgroup genome. If any sequences are absent a multiple sequence alignment will be created for curation.

2\) Manually curate the alignments in Geneious, JalView or the like. This is necessary as most *ab initio* aligners do not capture full repeats, and some families may be classified incorrectly. Additionally, this step can reveal redundant sequences not removed from clustering, e.g. non-autonomous DNA transposons derived from autonomous DNA transposons.

3\) Run HT_stage_2.sh with using curated repeats

`GENOME=<source_genome> OUTGROUP=<outgroup_genome> QUERY=<file_containing_curated_repeats> THREADS=<number of threads to use> bash HT_stage_2.sh`
