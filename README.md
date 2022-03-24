## Pipeline A

This simple pipeline is designed to identify potential horizontally transferred sequences in newly sequenced genomes based on their absence from the genome of a closely related species.

Requirements: MAFFT, NCBI BLAST, R with tidyverse, plyranges, BSgenome and optranges packages installed

All input files should be FASTA files. Raw repeats can be from any *ab initio* software package.

Usage:

0\) *Ab initio* TE annotation using RepeatModeller. Run HT_stage_0.sh

```bash
GENOME=<source_genome> THREADS=<number of threads to use> bash HT_stage_0.sh 
```

1\) Run HT_stage_1.sh

```bash
GENOME=<source_genome> OUTGROUP=<outgroup_genome> QUERY=<file_containing_raw_repeats> THREADS=<number of threads to use> bash HT_stage_1.sh
```

This script performs an initial search for sequences which have 2 or more copies in the query genome and are absent from the outgroup genome. If any sequences are absent a multiple sequence alignment will be created for curation.

2\) Manually curate the alignments in Geneious, JalView or the like. This is necessary as most *ab initio* aligners do not capture full repeats, and some families may be classified incorrectly. Additionally, this step can reveal redundant sequences not removed from clustering, e.g. non-autonomous DNA transposons derived from autonomous DNA transposons.

3\) Run HT_stage_2.sh with using curated repeats

```bash
GENOME=<source_genome> OUTGROUP=<outgroup_genome> QUERY=<file_containing_curated_repeats> THREADS=<number of threads to use> bash HT_stage_2.sh
```

This script carries out a validation of the initial search using consensus sequences generated from the curation step. This is necessary as fragmented repeats which appeared to be mssing from an outgroup species may in fact be present. For example, when searching with TEs from a seal genome using a mustelid as the outgroup, stage 1 identified 4 L1 fragments which appeared to be absent from the mink. After curation it became clear those L1s were in fact present in the mustelid, just not identifed in the initial sweep, likely due to their truncation/fragmentation. The initial curation step fixes this problem by ensuring that searches of the outgroup genome are done with queries that contain the complete TE of interest.
