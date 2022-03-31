#!/bin/bash
# usage: GENOME=<source_genome> THREADS=<number of threads to use> bash HT_stage_0.sh 

# run RepeatModeler
### create database
BuildDatabase -name genomes/${GENOME} genomes/${GENOME}

### run repeatmodeler
RepeatModeler -pa ${THREADS} -database genomes/${GENOME}

### remove database
rm genomes/${GENOME}.n* genomes/${GENOME}.translation