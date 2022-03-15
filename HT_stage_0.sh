#!/bin/bash
# usage: GENOME=<source_genome> THREADS=<number of threads to use> bash HT_stage_0.sh 

# run RepeatModeler
### create database
BuildDatabase -name seq/${GENOME} seq/${GENOME}

### run repeatmodeler
RepeatModeler -pa ${THREADS} -database seq/${GENOME}

### remove database
rm seq/${GENOME}.n* seq/${GENOME}.translation