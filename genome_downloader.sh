#!/bin/bash

# make directory for genomes
mkdir -p genomes

# get genome information
esearch -db assembly -query '("Serpentes"[Organism] OR "Echinodermata"[Organism]) AND ("latest genbank"[filter] AND ("representative genome"[filter]  OR "representative genome"[filter]) AND all[filter] NOT anomalous[filter] AND ("100000"[ScaffoldN50] : "5000000000"[ScaffoldN50])' | \
esummary | xtract -pattern DocumentSummary -def "NA" -element Genbank,AssemblyName,Organism,FtpPath_GenBank | \
sed 's/ /_/g;s/_(.*)//' | awk '{if ($4 ~ "ftp") print $0 "/"$1"_"$2"_genomic.fna.gz"}' | sed 's|ftp://|rsync://|' > updated_genome_info.txt

# make blank document for new genome paths
> new_genome_info.txt

# check if genomes already downloaded
while read a b c d; do
FILE=genomes/${a}_${b}_genomic.fna.gz
if [ ! -f "$FILE" ]; then
    echo ${d} >> new_genome_info.txt
fi
done<updated_genome_info.txt

# download genomes not already downloaded
parallel -vv --jobs 16 -a new_genome_info.txt rsync --copy-links --recursive --times --verbose {} genomes/

# make species genome names txt
awk '{OFS="\t"}{print $3,$1"_"$2}' updated_genome_info.txt > species_genome_list.txt