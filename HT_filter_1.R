# read in packages
library("tidyverse")
library("plyranges")
library("BSgenome")

flank_len <- 1500
threads <- 64

# read in query seq
query_seq <- Biostrings::readDNAStringSet(filepath = "latCol_rm.fa.centroids")
query_tbl <- tibble(seqnames = base::names(query_seq), start = 0, end = BiocGenerics::width(query_seq))

# read in genome searches, remove small hits, satellites and simple repeats
self_blast_out <- readr::read_tsv(file = paste0("data/latCol_rm.fa.centroids_latCor_2.0.fasta.out"),
                                      col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                                    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen")) %>%
  filter(length >= 0.1 * qlen, !grepl("Simple_repeat", qseqid), !grepl("Satellite", qseqid))

outgroup_blast_out <- readr::read_tsv(file = paste0("data/latCol_rm.fa.centroids_TS10Xv2-PRI.fasta.out"),
                             col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                           "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen")) %>%
  filter(length >= 0.1 * qlen, !grepl("Simple_repeat", qseqid), !grepl("Satellite", qseqid))

# calculate average % identity of 10 best hits
self_blast_out_top_10 <- self_blast_out %>%
  dplyr::group_by(qseqid) %>%
  dplyr::arrange(-bitscore) %>%
  dplyr::slice(1:10) %>%
  dplyr::mutate(av_pident = base::mean(pident)) %>%
  dplyr::ungroup()

outgroup_blast_out_top_10 <- outgroup_blast_out %>%
  dplyr::group_by(qseqid) %>%
  dplyr::arrange(-bitscore) %>%
  dplyr::slice(1:10) %>%
  dplyr::mutate(av_pident = base::mean(pident)) %>%
  dplyr::ungroup()

# join top 10 self blast and outgroup searches to reveal missing repeats
search_joined <- inner_join(x = (self_blast_out_top_10 %>%
                  dplyr::select(qseqid, av_pident) %>%
                  dplyr::rename(self_pident = av_pident) %>%
                  base::unique()),
           y = (outgroup_blast_out_top_10 %>%
                  dplyr::select(qseqid, av_pident) %>%
                  dplyr::rename(outgroup_pident = av_pident) %>%
                  base::unique()))

# identify potentially too divergent
pident_to_low <- search_joined %>%
  filter(self_pident > 90, outgroup_pident < 75)

# select too divergent and missing from top 10
ht_candidates_top_10 <- self_blast_out_top_10 %>%
  filter((!qseqid %in% search_joined$qseqid | qseqid %in% pident_to_low$qseqid), av_pident >= 90)

# ensure at least 2 copies present to avoid anomolies
ht_candidates_tbl <- dplyr::as_tibble(base::as.data.frame(base::table(ht_candidates_top_10$qseqid))) %>%
  dplyr::mutate(Var1 = as.character(Var1)) %>%
  dplyr::rename(qseqid = Var1, n = Freq) %>%
  dplyr::filter(n > 1)

### STEP TO KILL IF NO CANDIDATES ###
if(nrow(ht_candidates_tbl) < 1){
  
}
  

# get sequence of repeats
ht_candidate_seq <- query_seq[sub(" .*", "", names(query_seq)) %in% ht_candidates_tbl$qseqid]

# prepare for manual curation
for_curation_out <- readr::read_tsv(file = paste0("data/latCol_rm.fa.centroids_latCor_2.0.fasta.out"),
                                  col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                                "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen")) %>%
  filter(length >= 0.1 * qlen, qseqid %in% ht_candidates_tbl$qseqid, pident >= 90) %>%
  dplyr::group_by(qseqid) %>%
  dplyr::arrange(-bitscore) %>%
  dplyr::slice(1:30) %>%
  dplyr::ungroup()

# read in genome
genome_seq <- readDNAStringSet("seq/latCor_2.0.fasta")
names(genome_seq) <- sub(" .*", "", names(genome_seq))

for(i in 1:nrow(ht_candidates_tbl)){

  for_curation_ranges <- for_curation_out %>%
    dplyr::filter(qseqid == ht_candidates_tbl$qseqid[i]) %>%
    dplyr::mutate(start = ifelse(sstart < send, sstart - flank_len, send - flank_len),
                  end = ifelse(sstart > send, sstart + flank_len, send + flank_len),
                  strand = ifelse(sstart < send, "+", "-"),
                  start = ifelse(start < 1, 1, start),
                  end = ifelse(end > slen, slen, end)) %>%
    dplyr::select(sseqid, start, end, strand) %>%
    dplyr::rename(seqnames = sseqid) %>%
    plyranges::as_granges()
  
  # get sequence for alignment and name
  for_curation_seq <- getSeq(genome_seq, for_curation_ranges)
  names(for_curation_seq) <- paste0(seqnames(for_curation_ranges), ":", ranges(for_curation_ranges), "(", strand(for_curation_ranges), ")")
  
  # add original consensus
  for_curation_seq <- c(query_seq[sub(" .*", "", names(query_seq)) == ht_candidates_tbl$qseqid[i]], for_curation_seq)
  
  # write to file
  Biostrings::writeXStringSet(for_curation_seq, "out/temp.fa")
  
  # align
  system(paste0("mafft --thread ", threads, " --localpair out/temp.fa > out/aligned/", sub("/", "_", ht_candidates_tbl$qseqid[i]), ".fasta"))
}