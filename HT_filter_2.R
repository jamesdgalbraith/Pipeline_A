# read in packages
library("tidyverse")
library("plyranges")
library("BSgenome")

# read in genome searches, remove small hits, satellites and simple repeats
self_blast_out <- readr::read_tsv(file = paste0("data/latCol_rm.fa.curated_latCor_2.0.fasta.out"),
                                  col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                                "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen")) %>%
  filter(length >= 0.1 * qlen, !grepl("Simple_repeat", qseqid), !grepl("Satellite", qseqid))

outgroup_blast_out <- readr::read_tsv(file = paste0("data/latCol_rm.fa.curated_TS10Xv2-PRI.fasta.out"),
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