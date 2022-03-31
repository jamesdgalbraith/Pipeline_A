# read in packages
base::suppressPackageStartupMessages(library("tidyverse"))
base::suppressPackageStartupMessages(base::library("plyranges"))
base::suppressPackageStartupMessages(base::library("BSgenome"))
base::suppressPackageStartupMessages(base::library("optparse"))

# parse input variables
option_list <- list(
  make_option(c("-q", "--query"), type = "character", default = NULL,
              help = "repeat files", metavar = "character"),
  make_option(c("-g", "--genome"), type = "character", default = "NULL",
              help = "genome file name", metavar = "character"),
  make_option(c("-o", "--outgroup"), type = "character", default = "NULL",
              help = "outgroup genome file name", metavar = "character"),
  make_option(c("-s", "--species"), type = "character", default = "NULL",
              help = "name of query species for final output", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

query <- opt$query
genome_name <- opt$genome
outgroup_name <- opt$outgroup

candidate_seq <- readDNAStringSet(query)

# read in genome searches, remove small hits, satellites and simple repeats
self_blast_out <- readr::read_tsv(file = paste0("data/", query, "_", genome_name, ".out"),
                                  col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                                "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"),
                                  show_col_types = F) %>%
  filter(length >= 0.1 * qlen, !grepl("Simple_repeat", qseqid), !grepl("Satellite", qseqid))

outgroup_blast_out <- readr::read_tsv(file = paste0("data/", query, "_", outgroup_name, ".out"),
                                      col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                                    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"),
                                      show_col_types = F) %>%
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

# step to kill if no candidates
if (nrow(pident_to_low) == 0 & nrow(search_joined) == 0) {
  stop("No HTT candidates found")
}

# select too divergent and missing from top 10
ht_candidates_top_10 <- self_blast_out_top_10 %>%
  filter((!qseqid %in% search_joined$qseqid | qseqid %in% pident_to_low$qseqid), av_pident >= 90)

# ensure at least 2 copies present to avoid anomalies
ht_candidates_tbl <- dplyr::as_tibble(base::as.data.frame(base::table(ht_candidates_top_10$qseqid))) %>%
  dplyr::mutate(Var1 = as.character(Var1)) %>%
  dplyr::rename(qseqid = Var1, n = Freq) %>%
  dplyr::filter(n > 1)

# select true candidates
actual_candidate_seq <- candidate_seq[sub(" .*", "", names(candidate_seq)) %in% ht_candidates_tbl$qseqid]
writeXStringSet(x = actual_candidate_seq, filepath = paste0("out/final_HTT_", opt$species, "_candiates.fasta"))