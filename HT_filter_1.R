# read in packages
message("Loading libraries")
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("plyranges"))
suppressPackageStartupMessages(library("BSgenome"))
suppressPackageStartupMessages(library("optparse"))

# set R sh path to be the same as computers (ensures access to mafft and blastn)
con <- file("path.txt")
Sys.setenv(PATH = as.character(readLines(con)))
close(con)
# parse input variables
option_list <- list(
  make_option(c("-q", "--query"), type = "character", default = NULL,
              help = "repeat files", metavar = "character"),
  make_option(c("-t", "--threads"), type = "integer", default = parallel::detectCores(all.tests = FALSE, logical = TRUE),
              help = "number of threads to use (default is maximum available)", metavar = "integer"),
  make_option(c("-g", "--genome"), type = "character", default = NULL,
              help = "genome file name", metavar = "character"),
  make_option(c("-f", "--flank"), type = "integer", default = 3000,
              help = "outgroup genome file name", metavar = "integer"),
  make_option(c("-o", "--outgroup"), type = "character", default = NULL,
              help = "outgroup genome file name", metavar = "character")
)

message("Parsing variables")
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# check variables are set
if (is.null(opt$query) | is.null(opt$genome) | is.null(opt$outgroup)) {
  stop("Necessary variables have not been set.")
}

# check files exist
if (!file.exists(paste0("data/", opt$query, ".centroids")) |
    !file.exists(paste0("data/", opt$query, ".centroids_", opt$genome, ".out")) |
    !file.exists(paste0("data/", opt$query, ".centroids_", opt$outgroup, ".out"))
    ) {
  stop("Input file(s) are missing, ensure you have run the previous steps correctly")
}

# set variables
query <- opt$query
genome_name <- opt$genome
outgroup_name <- opt$outgroup
threads <- opt$threads
flank_len <- opt$flank

message("Reading data")
# read in query seq
query_seq <- Biostrings::readDNAStringSet(filepath = paste0("data/",query, ".centroids"))
query_tbl <- tibble(seqnames = base::names(query_seq), start = 0, end = BiocGenerics::width(query_seq))

# read in genome searches, remove small hits, satellites and simple repeats
self_blast_out <- readr::read_tsv(file = paste0("data/", query, ".centroids_", genome_name, ".out"),
                                      col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                                    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"),
                                  num_threads = threads, show_col_types = FALSE) %>%
  filter(length >= 0.1 * qlen, !grepl("Simple_repeat", qseqid), !grepl("Satellite", qseqid))

outgroup_blast_out <- readr::read_tsv(file = paste0("data/", query, ".centroids_", outgroup_name, ".out"),
                             col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                           "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"),
                             num_threads = threads, show_col_types = FALSE) %>%
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
search_joined <- suppressMessages(inner_join(x = (self_blast_out_top_10 %>%
                                                    dplyr::select(qseqid, av_pident) %>%
                                                    dplyr::rename(self_pident = av_pident) %>%
                                                    base::unique()),
                                             y = (outgroup_blast_out_top_10 %>%
                                                    dplyr::select(qseqid, av_pident) %>%
                                                    dplyr::rename(outgroup_pident = av_pident) %>%
                                                    base::unique())))

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
if (nrow(ht_candidates_tbl) < 1) {
  stop("No HT candidates found")
}


# Prepare MSA for curation
message("Preparing MSA for curation")

# get sequence of repeats
ht_candidate_seq <- query_seq[sub(" .*", "", names(query_seq)) %in% ht_candidates_tbl$qseqid]

# prepare for manual curation
for_curation_out <- readr::read_tsv(file = paste0("data/", query, ".centroids_", genome_name, ".out"),
                                  col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                                "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"),
                                  num_threads = threads, show_col_types = FALSE) %>%
  filter(length >= 0.5 * qlen, qseqid %in% ht_candidates_tbl$qseqid, pident >= 90) %>%
  dplyr::group_by(qseqid) %>%
  dplyr::arrange(-bitscore) %>%
  dplyr::slice(1:30) %>%
  dplyr::ungroup()

# read in genome
genome_seq <- readDNAStringSet(paste0("genomes/", genome_name))
names(genome_seq) <- sub(" .*", "", names(genome_seq))

for (i in 1:nrow(ht_candidates_tbl)) {

  message(paste0("Preparing ", i, " of ", nrow(ht_candidates_tbl)))
  for_curation_ranges <- for_curation_out %>%
    dplyr::filter(qseqid == ht_candidates_tbl$qseqid[i]) %>%
    dplyr::mutate(start = ifelse(sstart < send, sstart - flank_len, send - flank_len),
                  end = ifelse(sstart > send, sstart + flank_len, send + flank_len),
                  strand = ifelse(sstart < send, "+", "-"),
                  start = ifelse(start < 1, 1, start),
                  end = ifelse(end > slen, slen, end)) %>%
    dplyr::select(sseqid, start, end, strand) %>%
    dplyr::rename(seqnames = sseqid) %>%
    plyranges::as_granges() %>%
    IRanges::reduce()

  # get sequence for alignment and name
  for_curation_seq <- getSeq(genome_seq, for_curation_ranges)
  names(for_curation_seq) <- paste0(seqnames(for_curation_ranges), ":", ranges(for_curation_ranges), "(", strand(for_curation_ranges), ")")

  # blast for initial trim
  writeXStringSet(for_curation_seq, "out/temp.fa")
  system("blastn -query out/temp.fa -subject out/temp.fa -task dc-megablast -out out/temp.out -outfmt \"6 std qlen slen\"")
  recip_blast <- read_tsv("out/temp.out",
                          col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"),
                          show_col_types = FALSE)
  recip_blast_ranges <- recip_blast %>%
    filter(qseqid != sseqid) %>%
    dplyr::group_by(qseqid) %>%
    dplyr::arrange(-bitscore) %>%
    dplyr::slice(1:5) %>%
    dplyr::mutate(start = min(qstart), end = max(qend)) %>%
    dplyr::ungroup() %>%
    dplyr::select(qseqid, start, end) %>%
    dplyr::rename(seqnames = qseqid) %>%
    base::unique() %>%
    plyranges::as_granges()

  trimmed_curation_seq <- getSeq(for_curation_seq, recip_blast_ranges)
  names(trimmed_curation_seq) <- paste0(seqnames(recip_blast_ranges), ":", ranges(recip_blast_ranges))

  # add original consensus
  trimmed_curation_seq <- c(query_seq[sub(" .*", "", names(query_seq)) == ht_candidates_tbl$qseqid[i]], trimmed_curation_seq)

  # write to file
  Biostrings::writeXStringSet(trimmed_curation_seq, "out/temp.fa")

  # align
  system(paste0("mafft --thread ", threads, " --localpair out/temp.fa > out/aligned/", sub("/", "_", ht_candidates_tbl$qseqid[i]), ".fasta"))

}

message("All done!")