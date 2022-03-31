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

# working opt variable
# parse input variables
option_list <- list(
  make_option(c("-t", "--threads"), type = "integer", default = parallel::detectCores(all.tests = FALSE, logical = TRUE),
              help = "number of threads to use (default is maximum available)", metavar = "integer"),
  make_option(c("-g", "--genome"), type = "character", default = NULL,
              help = "name of original genome", metavar = "character"),
  make_option(c("-f", "--flank"), type = "integer", default = 1500,
              help = "length of flanks t use", metavar = "integer"),
  make_option(c("-o", "--outgroups"), type = "character", default = NULL,
              help = "list of genomes to curate TEs from", metavar = "character")
)

message("Parsing variables")
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# check variables are set
if ( is.null(opt$genome) | is.null(opt$outgroups)) {
  stop("Necessary variables have not been set.")
}

genome_list <- readr::read_tsv(file = opt$outgroups, col_names = c("species_name", "genome_name"))

for( i in 1:nrow(genome_list) ){

  # Check file exists and has is not empty
  if(file.exists(paste0("data/", opt$genome, "_candiates_in_", genome_list$species_name[i], ".out"))){
    if(file.size(paste0("data/", opt$genome, "_candiates_in_", genome_list$species_name[i], ".out")) > 0){
      
      for_curation <- readr::read_tsv(file = paste0("data/", opt$genome, "_candiates_in_", genome_list$species_name[i], ".out"),
                                      col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                                    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"),
                                      num_threads = opt$threads, show_col_types = FALSE) %>%
        dplyr::group_by(qseqid) %>%
        dplyr::arrange(-bitscore) %>%
        dplyr::slice(1:30) %>%
        dplyr::ungroup()
      
      curation_list <- for_curation %>% dplyr::select(qseqid) %>% base::unique()
      
      genome_seq <- readDNAStringSet(paste0("genomes/", genome_list$genome_name[i], ".fasta"))
      names(genome_seq) <- sub(" .*", "", names(genome_seq))
      
      for(j in 1:nrow(curation_list)){
        
        for_curation_ranges <- for_curation %>%
          dplyr::filter(qseqid == curation_list$qseqid[j], length/qlen >= 0.3) %>%
          dplyr::mutate(start = ifelse(sstart < send, sstart - opt$flank, send - opt$flank),
                        end = ifelse(sstart > send, sstart + opt$flank, send + opt$flank),
                        strand = ifelse(sstart < send, "+", "-"),
                        start = ifelse(start < 1, 1, start),
                        end = ifelse(end > slen, slen, end)) %>%
          dplyr::select(sseqid, start, end, strand) %>%
          dplyr::rename(seqnames = sseqid) %>%
          plyranges::as_granges() %>%
          IRanges::reduce()
        
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
        
        Biostrings::writeXStringSet(trimmed_curation_seq, "out/temp.fa")
        
        # align
        system(paste0("mafft --thread ", opt$threads, " --localpair out/temp.fa > out/other_alignments/",
                      genome_list$species_name[i], "_", sub("/", "_", curation_list$qseqid[j]), ".fasta"))
        
      }
    }
  }
}
