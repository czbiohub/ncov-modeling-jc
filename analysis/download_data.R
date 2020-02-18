library(rentrez)
library(reutils)
library(seqinr)
library(stringi)
library(stringr)
library(readxl)
library(parsedate)
library(magrittr)
library(dplyr)

curr_date <- gsub("-", "", as.character(Sys.Date()))



# Download sequences from Genbank -----------------------------------------

coronavirus_taxid <- 2697049

id_search <- entrez_search(db="nuccore", term=paste0('txid', coronavirus_taxid, '[Organism:noexp] AND ("20000"[SLEN] : "35000"[SLEN])'))$ids
seq_str <- entrez_fetch(db="nuccore", id=id_search, rettype="fasta") %>% 
  strsplit("\n") %>%
  `[[`(1)
seq_name_pos <- grep(">", seq_str)
seq_str[seq_name_pos] <- word(seq_str[seq_name_pos])
cat(seq_str, sep="\n", file=paste0("../data/gb_", curr_date, "_nCov_genomes.fasta"))


# Download metadata associated with Genbank sequences ---------------------

gb <- efetch(id_search, db="nuccore", rettype="gb", retmode="text")$content
gb_split <- strsplit(gb, "LOCUS       ")[[1]][-1]
gb_line_split <- strsplit(gb_split, "\n")

acc <- gsub("([A-Za-z0-9]+).*", "\\1", gb_split)

gb_country <- gb_line_split %>% 
  mapply(grep, "/country", ., value=TRUE, SIMPLIFY=FALSE) %>%
  lapply(stri_extract_all_regex, '(?<=").*?(?=")') %>%
  lapply(unlist) %>%
  lapply(function (x) ifelse(length(x)>0, x, NA))

gb_coll_date_raw <- gb_line_split %>% 
  mapply(grep, "/collection_date", ., value=TRUE, SIMPLIFY=FALSE) %>%
  lapply(stri_extract_all_regex, '(?<=").*?(?=")') %>%
  lapply(unlist) %>%
  setNames(., acc)

gb_coll_date <- gb_coll_date_raw %>%
  lapply(parse_date) %>%
  lapply(function (x) ifelse(length(x)>0, as.character(x), NA)) %>%
  unlist() %>%
  as.Date()

gb_coll_date_month_only <- gb_coll_date_raw %>% 
  lapply(str_count, "-") %>% 
  sapply(function (x) ifelse(length(x)==0, NA, x==1))


metadata <- lapply(1:length(acc), function (i) {
  data.frame(acc=acc[i], country=gb_country[[i]], 
             collection_date=gb_coll_date[[i]], 
             collection_date_month_only=gb_coll_date_month_only[[i]],
             stringsAsFactors=TRUE)
}) %>%
  do.call(what=rbind, .)

write.table(metadata, paste0("../data/gb_", curr_date, "_metadata.tsv"), sep="\t", row.names = FALSE, quote = FALSE)


# Read in sequences from GISAID -------------------------------------------

gisaid_seq <- list.files("../data/GISAID", full.names=TRUE) %>%
  lapply(read.fasta, forceDNAtolower=FALSE, as.string=TRUE) %>%
  unlist(recursive=FALSE) %>%
  lapply(as.SeqFastadna)
gisaid_seq_names <- names(gisaid_seq) %>% word(., start=2, sep=fixed("|"))

write.fasta(gisaid_seq, names=gisaid_seq_names, 
            file.out=paste0("../data/gisaid_", curr_date, ".fasta"))

gisaid_metadata <- read_excel("../data/gisaid_cov2020_acknowledgement_table.xls", skip=2) %>%
  filter(rowMeans(is.na(.))<1) %>%
  rename(acc=`Accession ID`, virus_name=`Virus name`, country=Location,
         collection_date=`Collection date`) %>%
  mutate(collection_date_month_only=str_count(collection_date, "-")==1) %>%
  mutate(collection_date=parse_date(collection_date))


write.table(gisaid_metadata, paste0("../data/gisaid_", curr_date, "_metadata.tsv"), sep="\t", row.names = FALSE, quote = FALSE)








