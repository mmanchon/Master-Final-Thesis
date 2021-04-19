library(data.table)
library(tibble)
library(foreach)

SERIES_FILE <- "~/Proyecto TFM/Master-Final-Thesis/data/geo/GSE50911.txt/GSE50911_series_matrix.txt"

# Read characteristics
con <- file(SERIES_FILE, "r")
characteristics <- c()
while(TRUE) {
  line <- readLines(con, n=1)
  if(length(line) == 0) {
    break
  } else if(startsWith(line, "!Sample_title")) {
    titles <- unlist(strsplit(line, "\t"))[-1]
    titles <- gsub("\\\"", "", titles)
  } else if(startsWith(line, "!Sample_characteristics")) {
    characteristics <- c(characteristics, line)
  } else if(startsWith(line, "!Sample_geo_accession")) {
    accession <- unlist(strsplit(line, "\t"))[-1]
    accession <- gsub("\\\"", "", accession)
  }
}
close(con)



