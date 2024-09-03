library(tidyverse)
library(here)

read_fai <- \(file){ read_tsv(file, col_names = c("chr", "length", "offset", "linebases", "linewidth")) }

scaf_on_x <- read_tsv(here("data/mirang/mirang_sex_chrom.bed"))
genome <- read_fai(here("data/mirang/mirang_filt.fa.gz.fai")) |>
  arrange( -length ) |>
  filter( !(chr %in% scaf_on_x$chr) ) |>
  mutate(cum_length = cumsum(length),
         genome_share = cum_length / sum(length),
         genome_gain = genome_share - lag(genome_share, default = 0))

genome |> 
  filter(grepl("NC_", chr)) |> 
  pluck("chr") %>%
  str_c('"', . ,'"', collapse = ",") |> cat()
