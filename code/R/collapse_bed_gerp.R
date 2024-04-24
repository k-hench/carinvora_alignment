library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
bed_in <- args[[1]]
bed_out <- args[[2]]

bed <- read_tsv(bed_in, col_names = c("seq", "start", "end", "gerp"))

bed |>
  group_by(seq) |>
  mutate(gerp_bin = cumsum(!(gerp == lag(gerp, default = NA) & start == lag(end, default = 0) ))) |>
  group_by(seq, cov_bin) |>
  summarize(start = min(start),
            end = max(end),
            gerp = gerp[[1]]) |>
  ungroup() |>
  select(-gerp_bin) |>
  write_tsv( bed_out, col_names = FALSE )