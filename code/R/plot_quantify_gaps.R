library(tidyverse)
library(ggstance)
fnt_sel <- "ubuntu mono"
data <- read_tsv("~/work/software/python/gapped_maf2fasta/data/base_counts.tsv",
                 skip = 1, col_types = "cidiiiiii") |>
  rename(sample = "# sample") |>
  mutate(sample = str_remove(sample, "# "))

data |>
  mutate(sample = fct_reorder(sample, -gaps)) |>
  pivot_longer(c(gaps,A:N)) |>
  mutate(name = factor(name, levels = c("gaps", "N", "A", "C", "G", "T"))) |>
  ggplot(aes(y = sample, x = value / (5940 * 1e3), fill = name)) +
  geom_barh(stat = "identity") +
  scale_fill_viridis_d(option = "C") +
  theme_minimal(base_family = fnt_sel)
