#!/usr/bin/env Rscript
# Rscript --vanilla gerp_fixation_levels.R <anc.tree> <clade_tree.pdf> <out.vcf> <in_tree>
args <- commandArgs(trailingOnly = TRUE)
print(args)

library(tidyverse)
library(glue)
library(ggtree)
# head -n 1 carnivora_set.txt

tree_in <- args[1]
pdf_out <- args[2]
vcf_out <- args[3]

fnt_sel = "Arial"

tree_anc <- read.tree(tree_in)
all_tips <- sort(tree_anc$tip.label)

tree_anc_data <- (tree_anc |>
  ggtree::groupClade(.node = c(94,95,102, 103, 104, 119,129,
                               69,88,70, 84, 71)) |>
  ggtree(aes(color = group)))$data

group_labels <- tibble(
  label = str_c("Anc", c("00", "01", "02", "05",
                         "06", "13", "12","22", "21",
                         "04", "03","08","07")),
  gr_label = c("Carnivora", "Feliformia","Caniformia","Canidae",
               "Arctoidea", "Ursidae", "*unnamed*","Pinnipedia", "Musteloidea",
               "Viverroidea","Felidae","Pantherinae","Felinae"),
  gr_idx = c(0,8,1,2,
             3,7,4,6,5,
             9,10,11,12)
)

p <- tree_anc_data |>
  left_join(group_labels) |>
  ggtree(aes(color = group)) +
  geom_nodelab(aes(label = str_c(str_remove(label, "Anc"), ":", node), filter = !isTip),
                 color = "red", family = fnt_sel) +
  geom_nodelab(aes(label =gr_label, filter = !isTip,
                 color = factor(gr_idx)), vjust = 1.8, family = fnt_sel) +
  geom_tiplab(family = fnt_sel) +
  scale_color_manual(values = c("black",
                                RColorBrewer::brewer.pal(4,"Purples"),
                                RColorBrewer::brewer.pal(3,"Reds"),
                                RColorBrewer::brewer.pal(3,"Greens"),
                                RColorBrewer::brewer.pal(3,"Blues")[2:3]))

ggsave(plot = p,
       filename = pdf_out,
       height = 9,
       width = 8,
       device = cairo_pdf)

# group_labels |> arrange(gr_idx)

# helper function to collect all species belonging to a clade
extract_species_for_group <- \(gr_contains){
  tree_anc_data |>
    filter(isTip,
           group %in% gr_contains) |>
    pluck("label")
}

# helper function to create genotypes based on assignment to a clade
create_fixed_genotypes <- \(spec_in){
  tibble(spec = all_tips,
         GT = if_else(spec %in% spec_in, "1/1", "0/0")) |>
    pivot_wider(names_from = spec, values_from = GT)
}

# compile the genotypes fixed within clades
group_structure <- tibble(
  gr_idx = 0:12,
  contains = list(1:12, c(2,3,4,5,6,7),
                  2, c(4,5,6,7),
                  c(5,6), 5, 6, 7,
                  c(9,10,11,12),
                  9, c(11,12), 11, 12) ) |>
  left_join(group_labels) |>
  mutate(species = map(contains, extract_species_for_group),
         genotypes = map(species, create_fixed_genotypes)
         )

# group_structure$species[[which(group_structure$gr_idx == 9)]]

# piecing together the elements of the vcf file
vcf_head <- glue('##fileformat=VCFv4.2
##fileDate={Sys.Date()}
##contig=<ID=dummy,length=12,species="none">
##INFO=<ID=FI,Number=1,Type=String,Description="Clade in which the SNP is fixed.">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
vcf_cols <- str_c(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", all_tips),
                  collapse = "\t")

# exporting the vcf
vcf_head |> write_lines(file = vcf_out)
vcf_cols |> write_lines(file = vcf_out, append = TRUE)

group_structure |>
  mutate(chr = "dummy",
         pos = gr_idx,
         id = ".",
         ref = "A",
         alt = "T",
         qual = ".",
         filter = "PASS",
         info = str_c("FI=",gr_label),
         format = "GT") |>
  filter(gr_idx > 0) |>
  select(chr:format,genotypes) |>
  unnest(genotypes) |>
  write_tsv(file = vcf_out,
            col_names = FALSE,
            append = TRUE)
