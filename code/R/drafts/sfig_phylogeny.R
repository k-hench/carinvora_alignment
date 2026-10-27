library(tidyverse)
library(ggtree)
library(glue)
library(prismatic)
library(ape)
library(here)

tree <- read.tree(here("data/carnivora_short_labels.nwk"))
data <- read_tsv(here("data/carnivora_genome_and_timetree.tsv"))
#tree <- read.tree("results/neutral_tree/rerooted.tree")

treedata <- ggtree(tree, layout = "circular")@data |> 
  left_join(data |> rename(label = "spec")) |> 
  mutate(lab = glue("{str_to_title(str_replace(organism_name, '_', ' '))} [ {label} ]"))

clr2 <- "gray50"
mr <- 1.2
ggtree(treedata, layout = "circular", color = clr2) +
  geom_tiplab(aes(label = glue("[ {label} ]"), color = label == "arcgaz"), family = "Ubuntu mono")+
  geom_tiplab(aes(x = x + 17, label = glue("{str_to_title(str_replace(organism_name, '_', ' '))}")),
              color = clr2, fontface = "italic") +
  scale_color_manual(values = c(`TRUE` = "#6ca299", `FALSE` = "black"), guide = FALSE)+
  theme(plot.margin = margin(mr,mr,mr,mr, unit = "in"))


side <- 12
ggsave(here("results/img/ref_phylo.pdf"), width = side, height = side, device = cairo_pdf)
side#geom_tiplab(aes(label = lab), family = "Ubuntu mono")
