library(tidyverse)
library(patchwork)
library(ggtree)
library(phytools)

tree <- ape::read.tree(file = "results/neutral_tree/multifa/combined_windows.fa.treefile")
orgtree <- ape::read.tree(file = "data/topology.tree")

p0 <- ggtree(orgtree)+
  geom_tiplab(aes(label = label))+
  geom_nodelab(aes(label = node))+
  coord_cartesian(xlim = c(0, 15))

tree_rerooted <- reroot(tree, node.number = as.numeric(args[[2]]))

root_node <- 75

purrr::reduce(c(68,130,133,121,108,101,78,81,92,95),
                    ape::rotate,
                    .init = reroot(tree, node.number = as.numeric(root_node)))  |>
  ggtree()+
  geom_tiplab(aes(label = label))+
  geom_nodelab(aes(label = node))+
  coord_cartesian(xlim = c(0, 0.25)) + p0
