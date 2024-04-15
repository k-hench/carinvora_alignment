library(phytools)
options(scipen = 6)
args <- commandArgs(trailingOnly = TRUE)

tree <- ape::read.tree(file = args[[1]])
tree_rerooted <- reroot(tree, node.number = as.numeric(args[[2]]))

ape::write.tree(phy = purrr::reduce(c(68,130,133,121,108,101,78,81,92,95),
                                    ape::rotate,
                                    .init = tree_rerooted),
                file = args[[3]], digits = 6)
