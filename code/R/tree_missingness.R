library(tidyverse)
library(ggtree)


tree <- read.tree(file = "data/carnivora_short_labels.nwk")


data <- read_tsv("results/neutral_tree/multifa/basereport.tsv", skip = 1) |>
  rename(sample = "# sample") |>
  mutate(sample = str_remove(sample, "# ")) |>
  left_join(ggtree(tree)$data |>
              filter(isTip) |>
              select(sample = label, x, y)) |>
  mutate( a = A/n,
          c = C/n,
          g = G/n,
          t = `T`/n)

wdth <- .9/2
xmn <- 65
xwdth <- 12.5
scl <- \(x){scales::rescale(x, from = c(0,1), to = c(0,xwdth))}
xbse2 <- xmn
ggtree(tree) +
  geom_tiplab() +
  geom_vline(xintercept = seq(from = xmn, to = xmn + xwdth, length.out = 5),
             linewidth = .2, color = rgb(0,0,0,.4)) +
  geom_rect(data = data,
            aes(ymin = y - wdth, ymax = y + wdth,
                xmin = xbse2, xmax = xbse2 + scl(a),
                fill = "A")) +
  geom_rect(data = data,
            aes(ymin = y - wdth, ymax = y + wdth,
                xmin = xbse2 + scl(a),
                xmax = xbse2  + scl(a) + scl(c),
                fill = "C")) +
  geom_rect(data = data,
            aes(ymin = y - wdth, ymax = y + wdth,
                xmin = xbse2 + scl(a) + scl(c),
                xmax = xbse2 + scl(a) + scl(c) + scl(g),
                fill = "G")) +
  geom_rect(data = data,
            aes(ymin = y - wdth, ymax = y + wdth,
                xmin = xbse2 + scl(a) + scl(c)  + scl(g),
                xmax = xbse2 + scl(a) + scl(c)  + scl(g) + scl(t),
                fill = "T")) +
  geom_rect(data = data,
            aes(ymin = y - wdth, ymax = y + wdth,
                xmin = xbse2 + scl(a) + scl(c)  + scl(g),
                xmax = xbse2 + scl(a) + scl(c)  + scl(g) + scl(t),
                fill = "T")) +
  geom_rect(data = data,
            aes(ymin = y - wdth, ymax = y + wdth,
                xmin = xbse2 + scl(a) + scl(c)  + scl(g) + scl(t),
                xmax = xbse2 + scl(a) + scl(c)  + scl(g) + scl(t) + scl(gaps/n),
                fill = "gap")) +
  scale_fill_manual(values = prismatic::clr_alpha(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"),.75)) +
  coord_cartesian(xlim = c(0,  75))

