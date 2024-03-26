library(tidyverse)
library(prismatic)
library(ggtree)
library(geomtextpath)
library(ggnewscale)
library(patchwork)
fnt_sel <- "Arial"

clrs <- viridis::magma(5) |>
  clr_desaturate(.3) |>
  clr_darken(.2)
clr <- clrs[3]
read_bed <- \(file){
  read_tsv(file,
           col_names = c("seq", "start", "end", "cov")) |>
    mutate(length = end - start)
}


tree <- read.tree("data/carnivora_short_labels.nwk")

path <- "results/neutral_tree/cov/"

files <- dir(path, pattern = ".collapsed.bed.gz")

data <- str_c(path, files) |>
  map_dfr(read_bed)

data_summary <- data |>
  group_by(cov) |>
  summarise(length = sum(length)) |>
  ungroup() |>
  mutate(percent  = length / sum(length) * 100,
         cum_length = cumsum(length),
         cum_percent = cumsum(percent))

total_length <- sum(data_summary$length)
cov_steps <- c(0, 5, 11, 33, 41, 67)

p1 <- data_summary |>
  ggplot(aes(x = total_length-cum_length, y = cov + 1)) +
  geom_vline(xintercept = c(5e3, 50e3, 1e5)*1e3,
             linetype = 3,
             linewidth = .5) +
  # geom_hline(yintercept = 64) +
  geom_linerange(data = data_summary |>
                   filter((cov + 1) %in% cov_steps[2:5] ),
                 aes(xmin = 0, xmax = total_length-cum_length,
                     y = cov + 1),
                 linetype = 3, linewidth = .5) +
  geom_area(color = clr, fill = clr_alpha(clr)) +
  geom_linerange(data = data_summary |>
               filter((cov + 1) %in% cov_steps[2:5] ),
             aes(ymin = 0, ymax = cov + 1,
                 x = total_length-cum_length),
             linewidth = .25, color = clr_alpha("white")) +
  geom_text(data = tibble(y = (cov_steps[1:5] + cov_steps[2:6])/2,
                          n = diff(cov_steps)),
            aes(x = 2e8, label = n, y = y),
            color = "white",
            family = fnt_sel) +
  scale_x_continuous(labels = \(x){sprintf("%.1f",x*1e-9)}) +
  labs(x ="sequence length (Gb)", y = "coverage (n genomes)") +
  coord_cartesian(expand = 0) +
  theme_minimal(base_family = fnt_sel) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(linewidth = .5))

data_expansions <- tibble(y1 = 67 - c(22, 16, 0, 0, 0) +.5 ,
                          y2 = 67 - c(27, 27, 33, 41, 67)+.5,
                          x = 56 + (5:1)*3,
                          labs = str_c(rep(c("","+ "), c(1,4)), diff(cov_steps)),
                          group = seq_along(x))

p2 <- tree |>
  ggtree( layout='circular', color = "transparent") +
  geom_hilight(node = 120, fill = clrs[5]) +
  geom_hilight(node = 124, fill = clrs[4]) +
  geom_hilight(node = 104, fill = clrs[3]) +
  geom_hilight(node = 129, fill = clrs[3]) +
  geom_hilight(node = 95, fill = clrs[2]) +
  geom_hilight(node = 69, fill = clrs[1]) +
  geom_tree()+
  geom_tiplab(aes(x = x + 18.5, color = label == "arcgaz"), family = "ubuntu mono") +
  scale_color_manual(values = c(`TRUE` = clr, `FALSE` = "gray60"), guide = "none") +
  new_scale_color() +
  geom_textsegment(data = data_expansions,
                aes(x = x, xend = x, y = y1, yend = y2, label = labs,
                    color = factor(group)),
                family = fnt_sel) +
  scale_color_manual(values = rev(clrs), guide = "none") +
  scale_fill_manual(values = clrs |> set_names(nm = str_c("inner", seq_along(clrs))), guide = "none") +
  scale_y_continuous(limits = c(.25, 67.75))

pp <- p1 + p2 + plot_layout(widths = c(.75, 1))

ggsave(plot = pp,
      filename = "results/img/hal_coverage.pdf",
      width = 12, height = 6,
      device = cairo_pdf)
