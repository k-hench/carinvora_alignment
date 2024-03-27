library(tidyverse)
library(prismatic)

read_and_summarize <- \(file, lab){
  read_tsv(file, col_names = c("seq", "start", "end")) |>
    mutate(length = end - start,
           length_class = cut(length, breaks = c(0, 1e3, 5e3, 1e4, 5e4, 1e5, 5e5, 1e6, Inf))) |>
    group_by(length_class) |>
    summarize(n = n(),
              sum = sum(length)) |>
    ungroup() |>
    mutate(class = lab)
}

data <- tibble(file = str_c("~/Downloads/available_space", c("_no_cds.bed", ".bed")),
       lab = c("cds", "avail")) |>
  pmap_dfr(read_and_summarize)

data |>
  pivot_longer(n:sum) |>
  mutate(length_class = fct_reorder(length_class, -as.numeric(length_class))) |>
  ggplot(aes(x = class, y = value, color = length_class)) +
  geom_bar(stat = "identity", aes(fill = after_scale(prismatic::clr_alpha(color,.7))))+
  facet_wrap(name~., scales = "free")+
  scale_color_manual(values = viridis::plasma(8))

