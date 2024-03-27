library(tidyverse)

data <- read_tsv("~/Downloads/available_space.bed", col_names = c("seq", "start", "end")) |>
  mutate(length = end - start,
         length_class = cut(length, breaks = c(0, 1e3, 5e3, 1e4, 5e4, 1e5, 5e5, 1e6, Inf)))

data_summary <- data |>
  group_by(length_class) |>
  summarize(n = n(),
            sum = sum(length)) |>
  ungroup()

data_summary |>
  pivot_longer(n:sum) |>
  ggplot(aes(x = 1, y = log10(value), fill = length_class)) +
  geom_bar(stat = "identity")+
  facet_wrap(name~., scales = "free")
