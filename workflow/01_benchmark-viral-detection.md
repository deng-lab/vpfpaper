01_benchmark-viral-detection
================
Jinlong Ru

- <a href="#tasks" id="toc-tasks">Tasks</a>
- <a href="#files-written" id="toc-files-written">Files written</a>

Updated: 2023-01-28 23:29:34 CET.

## Tasks

``` r
library(conflicted)
library(here)
library(vpfpaper)
```

    Loading required package: tidyverse

    ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
    ✔ ggplot2 3.3.6      ✔ purrr   0.3.5 
    ✔ tibble  3.1.8      ✔ dplyr   1.0.10
    ✔ tidyr   1.2.1      ✔ stringr 1.4.1 
    ✔ readr   2.1.3      ✔ forcats 0.5.2 

``` r
data("vpf_prab", package = "vpfpaper")
data("vpf_relab", package = "vpfpaper")
```

``` r
# plot with threshold as x-axis, value as y-axis, facet by metric, and fill by method
p <- vpf_prab %>%
  mutate(value = replace_na(value, 0)) %>%
  ggplot(aes(x = threshold, y = value, color = method)) +
  geom_point(aes(group = method), size = 0.05, alpha = 0.05) +
  geom_smooth(aes(color = method), method = "loess") +
  scale_color_manual(values = c('#4daf4a','#377eb8','#e41a1c')) +
  facet_grid(metric ~ taxalevel) +
  labs(x = "Abundance threshold (%)", y = "") +
  # ggtitle("Precision, recall and F1 score of different methods") +
  theme_bw() +
  theme(axis.text = element_text(size = 6),
        legend.title = element_blank())

ggsave(p, filename = path_target("precision_recall_f1_score_contig.pdf"), width = 10, height = 3.8)
```

    `geom_smooth()` using formula 'y ~ x'

``` r
p2 <- ggplot(vpf_relab, aes(x = taxalevel, y = bcs, fill = method)) +
  geom_boxplot() +
  scale_fill_manual(values = alpha(c('#4daf4a','#377eb8','#e41a1c'), 0.6)) +
  labs(x = "Taxonomy level", y = "Bray-Curtis dissimilarity") +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        legend.title = element_blank())

ggsave(p2, filename = path_target("BC4abundance_contig.pdf"), width = 10, height = 2.3)
```

``` r
library(patchwork)

p / p2 +
  plot_layout(heights = c(5, 3)) &
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12))
```

    `geom_smooth()` using formula 'y ~ x'

![](01_benchmark-viral-detection_files/figure-gfm/Combine%20plot-1.png)

``` r
ggsave(path_target("fig2_benchmark_contig.png"), width = 10, height = 6, dpi = 300)
```

    `geom_smooth()` using formula 'y ~ x'

## Files written

These files have been written to the target directory,
`data/01_benchmark-viral-detection`:

``` r
projthis::proj_dir_info(path_target())
```

    # A tibble: 3 × 4
      path                                 type         size modification_time  
      <fs::path>                           <fct> <fs::bytes> <dttm>             
    1 BC4abundance_contig.pdf              file        6.78K 2023-01-28 22:29:41
    2 fig2_benchmark_contig.png            file      950.72K 2023-01-28 22:29:49
    3 precision_recall_f1_score_contig.pdf file        2.06M 2023-01-28 22:29:41
