# Run from the project directory: source("run.R") or Rscript run.R.
# Install dependencies once; this script does not install or update packages:
# install.packages(c("ggplot2", "ggrepel", "patchwork", "ragg", "svglite", "circlize"),
#                  repos = "https://cloud.r-project.org")
# install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap", ask = FALSE, update = FALSE)
# Tested: R 4.6.1; ggplot2 4.0.3; ggrepel 0.9.8; patchwork 1.3.2;
#         ragg 1.5.2; svglite 2.2.2; ComplexHeatmap 2.28.0; circlize 0.4.18.
required <- c('ggplot2', 'ggrepel', 'patchwork', 'ragg', 'svglite', 'ComplexHeatmap', 'circlize')
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop('Missing packages: ', paste(missing, collapse = ', '),
                          '. See installation commands at the top of run.R.', call. = FALSE)
source('plots.R')

# ----- Configuration -----
regenerate_example <- FALSE # TRUE overwrites the three example CSVs. Keep FALSE for real data.
lfc_cutoff <- 1
padj_cutoff <- .05
heatmap_gene_count <- 24
representative_gene_count <- 4
image_dpi <- 180
if (regenerate_example) make_example_data()

# ----- Load inputs -----
expression <- as.matrix(read.csv('data/expression.csv', row.names = 1, check.names = FALSE))
samples <- read.csv('data/samples.csv')
results <- read.csv('data/differential_expression.csv')
invisible(check_expression(expression, samples))
if (!setequal(rownames(expression), results$gene)) stop('Expression and result gene IDs must match.')
selected_genes <- select_heatmap_genes(results, heatmap_gene_count, padj_cutoff, lfc_cutoff)
# Choose two up and two down genes where possible, with deterministic fallback.
ranked <- results[match(select_heatmap_genes(results, nrow(results), padj_cutoff, lfc_cutoff), results$gene), ]
representatives <- unique(c(head(ranked$gene[ranked$log2FC > 0], ceiling(representative_gene_count/2)),
                            head(ranked$gene[ranked$log2FC < 0], floor(representative_gene_count/2)), ranked$gene))
representatives <- head(representatives, representative_gene_count)

plots <- list(
  pca = plot_pca(expression, samples),
  correlation = plot_sample_correlation(expression, samples),
  volcano = plot_volcano(results, lfc_cutoff, padj_cutoff, label_n = 4, highlight = representatives),
  heatmap = plot_expression_heatmap(expression, samples, genes = selected_genes),
  gene_expression = plot_gene_expression(expression, samples, representatives)
)
plots$combined <- compose_figure(plots$pca, plots$correlation, plots$volcano,
                                plots$heatmap, plots$gene_expression, nrow(expression), ncol(expression))

# ----- Export five individual plots and one composite figure -----
dir.create('figures', showWarnings = FALSE)
for (name in names(plots)) {
  size <- switch(name, combined = c(15, 15), gene_expression = c(12, 4.5), c(8.5, 6.5))
  ggplot2::ggsave(file.path('figures', paste0(name, '.png')), plots[[name]],
    width = size[1], height = size[2], dpi = image_dpi, bg = 'white', device = ragg::agg_png)
  ggplot2::ggsave(file.path('figures', paste0(name, '.svg')), plots[[name]],
    width = size[1], height = size[2], bg = 'white', device = svglite::svglite)
}
expected <- file.path('figures', as.vector(outer(names(plots), c('.png', '.svg'), paste0)))
if (any(!file.exists(expected)) || any(file.info(expected)$size <= 0)) stop('Figure export incomplete.')
message('Done: 5 individual plots + 1 combined figure, each in PNG and SVG.')
