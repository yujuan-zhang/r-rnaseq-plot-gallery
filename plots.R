# Reusable plotting functions. Sourcing this file does not install packages or generate data.
# CRAN: install.packages(c("ggplot2", "ggrepel", "patchwork", "ragg", "svglite", "circlize"))
# Bioconductor: install.packages("BiocManager")
#               BiocManager::install("ComplexHeatmap", ask = FALSE, update = FALSE)
# Tested versions are listed in run.R. Use normalized, appropriately log2-transformed expression.

GROUP_COLORS <- c(Control = '#287A8A', Treatment = '#C65A46')
DE_COLORS <- c(Down = '#287A8A', 'Not significant' = '#CCD3DA', Up = '#C65A46')

check_columns <- function(data, required) {
  if (!is.data.frame(data) || !all(required %in% names(data)))
    stop('Required columns: ', paste(required, collapse = ', '))
  if (!nrow(data)) stop('Input must contain at least one row.')
}
check_probability <- function(x, name = 'padj') {
  if (!is.numeric(x) || any(!is.finite(x) | x < 0 | x > 1))
    stop(name, ' must contain finite numbers in [0, 1]; handle NA explicitly.')
}
check_expression <- function(expression, samples) {
  check_columns(samples, c('sample', 'group'))
  if (!is.matrix(expression) || !is.numeric(expression) ||
      nrow(expression) < 2 || ncol(expression) < 3 || any(!is.finite(expression)))
    stop('expression must be a finite numeric matrix: genes x samples (at least 2 x 3).')
  valid_ids <- function(x) !is.null(x) && !anyNA(x) && !anyDuplicated(x) && all(nzchar(as.character(x)))
  if (!valid_ids(rownames(expression)) || !valid_ids(colnames(expression)))
    stop('Unique, nonempty gene and sample names are required.')
  if (!valid_ids(samples$sample) || anyNA(samples$group) ||
      any(!nzchar(as.character(samples$group))) || !setequal(samples$sample, colnames(expression)))
    stop('Sample metadata must match expression columns exactly; groups cannot be missing.')
  samples <- samples[match(colnames(expression), samples$sample), , drop = FALSE]
  samples$group <- factor(samples$group, levels = unique(as.character(samples$group)))
  samples
}
group_palette <- function(groups) {
  groups <- unique(as.character(groups))
  if (all(groups %in% names(GROUP_COLORS))) return(GROUP_COLORS[groups])
  stats::setNames(grDevices::hcl.colors(length(groups), 'Dark 3'), groups)
}
plot_theme <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(color = '#25354A'),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = '#EDF0F3', linewidth = 0.3),
      plot.title = ggplot2::element_text(face = 'bold', size = base_size + 2),
      plot.subtitle = ggplot2::element_text(color = '#66768A', size = base_size - 1),
      plot.title.position = 'plot',
      legend.title = ggplot2::element_text(size = base_size - 1),
      legend.text = ggplot2::element_text(size = base_size - 1),
      plot.margin = ggplot2::margin(10, 12, 10, 10))
}

# All panels share these synthetic data. Regeneration is explicit, not automatic.
make_example_data <- function(seed = 20260930, directory = 'data') {
  set.seed(seed)
  n <- 24L; g <- 800L
  samples <- data.frame(sample = sprintf('S%02d', seq_len(n)),
                         group = rep(c('Control', 'Treatment'), each = n / 2),
                         batch = rep(rep(c('B1', 'B2'), each = 6), 2))
  base <- runif(g, 5, 10)
  gene_sd <- runif(g, 0.35, 0.85)
  effects <- c(runif(100, 1.2, 2.1), -runif(80, 1.2, 2.1), rnorm(g - 180, 0, .08))
  latent <- rnorm(n, 0, .25)
  expression <- matrix(base, g, n) + matrix(rnorm(g*n), g, n) * gene_sd +
    outer(effects, as.numeric(samples$group == 'Treatment')) +
    outer(rnorm(g, 0, .6), latent) +
    outer(rnorm(g, 0, .12), as.numeric(samples$batch == 'B2'))
  dimnames(expression) <- list(sprintf('Gene%04d', seq_len(g)), samples$sample)
  control <- samples$group == 'Control'; treatment <- !control
  # Demonstration only: Welch t tests on synthetic log2 expression, BH across all genes.
  # This is NOT a replacement for an RNA-seq count model or a batch-adjusted design.
  p <- apply(expression, 1, function(y) stats::t.test(y[treatment], y[control])$p.value)
  results <- data.frame(gene = rownames(expression),
                        log2FC = rowMeans(expression[, treatment]) - rowMeans(expression[, control]),
                        padj = stats::p.adjust(p, method = 'BH'))
  dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(data.frame(gene = rownames(expression), expression, check.names = FALSE),
                   file.path(directory, 'expression.csv'), row.names = FALSE)
  utils::write.csv(samples, file.path(directory, 'samples.csv'), row.names = FALSE)
  utils::write.csv(results, file.path(directory, 'differential_expression.csv'), row.names = FALSE)
  invisible(list(expression = expression, samples = samples, results = results))
}

# Center genes without unit-variance scaling; PCA uses all nonconstant input genes.
# Ellipses describe the score distribution, not confidence intervals for group means.
plot_pca <- function(expression, samples, ellipse_level = 0.8, labels = TRUE) {
  samples <- check_expression(expression, samples)
  if (length(ellipse_level) != 1 || !is.finite(ellipse_level) || ellipse_level <= 0 || ellipse_level >= 1)
    stop('ellipse_level must be between 0 and 1.')
  keep <- apply(expression, 1, stats::sd) > 0
  if (sum(keep) < 2) stop('PCA requires at least two variable genes.')
  fit <- stats::prcomp(t(expression[keep, , drop = FALSE]), center = TRUE, scale. = FALSE)
  scores <- cbind(samples, as.data.frame(fit$x[, 1:2, drop = FALSE]))
  if (fit$sdev[2] < sqrt(.Machine$double.eps) * fit$sdev[1]) stop('PCA requires at least two nonzero components.')
  variance <- 100 * fit$sdev^2 / sum(fit$sdev^2)
  # Skip undersized or singular groups rather than claiming a meaningful ellipse.
  usable <- vapply(split(scores, scores$group), function(d)
    nrow(d) >= 4 && det(stats::cov(d[, c('PC1', 'PC2')])) > 1e-10, logical(1))
  ellipse_data <- scores[as.character(scores$group) %in% names(usable)[usable], , drop = FALSE]
  p <- ggplot2::ggplot(scores, ggplot2::aes(PC1, PC2, color = group)) +
    ggplot2::geom_hline(yintercept = 0, color = '#DDE3E9', linewidth = .3) +
    ggplot2::geom_vline(xintercept = 0, color = '#DDE3E9', linewidth = .3)
  if (nrow(ellipse_data)) p <- p + ggplot2::stat_ellipse(data = ellipse_data,
    ggplot2::aes(fill = group), geom = 'polygon', type = 'norm', level = ellipse_level,
    alpha = .09, linewidth = .45, show.legend = FALSE)
  if ('batch' %in% names(scores)) {
    if (anyNA(scores$batch) || any(!nzchar(as.character(scores$batch))) || length(unique(scores$batch)) > 6)
      stop('batch must be nonempty and have at most six levels for point shapes.')
    p <- p + ggplot2::geom_point(ggplot2::aes(shape = batch), size = 3)
  } else p <- p + ggplot2::geom_point(size = 3)
  if (labels) p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = sample),
    size = 2.5, seed = 42, max.overlaps = Inf, show.legend = FALSE)
  p + ggplot2::scale_color_manual(values = group_palette(scores$group)) +
    ggplot2::scale_fill_manual(values = group_palette(scores$group)) +
    ggplot2::labs(title = 'Sample structure',
      subtitle = sprintf('All variable genes | %s%% distribution ellipses where estimable', 100 * ellipse_level),
      x = sprintf('PC1 (%.1f%%)', variance[1]), y = sprintf('PC2 (%.1f%%)', variance[2]),
      color = 'Group', shape = 'Batch') + plot_theme() + ggplot2::theme(legend.position = 'bottom')
}

# Return a grid grob suitable for print via grid.draw() or patchwork::wrap_elements().
# Capture at a known physical size; no anno_zoom/anno_link or device-dependent annotations.
heatmap_grob <- function(ht, title = NULL, width = 7, height = 5.5) {
  grid::grid.grabExpr(ComplexHeatmap::draw(ht, heatmap_legend_side = 'right',
    annotation_legend_side = 'right', merge_legends = TRUE,
    column_title = title, column_title_gp = grid::gpar(fontsize = 13, fontface = 'bold', col = '#25354A')), width = width, height = height)
}
sample_annotations <- function(samples) {
  df <- data.frame(Group = samples$group)
  colors <- list(Group = group_palette(samples$group))
  if ('batch' %in% names(samples)) {
    if (anyNA(samples$batch) || any(!nzchar(as.character(samples$batch)))) stop('batch must be nonempty.')
    df$Batch <- samples$batch
    colors$Batch <- stats::setNames(grDevices::colorRampPalette(c('#C3CAD3', '#586980'))(length(unique(samples$batch))), unique(samples$batch))
  }
  ComplexHeatmap::HeatmapAnnotation(df = df, col = colors,
    simple_anno_size = grid::unit(2.5, 'mm'), annotation_name_gp = grid::gpar(fontsize = 8),
    annotation_legend_param = list(title_gp = grid::gpar(fontsize = 9), labels_gp = grid::gpar(fontsize = 8)))
}
plot_sample_correlation <- function(expression, samples) {
  # Some heatmap annotations measure text during construction. Avoid Rplots.pdf.
  if (grDevices::dev.cur() == 1L) {
    grDevices::pdf(file = NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  samples <- check_expression(expression, samples)
  if (any(apply(expression, 2, stats::sd) == 0)) stop('Constant sample profiles cannot be correlated.')
  correlation <- stats::cor(expression, method = 'pearson')
  hc <- stats::hclust(stats::as.dist(pmax(1 - correlation, 0)), method = 'average')
  lower <- min(.95, floor(min(correlation) * 20) / 20)
  breaks <- c(lower, (lower + 1) / 2, 1)
  ht <- ComplexHeatmap::Heatmap(correlation, name = 'Pearson r',
    col = circlize::colorRamp2(breaks, c('#F0F3F6', '#83B1B9', '#245368')),
    cluster_rows = hc, cluster_columns = hc, row_dend_reorder = FALSE, column_dend_reorder = FALSE,
    top_annotation = sample_annotations(samples),
    show_row_names = TRUE, show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 7), column_names_gp = grid::gpar(fontsize = 7),
    column_names_rot = 45, rect_gp = grid::gpar(col = 'white', lwd = .3),
    heatmap_legend_param = list(at = breaks, title_gp = grid::gpar(fontsize = 9), labels_gp = grid::gpar(fontsize = 8)),
    column_title = NULL,
    border = FALSE)
  heatmap_grob(ht, 'Sample correlation')
}

# Filter and rank an existing differential-expression table; no tests are run here.
select_heatmap_genes <- function(results, n = 24, padj_cutoff = .05, lfc_cutoff = 1) {
  check_columns(results, c('gene', 'log2FC', 'padj'))
  check_probability(results$padj)
  if (!is.numeric(results$log2FC) || any(!is.finite(results$log2FC)) ||
      anyNA(results$gene) || anyDuplicated(results$gene) || any(!nzchar(as.character(results$gene)))) stop('Invalid differential-expression table.')
  if (length(n) != 1 || !is.finite(n) || n < 2 || n != floor(n)) stop('n must be an integer >= 2.')
  if (length(padj_cutoff) != 1 || !is.finite(padj_cutoff) || padj_cutoff <= 0 || padj_cutoff >= 1 ||
      length(lfc_cutoff) != 1 || !is.finite(lfc_cutoff) || lfc_cutoff <= 0) stop('Invalid selection thresholds.')
  eligible <- results[results$padj <= padj_cutoff & abs(results$log2FC) >= lfc_cutoff, ]
  eligible <- eligible[order(eligible$padj, -abs(eligible$log2FC), eligible$gene), ]
  if (nrow(eligible) < 2) stop('Fewer than two genes meet the heatmap thresholds.')
  head(eligible$gene, n)
}
plot_expression_heatmap <- function(expression, samples, genes = rownames(expression)) {
  # Some heatmap annotations measure text during construction. Avoid Rplots.pdf.
  if (grDevices::dev.cur() == 1L) {
    grDevices::pdf(file = NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  samples <- check_expression(expression, samples)
  if (length(genes) < 2 || anyNA(genes) || anyDuplicated(genes) || !all(genes %in% rownames(expression)))
    stop('genes must contain at least two unique gene IDs present in the matrix.')
  x <- expression[genes, , drop = FALSE]
  if (any(apply(x, 1, stats::sd) == 0)) stop('Remove constant genes before row scaling.')
  z <- t(scale(t(x)))
  if (length(genes) < 3) split <- NULL else split <- factor(cutree(hclust(dist(z)), k = 2), labels = c('Module 1', 'Module 2'))
  side <- ComplexHeatmap::rowAnnotation(Mean = ComplexHeatmap::anno_barplot(rowMeans(x),
    gp = grid::gpar(fill = '#8797AB', col = NA), border = FALSE,
    axis_param = list(gp = grid::gpar(fontsize = 6)), width = grid::unit(11, 'mm')),
    annotation_name_gp = grid::gpar(fontsize = 8))
  ht <- ComplexHeatmap::Heatmap(z, name = 'Row Z',
    col = circlize::colorRamp2(c(-2, 0, 2), c('#287A8A', '#FAFAF8', '#C65A46')),
    cluster_rows = TRUE, cluster_columns = TRUE, row_split = split,
    column_split = samples$group, cluster_column_slices = FALSE, cluster_row_slices = FALSE,
    top_annotation = sample_annotations(samples), right_annotation = side,
    row_gap = grid::unit(2, 'mm'), column_gap = grid::unit(3, 'mm'),
    row_title_gp = grid::gpar(fontsize = 8), column_title_gp = grid::gpar(fontsize = 10, fontface = 'bold'),
    row_names_gp = grid::gpar(fontsize = 7), column_names_gp = grid::gpar(fontsize = 7),
    column_names_rot = 45, rect_gp = grid::gpar(col = 'white', lwd = .25),
    heatmap_legend_param = list(at = c(-2, 0, 2), title_gp = grid::gpar(fontsize = 9), labels_gp = grid::gpar(fontsize = 8)))
  heatmap_grob(ht, 'Selected expression modules')
}

plot_volcano <- function(results, lfc_cutoff = 1, padj_cutoff = .05, label_n = 8, highlight = NULL) {
  check_columns(results, c('gene', 'log2FC', 'padj'))
  check_probability(results$padj)
  if (!is.numeric(results$log2FC) || any(!is.finite(results$log2FC))) stop('log2FC must be finite.')
  if (anyNA(results$gene) || anyDuplicated(results$gene) || any(!nzchar(as.character(results$gene)))) stop('gene IDs must be unique and nonempty.')
  if (length(lfc_cutoff) != 1 || !is.finite(lfc_cutoff) || lfc_cutoff <= 0 ||
      length(padj_cutoff) != 1 || !is.finite(padj_cutoff) || padj_cutoff <= 0 || padj_cutoff >= 1 ||
      length(label_n) != 1 || !is.finite(label_n) || label_n < 0 || label_n != floor(label_n)) stop('Invalid plot parameters.')
  if (anyNA(highlight) || !all(highlight %in% results$gene)) stop('Unknown highlight gene.')
  d <- results
  d$status <- factor(ifelse(d$padj <= padj_cutoff & d$log2FC >= lfc_cutoff, 'Up',
    ifelse(d$padj <= padj_cutoff & d$log2FC <= -lfc_cutoff, 'Down', 'Not significant')),
    levels = names(DE_COLORS))
  d$neglog10 <- -log10(pmax(d$padj, .Machine$double.xmin))
  selected <- d[d$status != 'Not significant', , drop = FALSE]
  selected <- head(selected[order(selected$padj, -abs(selected$log2FC), selected$gene), ], label_n)
  selected <- d[d$gene %in% unique(c(selected$gene, highlight)), , drop = FALSE]
  n <- table(d$status)
  p <- ggplot2::ggplot(d, ggplot2::aes(log2FC, neglog10)) +
    ggplot2::geom_hline(yintercept = -log10(padj_cutoff), linetype = 'dashed', color = '#AAB4C1', linewidth = .35) +
    ggplot2::geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = 'dashed', color = '#AAB4C1', linewidth = .35) +
    ggplot2::geom_point(ggplot2::aes(color = status), size = 1.5, alpha = .8) +
    ggplot2::geom_point(data = d[d$gene %in% highlight, , drop = FALSE], shape = 21, size = 3, stroke = .7, fill = NA, color = '#25354A') +
    ggrepel::geom_text_repel(data = selected, ggplot2::aes(label = gene), size = 2.6,
      seed = 42, box.padding = .35, min.segment.length = 0, max.overlaps = Inf) +
    ggplot2::scale_color_manual(values = DE_COLORS, drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(.03, .18))) +
    ggplot2::labs(title = 'Differential expression',
      subtitle = sprintf('Down: %d   |   Not significant: %d   |   Up: %d', n[1], n[2], n[3]),
      x = 'log2 fold change (Treatment - Control)', y = expression(-log[10]('adjusted P')),
      color = NULL) + plot_theme() + ggplot2::theme(legend.position = 'bottom')
  p
}

plot_gene_expression <- function(expression, samples, genes) {
  samples <- check_expression(expression, samples)
  if (!length(genes) || anyNA(genes) || anyDuplicated(genes) || !all(genes %in% rownames(expression))) stop('Unknown or duplicated gene IDs.')
  if (any(table(samples$group) < 3)) stop('Violin plots require at least three observations per group.')
  x <- expression[genes, , drop = FALSE]
  d <- data.frame(gene = factor(rep(genes, ncol(x)), levels = genes),
    sample = rep(colnames(x), each = nrow(x)), expression = as.vector(x))
  d$group <- samples$group[match(d$sample, samples$sample)]
  ggplot2::ggplot(d, ggplot2::aes(group, expression, fill = group)) +
    ggplot2::geom_violin(alpha = .2, linewidth = .35, trim = FALSE, color = NA) +
    ggplot2::geom_boxplot(width = .22, outlier.shape = NA, fill = 'white', color = '#536477', linewidth = .4) +
    ggplot2::geom_point(ggplot2::aes(color = group),
      position = ggplot2::position_jitter(width = .09, seed = 42), size = 1.6) +
    ggplot2::facet_wrap(ggplot2::vars(gene), nrow = 1, scales = 'free_y') +
    ggplot2::scale_fill_manual(values = group_palette(samples$group)) +
    ggplot2::scale_color_manual(values = group_palette(samples$group)) +
    ggplot2::labs(title = 'Representative genes', subtitle = 'Each dot is one sample; selected from the same differential-expression table',
      x = NULL, y = 'Normalized log2 expression', fill = 'Group', color = 'Group') +
    plot_theme() + ggplot2::theme(legend.position = 'none', panel.grid.major.x = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = 'bold'), axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
}

# Unequal row heights + a full-width lower panel. Heatmap legends stay with each
# heatmap because Group/Batch/Row Z/Pearson r encode different variables.
compose_figure <- function(pca, correlation, volcano, heatmap, expression_plot,
                           n_genes, n_samples) {
  correlation_panel <- patchwork::wrap_elements(full = correlation)
  heatmap_panel <- patchwork::wrap_elements(full = heatmap)
  figure <- patchwork::wrap_plots(A = pca, B = correlation_panel, C = volcano,
    D = heatmap_panel, E = expression_plot,
    design = 'AB\nCD\nEE', heights = c(1, 1.15, .75), widths = c(1, 1.12)) +
    patchwork::plot_annotation(title = 'EXPRESSION ATLAS',
      subtitle = sprintf('Synthetic demonstration  /  %s genes  /  %s samples  /  one connected dataset', n_genes, n_samples),
      caption = 'SIMULATED DATA ONLY  |  Welch tests on log2 values + BH adjustment; not an RNA-seq inference workflow.\nHeatmap genes and representative genes are selected from this dataset, not independently validated.',
      tag_levels = 'A',
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 23, face = 'bold', color = '#25354A'),
        plot.subtitle = ggplot2::element_text(size = 11, color = '#66768A'),
        plot.caption = ggplot2::element_text(size = 9, color = '#66768A', hjust = 0),
        plot.background = ggplot2::element_rect(fill = 'white', color = NA),
        plot.margin = ggplot2::margin(18, 18, 14, 18)))
  figure & ggplot2::theme(plot.tag = ggplot2::element_text(size = 16, face = 'bold', color = '#25354A'))
}
# Description, GeneRatio (numeric or 'k/n'), Count, p.adjust.
# This visualizes an existing enrichment table; it does not run enrichment tests.
plot_enrichment <- function(results, top_n = 12) {
  check_columns(results, c('Description', 'GeneRatio', 'Count', 'p.adjust'))
  check_probability(results$p.adjust, 'p.adjust')
  if (length(top_n) != 1 || !is.finite(top_n) || top_n < 1 || top_n != floor(top_n)) stop('top_n must be a positive integer.')
  d <- results
  if (anyNA(d$Description) || anyDuplicated(d$Description) || any(!nzchar(as.character(d$Description)))) stop('Description must be unique and nonempty.')
  if (!is.numeric(d$GeneRatio)) {
    ratios <- strsplit(as.character(d$GeneRatio), '/', fixed = TRUE)
    d$GeneRatio <- vapply(ratios, function(x) {
      if (length(x) != 2) stop('GeneRatio must be numeric or k/n.')
      v <- suppressWarnings(as.numeric(x))
      if (any(!is.finite(v)) || v[2] <= 0 || v[1] < 0 || any(v != floor(v))) stop('Invalid k/n ratio.')
      v[1] / v[2]
    }, numeric(1))
  }
  check_probability(d$GeneRatio, 'GeneRatio')
  if (!is.numeric(d$Count) || any(!is.finite(d$Count) | d$Count <= 0 | d$Count != floor(d$Count))) stop('Count must be positive integers.')
  d <- head(d[order(d$p.adjust, -d$Count, d$Description), ], top_n)
  d$Description <- factor(d$Description, levels = rev(d$Description))
  ggplot2::ggplot(d, ggplot2::aes(GeneRatio, Description, size = Count, color = p.adjust)) +
    ggplot2::geom_point(alpha = 1) + ggplot2::scale_size_area(max_size = 9) +
    ggplot2::scale_color_gradient(low = '#B2182B', high = '#2166AC', name = 'Adjusted P') +
    ggplot2::scale_y_discrete(labels = function(x) vapply(x, function(s) paste(strwrap(s, 38), collapse = '\n'), character(1))) +
    ggplot2::labs(title = 'Enrichment results', subtitle = 'Top terms ranked by adjusted P value', x = 'Gene ratio', y = NULL) + plot_theme()
}
