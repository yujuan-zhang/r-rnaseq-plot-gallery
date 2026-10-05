# RNA-seq Visualization in R

## What it does

Reusable R templates for visualizing RNA-seq results and creating publication-style multi-panel figures. Five coordinated views connect sample structure, differential expression and gene-level patterns. The included example uses simulated log2 expression data.

![Gene expression figure](figures/combined.png)

## Input

### Use your own data

| File | Required structure |
|---|---|
| `expression.csv` | First column: unique gene IDs. Remaining columns: normalized log2 expression, one sample per column. |
| `samples.csv` | `sample`, `group`; optional `batch`. Sample IDs must match the expression matrix. |
| `differential_expression.csv` | Unique `gene`, finite `log2FC`, and `padj` in [0, 1]. Gene IDs must match the matrix for the full workflow. |

Adjust thresholds, selected gene counts and resolution in `run.R`. Positive log2FC is labelled as Treatment versus Control; update the axis label for other contrasts. Inputs must not contain NA or Inf.

To reuse an individual plot:

```r
source('plots.R')
results <- read.csv('data/differential_expression.csv')
p <- plot_volcano(results, lfc_cutoff = 1, padj_cutoff = 0.05, label_n = 6)
ggplot2::ggsave('figures/custom_volcano.svg', p,
                device = svglite::svglite, width = 9, height = 6.5)
```

PCA, volcano and gene-expression functions return ggplot objects. Heatmap functions return grid grobs, which support `grid::grid.draw()` and `ggplot2::ggsave()`.

Change the panel arrangement in `compose_figure()`:

```r
design = 'AB\nCD\nEE'      # E spans both columns
heights = c(1, 1.15, .75)
widths = c(1, 1.12)
```

Heatmap legends remain within their panels. Recheck text and legend placement after resizing the composite.

## Output

PCA, correlation, volcano, heatmap, gene-expression and combined figures in `figures/`, each as PNG and SVG.

## Try it

### Quick start

Install dependencies in R:

```r
install.packages(c('ggplot2', 'ggrepel', 'patchwork', 'ragg', 'svglite', 'circlize'),
                 repos = 'https://cloud.r-project.org')
install.packages('BiocManager')
BiocManager::install('ComplexHeatmap', ask = FALSE, update = FALSE)
```

From the project directory, run:

```r
source('run.R')
```

Alternatively, use `Rscript run.R`. Outputs are saved to `figures/` as `pca`, `correlation`, `volcano`, `heatmap`, `gene_expression` and `combined`, each in PNG and SVG format.

Tested on macOS with R 4.6.1. Tested package versions and installation commands are listed at the top of `run.R`.

### Verify the example

The bundled example completed in **16.62 seconds** on an Intel macOS machine with R 4.6.1 (800 synthetic genes and 24 samples); installation is excluded. This is a measured example, not a runtime guarantee. No additional data download is needed and normal runs preserve the input CSVs.

```bash
Rscript verify_outputs.R
```

### Features

- PCA with sample labels, group ellipses and batch-specific shapes.
- Clustered correlation and expression heatmaps with sample annotations, expression modules and side bar charts.
- Volcano plots with configurable thresholds, gene labels and highlighted genes linked to expression panels.
- Violin, box and sample-point layers for representative genes.
- A–E panel labels, unequal row heights and a full-width lower panel using patchwork.
- PNG and editable SVG export; input validation for identifiers, missing values and degenerate data.

### Project structure

```text
run.R       Configuration and figure generation
plots.R     Reusable plotting and simulation functions
data/       Expression matrix, sample metadata and differential results
figures/    Individual plots and composite figure
README.md   Usage and methods
```

### Example data and methods

The example contains **800 synthetic genes and 24 samples**, split across two groups and two balanced batches. Differential results use Welch t-tests on the same simulated log2 matrix, followed by BH adjustment across all genes. This demonstrates visualization, not a count-based RNA-seq analysis or batch-adjusted inference.

- **PCA:** centered, nonconstant genes without unit-variance scaling. The 80% normal-approximation ellipses describe score distributions, not confidence intervals for group means; undersized or singular groups are omitted.
- **Correlation:** Pearson r, average-linkage clustering on `1 − r`. The labelled color scale adapts to the observed range.
- **Expression heatmap:** top 24 qualifying genes ranked by adjusted P value; row Z-scores, Euclidean distance and complete linkage. Two clusters represent expression patterns, not pathways. Colors saturate at ±2; side bars show mean log2 expression. Samples cluster within groups.
- **Representative genes:** up to two upregulated and two downregulated genes, with fallback to other qualifying genes. Facets have independent y-axis scales; selection is exploratory, not independent validation.

Set `regenerate_example <- TRUE` only to rebuild the simulated inputs; this overwrites the three example CSVs. Normal runs leave inputs unchanged.

### References

[ggplot2](https://ggplot2.tidyverse.org/) · [ggrepel](https://ggrepel.slowkow.com/) · [ComplexHeatmap](https://jokergoo.github.io/ComplexHeatmap-reference/book/) · [patchwork](https://patchwork.data-imaginist.com/)

