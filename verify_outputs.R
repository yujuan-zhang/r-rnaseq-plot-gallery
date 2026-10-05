args <- commandArgs()
script <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
root <- dirname(normalizePath(script))
names <- c("pca","correlation","volcano","heatmap","gene_expression","combined")
files <- file.path(root, "figures", as.vector(outer(names, c(".png", ".svg"), paste0)))
stopifnot(all(file.exists(files)), all(file.info(files)$size > 0))
for (file in files[grepl("\\.png$", files)]) {
  con <- file(file, "rb"); signature <- readBin(con, "raw", 8); close(con)
  stopifnot(identical(signature, as.raw(c(137,80,78,71,13,10,26,10))))
}
for (file in files[grepl("\\.svg$", files)])
  stopifnot(any(grepl("<svg", readLines(file, warn = FALSE), fixed = TRUE)))
message("Verified all 12 non-empty PNG/SVG outputs.")
