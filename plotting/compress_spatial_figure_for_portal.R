#!/usr/bin/env Rscript

# Re-encode the publication PNG as a compact, page-sized PDF without altering
# the pixels or scientific content. This keeps the portal/source bundle small
# while preserving the high-resolution PNG as the archival figure.

suppressPackageStartupMessages({
  library(png)
  library(grid)
})

input <- file.path("figures", "final", "Figure5.png")
output <- file.path("figures", "final", "Figure5.pdf")

if (!file.exists(input)) {
  stop("Missing input figure: ", input)
}

img <- readPNG(input)
dims <- dim(img)
aspect <- dims[2] / dims[1]
width_in <- 7.4
height_in <- width_in / aspect

cairo_pdf(output, width = width_in, height = height_in, onefile = TRUE)
grid.newpage()
grid.raster(img, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)
dev.off()

message(normalizePath(output, winslash = "/"))
