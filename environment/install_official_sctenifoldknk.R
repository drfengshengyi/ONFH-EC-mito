#!/usr/bin/env Rscript

# Install the authors' official scTenifoldKnk implementation and the small
# set of packages used by virtual_knockout/run_official_vko.R.

options(
  repos = c(
    CaiLab = "https://cailab-tamu.r-universe.dev",
    CRAN = "https://cloud.r-project.org"
  )
)

packages <- c(
  "scTenifoldNet", "scTenifoldKnk", "Matrix", "jsonlite",
  "ggplot2", "patchwork", "ggrepel", "dplyr", "tidyr", "ragg", "svglite",
  "cowplot", "tibble"
)
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing) > 0L) {
  install.packages(missing, dependencies = TRUE)
}

still_missing <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(still_missing) > 0L) {
  stop(
    "Installation failed for: ", paste(still_missing, collapse = ", "),
    call. = FALSE
  )
}

message("Installed package versions:")
for (package in packages) {
  message("  ", package, " ", as.character(utils::packageVersion(package)))
}
