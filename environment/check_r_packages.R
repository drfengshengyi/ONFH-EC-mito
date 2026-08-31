#!/usr/bin/env Rscript

# Verify the exact R runtime and directly used package versions archived for
# the manuscript release.  This script does not install or update packages.

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
script_file <- if (length(script_arg)) {
  sub("^--file=", "", script_arg[[1]])
} else {
  "environment/check_r_packages.R"
}
root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
manifest_path <- file.path(root, "environment", "r-package-versions.tsv")
manifest <- read.delim(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)

problems <- character(0)
expected_r <- manifest$version[manifest$component == "R"]
if (length(expected_r) != 1L || as.character(getRversion()) != expected_r) {
  problems <- c(
    problems,
    sprintf("R: expected %s, found %s", paste(expected_r, collapse = "/"), getRversion())
  )
}

packages <- manifest[manifest$component != "R", , drop = FALSE]
for (index in seq_len(nrow(packages))) {
  package <- packages$component[[index]]
  expected <- packages$version[[index]]
  if (!requireNamespace(package, quietly = TRUE)) {
    problems <- c(problems, sprintf("%s: missing (expected %s)", package, expected))
    next
  }
  observed <- as.character(utils::packageVersion(package))
  if (!identical(observed, expected)) {
    problems <- c(problems, sprintf("%s: expected %s, found %s", package, expected, observed))
  }
}

if (length(problems)) {
  cat(paste0("FAIL  ", problems, "\n"), sep = "")
  stop("R release environment does not match environment/r-package-versions.tsv", call. = FALSE)
}

cat(sprintf("PASS  R %s and %d package versions match the release manifest.\n",
            getRversion(), nrow(packages)))
