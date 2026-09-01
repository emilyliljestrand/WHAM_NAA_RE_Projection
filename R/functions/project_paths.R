# Helpers for paths rooted at the project directory.

project_path <- function(...) {
  file.path(normalizePath(getwd(), mustWork = TRUE), ...)
}

create_output_dir <- function(...) {
  output_dir <- project_path("output", ...)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir
}