#' @title Project Path Utilities
#' @description Helper functions for constructing paths relative to the project root
#'   and creating output directory structures for WHAM / MSE projections.
#' @details Provides path normalization ensuring consistency across platforms and operating environments.
#' @name project_paths
NULL

#' Construct Absolute Project Path
#'
#' @description Generates a normalized absolute file path relative to the current project directory.
#' @param ... Optional character components or subdirectory names to append to the root path.
#' @return A character string representing the normalized absolute file path.
#' @export
project_path <- function(...) {
  # Normalize current working directory path to resolve symlinks and relative tokens
  base_dir <- normalizePath(getwd(), mustWork = TRUE)
  # Append any provided sub-paths or file names safely
  file.path(base_dir, ...)
}

#' Create Project Output Directory
#'
#' @description Constructs a directory path within the project's `output/` folder and creates
#'   the directory if it does not already exist.
#' @param ... Character strings representing subdirectories within `output/`.
#' @return A character string with the absolute path to the initialized output directory.
#' @export
create_output_dir <- function(...) {
  # Resolve path under project output directory
  output_dir <- project_path("output", ...)
  # Recursively create target folder, suppressing warnings if directory exists
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir
}
