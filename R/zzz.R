# This script is from the tidyverse package
load <- c("ggplot2", "tidyr", "readr", "dplyr", "qtl", "rrBLUP", "agricolae", "lattice")

# Function to run when the package is attached
.onAttach <- function(...) {

  # Find the packages that aren't already attached
  needed <- load[!is_attached(load)]

  if (length(needed) == 0)
    return()

  packageStartupMessage(paste0("Loading pbr: ", needed, collapse = "\n"))
  suppressPackageStartupMessages(
    lapply(needed, library, character.only = TRUE, warn.conflicts = FALSE)
  )

}

is_attached <- function(x) {
  paste0("package:", x) %in% search()
}
