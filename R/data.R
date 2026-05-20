# =============================================================================
# data.R — Lazy-loader for the toy CTBN dataset shipped with the package.
#
# The dataset is stored as a gzipped CSV in inst/extdata/toy_DT_wide.csv.gz
# (so it is human-inspectable). load_toy_DT_wide() reads it into a
# data.table on demand. If the corresponding .rds is present (created by
# data-raw/make_toy_DT_wide.R) it is preferred for speed.
# =============================================================================

#' Load the toy CTBN panel dataset shipped with the package
#'
#' @return A \code{data.table} with the toy panel data; see
#'   \code{\link{toy_DT_wide}} for the column dictionary.
#' @export
load_toy_DT_wide <- function() {
  rds <- system.file("extdata", "toy_DT_wide.rds",    package = "CTBNMultimorbidity")
  csv <- system.file("extdata", "toy_DT_wide.csv.gz", package = "CTBNMultimorbidity")

  # Stand-alone (non-installed) fallback paths used when sourcing R/ directly
  if (rds == "")    rds <- "inst/extdata/toy_DT_wide.rds"
  if (csv == "")    csv <- "inst/extdata/toy_DT_wide.csv.gz"

  if (file.exists(rds)) return(readRDS(rds))
  if (file.exists(csv)) {
    if (!requireNamespace("data.table", quietly = TRUE))
      stop("Package 'data.table' is required to load the toy dataset.")
    return(data.table::fread(csv))
  }
  stop("Toy dataset not found. Re-run data-raw/make_toy_DT_wide.R or ",
       "install the package.")
}
