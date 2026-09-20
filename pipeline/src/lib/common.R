# ==============================================================================
# src/lib/common.R
# Small shared helpers (package loading, manuscript checkpoints).
# All code comments in this repository are in English.
# ==============================================================================

# Attach a set of packages, stopping with an informative error message if any
# package is missing. Mirrors the library() calls that the original analysis
# scripts performed at the top of each script.
require_pkgs <- function(pkgs) {
    missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, logical(1))]
    if (length(missing) > 0L) {
        stop(
            "The following required R packages are not installed: ",
            paste(missing, collapse = ", "),
            ". Install them before running this script."
        )
    }
    invisible(lapply(pkgs, library, character.only = TRUE))
}

# Informative sanity check against numbers reported in the manuscript.
# The check NEVER aborts the pipeline; it only prints a warning. Turn it off by
# setting CHECK_MANUSCRIPT <- FALSE in config/params.R.
check_manuscript <- function(metric, observed, expected) {
    if (!isTRUE(CHECK_MANUSCRIPT)) return(invisible(NULL))
    if (!identical(as.integer(observed), as.integer(expected))) {
        warning(
            sprintf("Manuscript checkpoint mismatch: %s observed = %s, expected = %s.",
                    metric, observed, expected),
            immediate. = TRUE
        )
    } else {
        message(sprintf("Manuscript checkpoint OK: %s = %s.", metric, observed))
    }
    invisible(NULL)
}
