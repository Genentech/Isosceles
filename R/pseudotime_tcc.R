#' Merge TCC values using moving window over pseudotime
#'
#' Prepares a pseudotime window TCC SummarizedExperiment from TCC data and
#' pseudotime values.
#'
#' @param se_tcc A TCC SummarizedExperiment object returned by the
#' \code{\link{bam_to_tcc}} function.
#' @param pseudotime A numeric vector containing the pseudotime values for
#' each cell. Cells not belonging to the analyzed trajectory should be denoted
#' using NA values.
#' @param trim A numeric scalar specifying the fraction (0 to 0.5) of cells
#' to be trimmed from each end of the pseudotime spectrum.
#' @param window_size An integer scalar specifying the window size.
#' @param window_step An integer scalar specifying the window step.
#' @return A SummarizedExperiment object containing TCC data for pseudotime
#' windows.
#' @export
pseudotime_tcc <- function(se_tcc,
                           pseudotime,
                           trim = 0,
                           window_size = 30,
                           window_step = 15) {

    # Check arguments
    assertthat::assert_that(methods::is(se_tcc, "SummarizedExperiment"))
    assertthat::assert_that(is.element(
        "counts", SummarizedExperiment::assayNames(se_tcc)
    ))
    assertthat::assert_that(is.numeric(pseudotime))
    assertthat::assert_that(identical(length(pseudotime), ncol(se_tcc)))
    assertthat::assert_that(assertthat::is.number(trim))
    assertthat::assert_that(trim >= 0)
    assertthat::assert_that(trim < 0.5)
    assertthat::assert_that(assertthat::is.count(window_size))
    assertthat::assert_that(window_size <= ncol(se_tcc))
    assertthat::assert_that(assertthat::is.count(window_step))
    assertthat::assert_that(window_step <= window_size)

    # Remove cells with smallest & highest pseudotime values
    quantile_low <- stats::quantile(pseudotime, trim, na.rm = TRUE)
    quantile_high <- stats::quantile(pseudotime, 1 - trim, na.rm = TRUE)
    pseudotime[pseudotime < quantile_low] <- NA
    pseudotime[pseudotime > quantile_high] <- NA

    # Arrange cells by pseudotime
    pseudotime <- stats::setNames(pseudotime, colnames(se_tcc))
    pseudotime <- pseudotime[!is.na(pseudotime)]
    pseudotime <- sort(pseudotime)
    se_tcc <- se_tcc[, names(pseudotime)]

    # Prepare TCC data
    tcc_counts <- SummarizedExperiment::assay(se_tcc, "counts")
    is_bulk_rnaseq <- is.matrix(tcc_counts)

    # Merge data using a sliding window
    num_windows <- ceiling((1 + length(pseudotime) - window_size) / window_step)
    window_positions <- lapply(seq(num_windows), function(x) {
        c((x - 1) * window_step + 1, (x - 1) * window_step + window_size)
    })
    pseudotime_avg <- sapply(window_positions, function(x) {
        mean(pseudotime[x[1]:x[2]])
    })
    merged_tcc_counts <- lapply(window_positions, function(x) {
        window_counts <- Matrix::rowSums(tcc_counts[, x[1]:x[2]])
        window_counts <- as.matrix(window_counts)
        if (!is_bulk_rnaseq) {
            window_counts <- methods::as(window_counts, "dgCMatrix")
        }
    })
    merged_tcc_counts <- do.call(cbind, merged_tcc_counts)
    colnames(merged_tcc_counts) <- paste0(
        "Cells_", sapply(window_positions, paste0, collapse = "_")
    )

    # Prepare the SummarizedExperiment object
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = merged_tcc_counts),
        rowData = SummarizedExperiment::rowData(se_tcc)
    )
    SummarizedExperiment::colData(se)$pseudotime <- pseudotime_avg
    S4Vectors::metadata(se) <- S4Vectors::metadata(se_tcc)
    SummarizedExperiment::assay(se, "tpm") <- calculate_tpm(se)
    SummarizedExperiment::assay(se, "relative_expression") <-
        calculate_relative_expression(se)

    return(se)
}
