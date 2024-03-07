#' Filter PSI events
#'
#' Prepares a vector of filtered PSI events suitable to be used with the
#' \code{\link{psi_to_dexseq}} function.
#'
#' @param se_psi A PSI SummarizedExperiment object returned by a function from
#' the \code{\link{Isosceles-package}}.
#' @param gene_ids A character vector specifying the gene IDs to restrict
#' the results to (ignored if set to NULL).
#' @param exclude_tss A logical scalar specifying if TSS (transcription start
#' sites) regions should be excluded from the results.
#' @param exclude_tes A logical scalar specifying if TES (transcription end
#' sites) regions should be excluded from the results.
#' @param min_mean_psi A numeric scalar specifying the minimum value threshold
#' for mean PSI values across all cells.
#' @param max_mean_psi A numeric scalar specifying the maximum value threshold
#' for mean PSI values across all cells.
#' @param min_alt_psi_count An integer scalar specifying the required minimum
#' number of cells with PSI values different than 0, 0.5 or 1.
#' @param min_incl_psi A numeric scalar specifying the minimum PSI value for
#' cells considered to include given PSI event.
#' @param min_incl_psi_count An integer scalar specifying the required minimum
#' number of cells including given PSI event.
#' @return A vector of filtered PSI event identifiers.
#' @export
filter_psi_events <- function(se_psi,
                              gene_ids = NULL,
                              exclude_tss = TRUE,
                              exclude_tes = TRUE,
                              min_mean_psi = 0.025,
                              max_mean_psi = 0.975,
                              min_alt_psi_count = 30,
                              min_incl_psi = 0.1,
                              min_incl_psi_count = 30) {

    # Check arguments
    assertthat::assert_that(methods::is(se_psi, "SummarizedExperiment"))
    assertthat::assert_that(is.element(
        "psi", SummarizedExperiment::assayNames(se_psi)
    ))
    assertthat::assert_that(is.element(
        "gene_id", colnames(SummarizedExperiment::rowData(se_psi))
    ))
    if (!is.null(gene_ids)) {
        assertthat::assert_that(is.character(gene_ids))
    }
    assertthat::assert_that(assertthat::is.flag(exclude_tss))
    assertthat::assert_that(assertthat::is.flag(exclude_tes))
    assertthat::assert_that(assertthat::is.number(min_mean_psi))
    assertthat::assert_that(min_mean_psi >= 0)
    assertthat::assert_that(min_mean_psi <= 1)
    assertthat::assert_that(assertthat::is.number(max_mean_psi))
    assertthat::assert_that(max_mean_psi >= 0)
    assertthat::assert_that(max_mean_psi <= 1)
    assertthat::assert_that(assertthat::is.count(min_alt_psi_count))
    assertthat::assert_that(assertthat::is.number(min_incl_psi))
    assertthat::assert_that(min_incl_psi >= 0)
    assertthat::assert_that(min_incl_psi <= 1)
    assertthat::assert_that(assertthat::is.count(min_incl_psi_count))

    # Prepare filtered PSI events of interest
    psi_assay <- SummarizedExperiment::assay(se_psi, "psi")
    if (!is.null(gene_ids)) {
        psi_assay <- psi_assay[
            SummarizedExperiment::rowData(se_psi)$gene_id %in% gene_ids,,drop = FALSE
        ]

    }
    if (exclude_tss) {
        psi_assay <- psi_assay[
            !grepl(":TSS", rownames(psi_assay)),,drop = FALSE
        ]
    }
    if (exclude_tes) {
        psi_assay <- psi_assay[
            !grepl(":TES", rownames(psi_assay)),,drop = FALSE
        ]
    }
    psi_selector <- (apply(psi_assay, 1, mean) >= min_mean_psi) &
        (apply(psi_assay, 1, mean) <= max_mean_psi) &
        (apply(psi_assay, 1, function(x) { sum((x > 0) & (x < 0.999) & (x != 0.5)) }) >= min_alt_psi_count) &
        (apply(psi_assay, 1, function(x) { sum(x >= min_incl_psi) }) >= min_incl_psi_count)
    psi_assay <- psi_assay[psi_selector,,drop = FALSE]
    psi_events <- rownames(psi_assay)

    return(psi_events)
}
