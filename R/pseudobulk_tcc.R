#' Prepare a pseudobulk TCC SummarizedExperiment object
#'
#' Prepares a pseudobulk TCC SummarizedExperiment from TCC data and given cell
#' labels.
#'
#' @param se_tcc A TCC SummarizedExperiment object returned by the
#' \code{\link{prepare_tcc_se}} function.
#' @param cell_labels A vector or a factor containing cell labels acting as a
#' grouping variable.
#' @return A pseudobulk SummarizedExperiment object containing TCC annotation
#' and quantification data.
#' @export
pseudobulk_tcc <- function(se_tcc,
                           cell_labels) {

    # Check arguments
    assertthat::assert_that(methods::is(se_tcc, "SummarizedExperiment"))
    assertthat::assert_that(is.element(
        "counts", SummarizedExperiment::assayNames(se_tcc)
    ))
    assertthat::assert_that(length(cell_labels) == ncol(se_tcc))

    # Sum raw count values across cells based on their labels
    merged_count_matrix <- scuttle::sumCountsAcrossCells(
        SummarizedExperiment::assay(se_tcc, "counts"),
        cell_labels
    ) %>% SummarizedExperiment::assay("sum")

    # Create the SummarizedExperiment object containing pseudo-bulk expression data
    se_pseudobulk <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = merged_count_matrix)
    )
    SummarizedExperiment::rowRanges(se_pseudobulk) <-
        SummarizedExperiment::rowRanges(se_tcc)
    SummarizedExperiment::rowData(se_pseudobulk) <-
        SummarizedExperiment::rowData(se_tcc)
    S4Vectors::metadata(se_pseudobulk) <- S4Vectors::metadata(se_tcc)
    SummarizedExperiment::assay(se_pseudobulk, "tpm") <-
        calculate_tpm(se_pseudobulk)
    SummarizedExperiment::assay(se_pseudobulk, "relative_expression") <-
        calculate_relative_expression(se_pseudobulk)

    return(se_pseudobulk)
}
