#' Add count data to a PSI SummarizedExperiment object
#'
#' Adds two assays ('counts' and 'other_counts') to a PSI SummarizedExperiment
#' object, making it suitable for downstream analysis using the DEXSeq package.
#'
#' @param se_psi A PSI SummarizedExperiment object returned by the
#' \code{\link{transcript_to_psi}} function.
#' @param se_gene A gene-level SummarizedExperiment object returned by the
#' \code{\link{tcc_to_gene}} function. It must be compatible with the se_psi
#' object (i.e. they must originate from the same TCC data).
#' @return A copy of the se_psi PSI SummarizedExperiment object with count
#' assays added.
#' @export
add_psi_counts <- function(se_psi,
                           se_gene) {

    # Check arguments
    assertthat::assert_that(methods::is(se_psi, "SummarizedExperiment"))
    assertthat::assert_that(is.element(
        "psi", SummarizedExperiment::assayNames(se_psi)
    ))
    assertthat::assert_that(is.element(
        "gene_id", colnames(SummarizedExperiment::rowData(se_psi))
    ))
    assertthat::assert_that(methods::is(se_gene, "SummarizedExperiment"))
    assertthat::assert_that(is.element(
        "counts", SummarizedExperiment::assayNames(se_gene)
    ))
    assertthat::assert_that(is.element(
        "gene_id", colnames(SummarizedExperiment::rowData(se_gene))
    ))
    assertthat::assert_that(identical(ncol(se_psi), ncol(se_gene)))

    # Prepare the counts assay
    psi_values <- SummarizedExperiment::assay(se_psi, "psi")
    is_bulk_rnaseq <- is.matrix(psi_values)
    gene_counts <- SummarizedExperiment::assay(se_gene, "counts")[
        SummarizedExperiment::rowData(se_psi)$gene_id,
    ]
    psi_counts <- psi_values * gene_counts
    if (!is_bulk_rnaseq) {
        psi_counts <- methods::as(psi_counts, "dgCMatrix")
    }

    # Prepare the other_counts assay
    psi_other_values <- 1 - psi_values
    psi_other_counts <- psi_other_values * gene_counts
    if (!is_bulk_rnaseq) {
        psi_other_counts <- methods::as(psi_other_counts, "dgCMatrix")
    }

    # Add count assays to the PSI SE object
    se <- se_psi
    SummarizedExperiment::assay(se, "counts") <- psi_counts
    SummarizedExperiment::assay(se, "other_counts") <- psi_other_counts

    return(se)
}
