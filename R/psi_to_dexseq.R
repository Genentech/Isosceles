#' Prepare a PSI count DEXSeqDataSet object
#'
#' Prepares a PSI count DEXSeqDataSet object suitable for the analysis of PSI
#' count changes along the given variable(s) (categorical or continuous).
#'
#' @param se_psi A PSI SummarizedExperiment object returned by the
#' \code{\link{add_psi_counts}} function.
#' @param condition A vector or a factor containing the condition variable
#' (categorical or continuous) used in the design formula. Alternatively,
#' a data frame containing multiple variables in separate columns can be used,
#' in which case the design formula needs to be adjusted by the user.
#' @param design A formula which specifies the design of the experiment. See
#' the DEXSeq package documentation for more information.
#' @param psi_events A character vector specifying the PSI events to restrict
#' the analysis to (ignored if set to NULL).
#' @param remove_redundant_psi A logical scalar specifying if PSI events with
#' redundant count profiles should be removed from the analysis.
#' @return A DEXSeqDataSet object containing PSI count data, suitable for
#' further analysis using the DEXSeq package.
#' @export
psi_to_dexseq <- function(se_psi,
                          condition,
                          design = ~ sample + exon + condition:exon,
                          psi_events = NULL,
                          remove_redundant_psi = TRUE) {

    # Check arguments
    assertthat::assert_that(methods::is(se_psi, "SummarizedExperiment"))
    assertthat::assert_that(is.element(
        "counts", SummarizedExperiment::assayNames(se_psi)
    ))
    assertthat::assert_that(is.element(
        "other_counts", SummarizedExperiment::assayNames(se_psi)
    ))
    if (!is.data.frame(condition)) {
        assertthat::assert_that(length(condition) == ncol(se_psi))
    }
    if (is.data.frame(condition)) {
        assertthat::assert_that(nrow(condition) == ncol(se_psi))
    }
    assertthat::assert_that(methods::is(design, "formula"))
    if (!is.null(psi_events)) {
        assertthat::assert_that(is.character(psi_events))
        assertthat::assert_that(length(psi_events) > 1)
    }
    assertthat::assert_that(assertthat::is.flag(remove_redundant_psi))

    # Prepare PSI count data
    if (is.null(psi_events)) {
        psi_events <- rownames(se_psi)
    }
    psi_counts <- as.matrix(SummarizedExperiment::assay(
        se_psi[psi_events,], "counts"
    ))
    psi_other_counts <- as.matrix(SummarizedExperiment::assay(
        se_psi[psi_events,], "other_counts"
    ))

    # Remove PSI events with redundant count profiles if required
    if (remove_redundant_psi) {
        row_gene_ids <- sapply(strsplit(rownames(psi_counts), ":"), "[", 1)
        row_hash_ids <- apply(psi_counts, 1, rlang::hash)
        row_ids <- paste0(row_gene_ids, ".", row_hash_ids)
        row_selector <- !duplicated(row_ids)
        psi_counts <- psi_counts[row_selector,]
        psi_other_counts <- psi_other_counts[row_selector,]
    }

    # Prepare the DEXSeqDataSet object
    psi_counts <- round(psi_counts)
    psi_other_counts <- round(psi_other_counts)
    rownames(psi_counts) <- gsub(":", "|", rownames(psi_counts))
    rownames(psi_other_counts) <- gsub(":", "|", rownames(psi_other_counts))
    col_data <- data.frame(
        se_psi_colname = factor(colnames(psi_counts),
                                levels = colnames(psi_counts))
    )
    if (!is.data.frame(condition)) {
        col_data$condition <- condition
    } else {
        col_data <- cbind(col_data, condition)
    }
    dxd <- DEXSeq::DEXSeqDataSet(
        countData = psi_counts,
        alternativeCountData = psi_other_counts,
        sampleData = col_data,
        design = design,
        featureID = rownames(psi_counts),
        groupID = rownames(psi_counts)
    )
    SummarizedExperiment::rowData(dxd)$featureID <-
        gsub("\\|", ":", SummarizedExperiment::rowData(dxd)$featureID)
    SummarizedExperiment::rowData(dxd)$groupID <-
        gsub("\\|", ":", SummarizedExperiment::rowData(dxd)$groupID)
    rownames(dxd) <- SummarizedExperiment::rowData(dxd)$featureID

    return(dxd)
}
