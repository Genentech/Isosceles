#' Prepare a PSI count DEXSeqDataSet object
#'
#' Aggregates TCC values using pseudotime windows and creates a DEXSeqDataSet
#' object suitable for the analysis of PSI count changes along given pseudotime
#' trajectory.
#'
#' @param se_tcc A TCC SummarizedExperiment object returned by the
#' \code{\link{bam_to_tcc}} function.
#' @param pseudotime A numeric vector containing the pseudotime values for
#' each cell. Cells not belonging to the analyzed trajectory should be denoted
#' using NA values.
#' @param psi_events A character vector specifying the PSI events to restrict
#' the analysis to (ignored if set to NULL).
#' @param trim A numeric scalar specifying the fraction (0 to 0.5) of cells
#' to be trimmed from each end of the pseudotime spectrum.
#' @param window_size An integer scalar specifying the window size.
#' @param window_step An integer scalar specifying the window step.
#' @param remove_redundant_psi A logical scalar specifying if PSI events with
#' redundant count profiles should be removed from the analysis.
#' @param scale_pseudotime A logical scalar specifying if pseudotime values for
#' the windows should be scaled.
#' @param ncpu An integer scalar specifying the number of cores to use for
#' multicore parallelization.
#' @return A DEXSeqDataSet object containing PSI count data for pseudotime
#' windows, suitabe for further analysis using the DEXSeq package.
#' @export
dexseq_psi <- function(se_tcc,
                       pseudotime,
                       psi_events = NULL,
                       trim = 0,
                       window_size = 30,
                       window_step = 15,
                       remove_redundant_psi = TRUE,
                       scale_pseudotime = TRUE,
                       ncpu = 1) {

    # Check arguments
    assertthat::assert_that(methods::is(se_tcc, "SummarizedExperiment"))
    assertthat::assert_that(is.element(
        "counts", SummarizedExperiment::assayNames(se_tcc)
    ))
    assertthat::assert_that(is.numeric(pseudotime))
    assertthat::assert_that(identical(length(pseudotime), ncol(se_tcc)))
    if (!is.null(psi_events)) {
        assertthat::assert_that(is.character(psi_events))
        assertthat::assert_that(length(psi_events) > 1)
    }
    assertthat::assert_that(assertthat::is.number(trim))
    assertthat::assert_that(trim >= 0)
    assertthat::assert_that(trim < 0.5)
    assertthat::assert_that(assertthat::is.count(window_size))
    assertthat::assert_that(assertthat::is.count(window_step))
    assertthat::assert_that(assertthat::is.flag(remove_redundant_psi))
    assertthat::assert_that(assertthat::is.flag(scale_pseudotime))
    assertthat::assert_that(assertthat::is.count(ncpu))

    # Prepare pseudotime windows SE objects
    se_window_tcc <- pseudotime_tcc(
        se_tcc = se_tcc,
        pseudotime = pseudotime,
        trim = trim,
        window_size = window_size,
        window_step = window_step
    )
    se_window_gene <- tcc_to_gene(
        se_tcc = se_window_tcc
    )
    se_window_transcript <- tcc_to_transcript(
        se_tcc = se_window_tcc,
        use_length_normalization = FALSE,
        ncpu = ncpu
    )
    se_window_psi <- transcript_to_psi(
        se = se_window_transcript,
        ncpu = ncpu
    )

    # Prepare PSI count data
    if (is.null(psi_events)) {
        psi_events <- rownames(se_window_psi)
    }
    psi_values <- SummarizedExperiment::assay(
        se_window_psi[psi_events,], "psi"
    )
    psi_values <- as.matrix(psi_values)
    gene_counts <- SummarizedExperiment::assay(se_window_gene, "counts")[
        SummarizedExperiment::rowData(se_window_psi[psi_events,])$gene_id,
    ]
    gene_counts <- as.matrix(gene_counts)
    psi_counts <- psi_values * gene_counts

    # Remove PSI events with redundant count profiles if required
    if (remove_redundant_psi) {
        row_gene_ids <- sapply(strsplit(rownames(psi_counts), ":"), "[", 1)
        row_hash_ids <- apply(psi_counts, 1, rlang::hash)
        row_ids <- paste0(row_gene_ids, ".", row_hash_ids)
        row_selector <- !duplicated(row_ids)
        psi_counts <- psi_counts[row_selector,]
    }

    # Prepare PSI other count data
    psi_other_values <- 1 - psi_values
    psi_other_counts <- psi_other_values * gene_counts
    psi_other_counts <- psi_other_counts[rownames(psi_counts),]

    # Prepare the DEXSeqDataSet object
    psi_counts <- round(psi_counts)
    psi_other_counts <- round(psi_other_counts)
    rownames(psi_counts) <- gsub(":", "|", rownames(psi_counts))
    rownames(psi_other_counts) <- gsub(":", "|", rownames(psi_other_counts))
    window_pseudotime <- se_window_tcc$pseudotime
    if (scale_pseudotime) {
        window_pseudotime <- scale(window_pseudotime)
    }
    col_data <- data.frame(
        window_name = factor(colnames(psi_counts),
                             levels = colnames(psi_counts)),
        pseudotime = window_pseudotime
    )
    dxd <- DEXSeq::DEXSeqDataSet(
        countData = psi_counts,
        alternativeCountData = psi_other_counts,
        sampleData = col_data,
        design = ~sample + exon + pseudotime:exon,
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
