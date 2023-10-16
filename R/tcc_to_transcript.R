#' Prepare a transcript-level SummarizedExperiment object
#'
#' Prepares a transcript-level SummarizedExperiment from TCC data using the EM
#' algorithm.
#'
#' @param se_tcc A TCC SummarizedExperiment object returned by a function from
#' the \code{\link{Isosceles-package}}.
#' @param em.maxiter An integer scalar specifying the maximum number of EM
#' iterations.
#' @param em.conv A numeric scalar specifying the EM convergence threshold.
#' @param use_length_normalization A logical scalar specifying if normalization
#' using effective transcript lengths should be used during EM.
#' @param ncpu An integer scalar specifying the number of cores to use for
#' multicore parallelization.
#' @return A SummarizedExperiment object containing transcript annotation and
#' quantification data.
#' @export
tcc_to_transcript <- function(se_tcc,
                              em.maxiter = 250,
                              em.conv = 0.01,
                              use_length_normalization = TRUE,
                              ncpu = 1) {

    # Check arguments
    assertthat::assert_that(methods::is(se_tcc, "SummarizedExperiment"))
    assertthat::assert_that(is.element(
        "counts", SummarizedExperiment::assayNames(se_tcc)
    ))
    assertthat::assert_that(is.element(
        "gene_id", colnames(SummarizedExperiment::rowData(se_tcc))
    ))
    assertthat::assert_that(assertthat::is.count(em.maxiter))
    assertthat::assert_that(assertthat::is.number(em.conv))
    assertthat::assert_that(assertthat::is.count(ncpu))
    assertthat::assert_that(assertthat::is.flag(use_length_normalization))

    # Prepare the parallelization backend
    BPPARAM <- BiocParallel::MulticoreParam(ncpu)

    # Prepare input data
    ec_count_matrix <- SummarizedExperiment::assay(se_tcc, "counts")
    is_bulk_rnaseq <- is.matrix(ec_count_matrix)
    sample_ids <- BiocGenerics::colnames(se_tcc)
    gene_ids <- unique(S4Vectors::metadata(se_tcc)$transcript_df$gene_id)
    mean_read_length <- S4Vectors::metadata(se_tcc)$mean_read_length
    tx_exon_granges_list <- S4Vectors::metadata(se_tcc)$transcript_exon_granges_list
    tx_effective_lengths <- BiocGenerics::width(tx_exon_granges_list) %>%
        BiocGenerics::sapply(sum) %>%
        unname() %>%
        as.numeric()
    tx_effective_lengths <- ifelse(tx_effective_lengths < mean_read_length,
                                   mean_read_length, tx_effective_lengths)

    # Run EM for each sample
    tx_count_matrix <- BiocParallel::bplapply(sample_ids, function(sample_id) {
        ec_counts_sample <- unname(ec_count_matrix[, sample_id])
        ec_selector_em_genes <- ec_counts_sample > 0
        em_gene_ids <- unique(
            SummarizedExperiment::rowData(se_tcc)$gene_id[ec_selector_em_genes]
        )
        tx_selector_em_genes <-
            S4Vectors::metadata(se_tcc)$transcript_df$gene_id %in% em_gene_ids
        tx_counts_em_genes <- lapply(em_gene_ids, function(gene_id) {
            ec_selector <- SummarizedExperiment::rowData(se_tcc)$gene_id == gene_id
            tx_selector <- S4Vectors::metadata(se_tcc)$transcript_df$gene_id == gene_id
            ec_counts_gene <- ec_counts_sample[ec_selector]
            if (length(ec_counts_gene) == 0) {
                return(rep(0, sum(tx_selector)))
            }
            if (length(ec_counts_gene) == 1) {
                ec_counts_gene <- as.matrix(ec_counts_gene)
            } else {
                ec_counts_gene <- diag(ec_counts_gene)
            }
            compatibility_matrix <- as.matrix(
                S4Vectors::metadata(se_tcc)$compatibility_matrix[ec_selector, tx_selector, drop = FALSE]
            )
            tx_effective_lengths_gene <- tx_effective_lengths[tx_selector]
            tx_counts_gene <- EM(
                counts = ec_counts_gene,
                compatibility_matrix = compatibility_matrix,
                use_length_normalization = use_length_normalization,
                tx_effective_lengths = tx_effective_lengths_gene,
                maxiter = em.maxiter,
                conv = em.conv
            )
            return(tx_counts_gene)
        })
        tx_counts_em_genes <- unlist(tx_counts_em_genes)
        tx_counts_sample <- rep(
            0, nrow(S4Vectors::metadata(se_tcc)$transcript_df)
        )
        tx_counts_sample[tx_selector_em_genes] <- tx_counts_em_genes
        tx_counts_sample <- as.matrix(tx_counts_sample)
        if (!is_bulk_rnaseq) {
            tx_counts_sample <- methods::as(tx_counts_sample, "dgCMatrix")
        }
        return(tx_counts_sample)
    }, BPPARAM = BPPARAM)

    # Prepare the transcript count matrix
    tx_count_matrix <- do.call(cbind, tx_count_matrix)
    colnames(tx_count_matrix) <- sample_ids
    rownames(tx_count_matrix) <-
        S4Vectors::metadata(se_tcc)$transcript_df$transcript_id

    # Prepare the SummarizedExperiment object
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = tx_count_matrix),
        rowRanges = S4Vectors::metadata(se_tcc)$transcript_exon_granges_list,
        colData = SummarizedExperiment::colData(se_tcc)
    )
    SummarizedExperiment::rowData(se) <- S4Vectors::metadata(se_tcc)$transcript_df
    SummarizedExperiment::assay(se, "tpm") <- calculate_tpm(se)
    SummarizedExperiment::assay(se, "relative_expression") <-
        calculate_relative_expression(se)

    return(se)
}
