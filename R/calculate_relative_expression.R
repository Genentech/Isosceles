#' Calculate relative expression
#'
#' Calculates relative expression from TPM values stored in a SummarizedExperiment
#' or SingleCellExperiment object
#'
#' @param se A SummarizedExperiment or SingleCellExperiment object containing
#' a 'tpm' assay and a 'gene_id' column in its row metadata
#' @return A matrix (regular or sparse, depending on the input data) containing
#' the calculated relative expression values
#' @keywords internal
calculate_relative_expression <- function(se) {

    # Check arguments
    assertthat::assert_that(methods::is(se, "SummarizedExperiment"))
    assertthat::assert_that(is.element(
        "tpm", SummarizedExperiment::assayNames(se)
    ))
    assertthat::assert_that(is.element(
        "gene_id", colnames(SummarizedExperiment::rowData(se))
    ))

    # Prepare input data
    is_bulk_rnaseq <- is.matrix(SummarizedExperiment::assay(se, "tpm"))
    gene_ids <- SummarizedExperiment::rowData(se)$gene_id
    gene_ids <- as.numeric(as.factor(gene_ids))

    # Convert the TPM matrix to a triplet representation data frame
    tpm_matrix <- SummarizedExperiment::assay(se, "tpm")
    if (is_bulk_rnaseq) {
        relative_expression_df <- data.frame(
            i = rep(seq(nrow(tpm_matrix)), times = ncol(tpm_matrix)),
            j = rep(seq(ncol(tpm_matrix)), each = nrow(tpm_matrix)),
            x = as.numeric(tpm_matrix)
        )
    } else {
        relative_expression_df <- as.data.frame(Matrix::summary(tpm_matrix))
    }

    # Calculate relative expression values
    relative_expression_df$gene_id <- gene_ids[relative_expression_df$i]
    relative_expression_df <- relative_expression_df %>%
        dplyr::group_by(.data$gene_id, .data$j) %>%
        dplyr::mutate(x = .data$x / sum(.data$x)) %>%
        dplyr::ungroup()
    relative_expression_df$x[is.nan(relative_expression_df$x)] <- 0

    # Convert the data frame back to a matrix
    if (is_bulk_rnaseq) {
        relative_expression_matrix <- matrix(relative_expression_df$x,
                                             nrow = nrow(tpm_matrix),
                                             ncol = ncol(tpm_matrix))
    } else {
        relative_expression_matrix <- Matrix::sparseMatrix(
            i = relative_expression_df$i,
            j = relative_expression_df$j,
            x = relative_expression_df$x,
            dims = c(nrow(tpm_matrix), ncol(tpm_matrix))
        )
        relative_expression_matrix <- Matrix::drop0(relative_expression_matrix)
    }
    rownames(relative_expression_matrix) <- rownames(tpm_matrix)
    colnames(relative_expression_matrix) <- colnames(tpm_matrix)

    return(relative_expression_matrix)
}
