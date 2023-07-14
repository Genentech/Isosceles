#' Calculate TPM values
#'
#' Calculates TPM values from counts stored in a SummarizedExperiment
#' or SingleCellExperiment object
#'
#' @param se A SummarizedExperiment or SingleCellExperiment object containing
#' a 'counts' assay
#' @return A matrix (regular or sparse, depending on the input data) containing
#' the calculated TPM values
#' @keywords internal
calculate_tpm <- function(se) {

    # Check arguments
    assertthat::assert_that(
        grepl("SummarizedExperiment", class(se)) ||
            grepl("SingleCellExperiment", class(se))
    )
    assertthat::assert_that(is.element(
        "counts", SummarizedExperiment::assayNames(se)
    ))

    # Prepare input data
    is_bulk_rnaseq <- is.matrix(SummarizedExperiment::assay(se, "counts"))

    # Convert the count matrix to a triplet representation data frame
    count_matrix <- SummarizedExperiment::assay(se, "counts")
    if (is_bulk_rnaseq) {
        tpm_df <- data.frame(
            i = rep(seq(nrow(count_matrix)), times = ncol(count_matrix)),
            j = rep(seq(ncol(count_matrix)), each = nrow(count_matrix)),
            x = as.numeric(count_matrix)
        )
    } else {
        tpm_df <- as.data.frame(Matrix::summary(count_matrix))
    }

    # Calculate TPM values
    tpm_df <- tpm_df %>%
        dplyr::group_by(.data$j) %>%
        dplyr::mutate(x = .data$x / sum(.data$x) * 1e6) %>%
        dplyr::ungroup()
    tpm_df$x[is.nan(tpm_df$x)] <- 0

    # Convert the data frame back to a matrix
    if (is_bulk_rnaseq) {
        tpm_matrix <- matrix(tpm_df$x,
                             nrow = nrow(count_matrix),
                             ncol = ncol(count_matrix))
    } else {
        tpm_matrix <- Matrix::sparseMatrix(
            i = tpm_df$i,
            j = tpm_df$j,
            x = tpm_df$x,
            dims = c(nrow(count_matrix), ncol(count_matrix))
        )
        tpm_matrix <- Matrix::drop0(tpm_matrix)
    }
    rownames(tpm_matrix) <- rownames(count_matrix)
    colnames(tpm_matrix) <- colnames(count_matrix)

    return(tpm_matrix)
}
