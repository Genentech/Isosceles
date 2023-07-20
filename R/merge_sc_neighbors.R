#' Merging neighboring cell TCC values in scRNA-Seq data
#'
#' Prepares a TCC SummarizedExperiment object where count values from the
#' nearest k neighbors are added to the count values of each cell
#'
#' @param se_tcc A TCC SummarizedExperiment object returned by the
#' \code{\link{prepare_tcc_se}} function
#' @param pca_mat A matrix containing PCA coordinates of each cell
#' @param k An integer scalar specifying the number of nearest neighbors to use
#' @param use_annoy A logical scalar indicating whether to use the Annoy
#' algorithm for approximate nearest neighbor identification (recommended for
#' big datasets)
#' @param ncpu An integer scalar specifying the number of cores to use for
#' multicore parallelization
#' @return A SummarizedExperiment object containing merged TCC data
#' @export
merge_sc_neighbors <- function(se_tcc,
                               pca_mat,
                               k = 10,
                               use_annoy = FALSE,
                               ncpu = 1) {

    # Check arguments
    assertthat::assert_that(methods::is(se_tcc, "SummarizedExperiment"))
    assertthat::assert_that(is.element(
        "counts", SummarizedExperiment::assayNames(se_tcc)
    ))
    assertthat::assert_that(is.matrix(pca_mat))
    assertthat::assert_that(identical(colnames(se_tcc), rownames(pca_mat)))
    assertthat::assert_that(assertthat::is.count(k))
    assertthat::assert_that(assertthat::is.flag(use_annoy))
    assertthat::assert_that(assertthat::is.count(ncpu))

    # Prepare the parallelization backend
    BPPARAM <- BiocParallel::MulticoreParam(ncpu)

    # Prepare the nearest neighbor search backend
    BNPARAM <- BiocNeighbors::KmknnParam()
    if (use_annoy) {
        BNPARAM <- BiocNeighbors::AnnoyParam()
    }

    # Calculate the k nearest neighbors
    knn_results <- BiocNeighbors::findKNN(pca_mat, k, get.distance = FALSE,
                                          BNPARAM = BNPARAM, BPPARAM = BPPARAM)
    knn_results <- knn_results$index
    knn_results <- cbind(seq(ncol(se_tcc)), knn_results)

    # Merge TCC values from the neighboring cells
    count_matrix <- SummarizedExperiment::assay(se_tcc, "counts")
    merged_count_matrix <- scuttle::sumCountsAcrossFeatures(
        Matrix::t(count_matrix),
        ids = as.list(as.data.frame(t(knn_results))),
        average = FALSE
    )
    merged_count_matrix <- methods::as(merged_count_matrix, "dgCMatrix")
    merged_count_matrix <- Matrix::t(merged_count_matrix)
    colnames(merged_count_matrix) <- colnames(count_matrix)

    # Prepare the SummarizedExperiment object
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = merged_count_matrix),
        rowData = SummarizedExperiment::rowData(se_tcc),
        colData = SummarizedExperiment::colData(se_tcc)
    )
    S4Vectors::metadata(se) <- S4Vectors::metadata(se_tcc)
    SummarizedExperiment::assay(se, "tpm") <- calculate_tpm(se)
    SummarizedExperiment::assay(se, "relative_expression") <-
        calculate_relative_expression(se)

    return(se)
}
