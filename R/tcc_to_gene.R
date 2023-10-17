#' Prepare a gene-level SummarizedExperiment object
#'
#' Prepares a gene-level SummarizedExperiment from TCC data.
#'
#' @param se_tcc A TCC SummarizedExperiment object returned by a function from
#' the \code{\link{Isosceles-package}}.
#' @return A SummarizedExperiment object containing gene annotation and
#' quantification data.
#' @export
tcc_to_gene <- function(se_tcc) {

    # Check arguments
    assertthat::assert_that(methods::is(se_tcc, "SummarizedExperiment"))
    assertthat::assert_that(is.element(
        "counts", SummarizedExperiment::assayNames(se_tcc)
    ))
    assertthat::assert_that(is.element(
        "tpm", SummarizedExperiment::assayNames(se_tcc)
    ))
    assertthat::assert_that(is.element(
        "relative_expression", SummarizedExperiment::assayNames(se_tcc)
    ))
    assertthat::assert_that(is.element(
        "gene_id", colnames(SummarizedExperiment::rowData(se_tcc))
    ))
    assertthat::assert_that(is.element(
        "gene_name", colnames(SummarizedExperiment::rowData(se_tcc))
    ))

    # Prepare input data
    is_bulk_rnaseq <- is.matrix(SummarizedExperiment::assay(se_tcc, "counts"))
    gene_ids <- SummarizedExperiment::rowData(se_tcc)$gene_id
    transcript_df <- S4Vectors::metadata(se_tcc)$transcript_df
    gene_df <- unique(transcript_df[, c("gene_id", "gene_name")])
    rownames(gene_df) <- NULL

    # Prepare missing gene matrix
    missing_gene_ids <- setdiff(gene_df$gene_id, unique(gene_ids))
    missing_gene_matrix <- matrix(
        0, nrow = length(missing_gene_ids), ncol = ncol(se_tcc)
    )
    rownames(missing_gene_matrix) <- missing_gene_ids

    # Prepare gene-level assays
    gene_counts <- rowsum(SummarizedExperiment::assay(se_tcc, "counts"),
                          gene_ids)
    gene_counts <- rbind(gene_counts, missing_gene_matrix)
    gene_counts <- gene_counts[gene_df$gene_id,,drop = FALSE]
    if (!is_bulk_rnaseq) {
        gene_counts <- methods::as(gene_counts, "dgCMatrix")
    }
    gene_tpm <- rowsum(SummarizedExperiment::assay(se_tcc, "tpm"),
                       gene_ids)
    gene_tpm <- rbind(gene_tpm, missing_gene_matrix)
    gene_tpm <- gene_tpm[gene_df$gene_id,,drop = FALSE]
    if (!is_bulk_rnaseq) {
        gene_tpm <- methods::as(gene_tpm, "dgCMatrix")
    }
    gene_rel_expr <- rowsum(
        SummarizedExperiment::assay(se_tcc, "relative_expression"), gene_ids
    )
    gene_rel_expr <- rbind(gene_rel_expr, missing_gene_matrix)
    gene_rel_expr <- gene_rel_expr[gene_df$gene_id,,drop = FALSE]
    if (!is_bulk_rnaseq) {
        gene_rel_expr <- methods::as(gene_rel_expr, "dgCMatrix")
    }

    # Prepare the SummarizedExperiment object
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            counts = gene_counts,
            tpm = gene_tpm,
            relative_expression = gene_rel_expr
        ),
        rowData = gene_df,
        colData = SummarizedExperiment::colData(se_tcc)
    )

    return(se)
}
