#' Plot a PSI heatmap
#'
#' Creates a heatmap of PSI (Percent Spliced In) values for the regions of a
#' given gene across samples or cells
#'
#' @param se_psi A PSI SummarizedExperiment object returned by the
#' \code{\link{prepare_psi_se}} function
#' @param gene_id A string containing the identifier of the gene to plot
#' @param heatmap_colors A character vector containing the color palette used
#' in the heatmap
#' @param region_colors A named character vector of colors for the region type
#' annotations
#' @param \dots Additional parameters for the plot, passed to the
#' \code{\link{pheatmap}} function
#' @return A plot object
#' @export
plot_psi_heatmap <- function(se_psi,
                             gene_id,
                             heatmap_colors = viridis::cividis(100),
                             region_colors = NULL,
                             ...) {

    # Check arguments
    assertthat::assert_that(methods::is(se_psi, "SummarizedExperiment"))
    assertthat::assert_that(is.element(
        "psi", SummarizedExperiment::assayNames(se_psi)
    ))
    assertthat::assert_that(is.element(
        "gene_id", colnames(SummarizedExperiment::rowData(se_psi))
    ))
    assertthat::assert_that(assertthat::is.string(gene_id))
    assertthat::assert_that(is.element(
        gene_id, SummarizedExperiment::rowData(se_psi)$gene_id
    ))
    assertthat::assert_that(is.character(heatmap_colors))
    if (!is.null(region_colors)) {
        assertthat::assert_that(is.character(region_colors))
        assertthat::assert_that(identical(
            sort(names(region_colors)),
            c("A3", "A5", "CE", "RI", "TES", "TSS")
        ))
    }

    # Filter the PSI SummarizedExperiment object
    se <- se_psi[
        SummarizedExperiment::rowData(se_psi)$gene_id == gene_id,
    ]

    # Prepare a PSI values matrix and remove cells with no gene expression
    psi_matrix <- se %>%
        SummarizedExperiment::assay("psi") %>%
        Matrix::t() %>%
        as.matrix()
    psi_matrix <- psi_matrix[rowSums(psi_matrix) > 0,]

    # Prepare the heatmap column annotations data frame
    col_anno_df <- SummarizedExperiment::rowData(se) %>%
        as.data.frame %>%
        dplyr::transmute(region = factor(
            .data$type, levels = c("TSS", "CE", "RI", "A5", "A3", "TES")
        ))

    # Prepare region colors
    annotation_colors <- NA
    if (!is.null(region_colors)) {
        annotation_colors <- list(region = region_colors)
    }

    # Create the PSI value heatmap
    psi_heatmap <- pheatmap::pheatmap(
        mat = psi_matrix,
        annotation_col = col_anno_df,
        color = heatmap_colors,
        cluster_rows = TRUE, cluster_cols = FALSE,
        show_rownames = FALSE, show_colnames = FALSE,
        treeheight_row = 0, annotation_names_col = FALSE,
        annotation_colors = annotation_colors,
        ...
    )

    return(psi_heatmap)
}
