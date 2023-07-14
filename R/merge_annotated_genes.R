#' Merging genes by shared introns
#'
#' Merge IDs and names of genes sharing annotated introns
#'
#' @param gene_df A data frame containing gene IDs and gene names
#' @param intron_df A data frame containing gene IDs and intron positions
#' @return A copy of the gene_df data frame with merged gene IDs and gene names
#' added to it
#' @keywords internal
merge_annotated_genes <- function(gene_df, intron_df) {

    # Check arguments
    assertthat::assert_that(is.data.frame(gene_df))
    assertthat::assert_that(is.data.frame(intron_df))
    assertthat::assert_that(assertthat::has_name(gene_df, "gene_id"))
    assertthat::assert_that(assertthat::has_name(gene_df, "gene_name"))
    assertthat::assert_that(assertthat::has_name(intron_df, "gene_id"))
    assertthat::assert_that(assertthat::has_name(intron_df, "position"))

    # Indexes of genes sharing an intron
    intron_pairs <- data.frame(
        intron_idx = as.numeric(as.factor(intron_df$position)),
        gene_idx = match(intron_df$gene_id, gene_df$gene_id)
    )
    intron_pairs <- dplyr::distinct(intron_pairs)
    intron_pairs <- intron_pairs %>%
        dplyr::full_join(intron_pairs, by = "intron_idx") %>%
        dplyr::select(-.data$intron_idx) %>%
        dplyr::filter(.data$gene_idx.x < .data$gene_idx.y) %>%
        dplyr::distinct()

    # Creating a shared intron graph & gene clustering
    gene_graph <- igraph::graph.data.frame(
        intron_pairs, directed = FALSE,
        vertices = seq_along(gene_df$gene_id)
    )
    gene_cluster_id <- igraph::clusters(gene_graph)$membership

    # Preparing merged gene data
    merged_gene_df <- gene_df %>%
        dplyr::mutate(gene_cluster_id = gene_cluster_id) %>%
        dplyr::group_by(.data$gene_cluster_id) %>%
        dplyr::mutate(
            merged_gene_id = paste0(.data$gene_id, collapse = ","),
            merged_gene_name = paste0(.data$gene_name, collapse = ",")
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(-.data$gene_cluster_id) %>%
        as.data.frame()

    return(merged_gene_df)
}
