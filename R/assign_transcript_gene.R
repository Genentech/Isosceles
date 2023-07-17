#' Assign transcripts detected from the BAM file(s) to genes
#' @noRd
assign_transcript_gene <- function(tx_df,
                                   anno_data) {

    # Check arguments
    assertthat::assert_that(is.data.frame(tx_df))
    assertthat::assert_that(is.list(anno_data))
    assertthat::assert_that(assertthat::has_name(anno_data, "gene_df"))
    assertthat::assert_that(is.data.frame(anno_data$gene_df))
    assertthat::assert_that(assertthat::has_name(anno_data, "intron_df"))
    assertthat::assert_that(is.data.frame(anno_data$intron_df))

    # Prepare annotated gene & splice data data
    gene_id_to_gene_name <- stats::setNames(
        anno_data$gene_df$gene_name, anno_data$gene_df$gene_id
    )
    gene_anno_intron_df <- anno_data$intron_df %>%
        dplyr::select("gene_id", "position") %>%
        dplyr::distinct()
    gene_anno_splice_5p_df <- gene_anno_intron_df %>%
        dplyr::mutate(position = gsub("-\\d+", "", .data$position))
    gene_anno_splice_3p_df <- gene_anno_intron_df %>%
        dplyr::mutate(position = gsub("\\d+-", "", .data$position))

    # Get transcript splice sites
    tx_gene_intron_df <- tx_df %>%
        dplyr::select("hash_id", "intron_positions") %>%
        dplyr::transmute(
            hash_id = .data$hash_id,
            position = strsplit(.data$intron_positions, ",")
        ) %>%
        tidyr::unchop("position")
    tx_gene_splice_5p_df <- tx_gene_intron_df %>%
        dplyr::mutate(position = gsub("-\\d+", "", .data$position))
    tx_gene_splice_3p_df <- tx_gene_intron_df %>%
        dplyr::mutate(position = gsub("\\d+-", "", .data$position))

    # Detect transcripts sharing splice sites with annotated genes
    tx_gene_splice_5p_df <- tx_gene_splice_5p_df %>%
        dplyr::inner_join(gene_anno_splice_5p_df,
                          relationship = "many-to-many")
    tx_gene_splice_3p_df <- tx_gene_splice_3p_df %>%
        dplyr::inner_join(gene_anno_splice_3p_df,
                          relationship = "many-to-many")
    tx_gene_splice_df <- rbind(tx_gene_splice_5p_df, tx_gene_splice_3p_df)
    tx_gene_splice_df <- tx_gene_splice_df %>%
        dplyr::group_by(.data$hash_id) %>%
        dplyr::summarise(
            compatible_gene_count = length(unique(.data$gene_id)),
            compatible_gene_ids = paste0(unique(.data$gene_id), collapse = "|"),
            compatible_gene_names = paste0(gene_id_to_gene_name[unique(.data$gene_id)],
                                           collapse = "|")
        )

    # Assign transcripts to annotated genes
    tx_gene_df <- tx_df %>%
        dplyr::select("hash_id") %>%
        dplyr::left_join(tx_gene_splice_df) %>%
        tidyr::replace_na(list(compatible_gene_count = 0))
    tx_gene_df$gene_id <- ifelse(
        tx_gene_df$compatible_gene_count == 1, tx_gene_df$compatible_gene_ids, NA
    )
    tx_gene_df$gene_name <- ifelse(
        tx_gene_df$compatible_gene_count == 1, tx_gene_df$compatible_gene_names, NA
    )
    tx_gene_df <- tx_gene_df[, c("hash_id", "gene_id", "gene_name",
                                 "compatible_gene_ids", "compatible_gene_names",
                                 "compatible_gene_count")]

    return(tx_gene_df)
}
