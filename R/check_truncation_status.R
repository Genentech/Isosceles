#' Check truncation status of transcripts detected from the BAM file(s)
#' @noRd
check_truncation_status <- function(tx_df,
                                    tx_exon_granges_list,
                                    tx_intron_granges_list,
                                    anno_data) {

    # Check arguments
    assertthat::assert_that(is.data.frame(tx_df))
    assertthat::assert_that(assertthat::has_name(tx_df, "hash_id"))
    assertthat::assert_that(grepl("GRangesList", class(tx_exon_granges_list)))
    assertthat::assert_that(grepl("GRangesList", class(tx_intron_granges_list)))
    assertthat::assert_that(is.list(anno_data))
    assertthat::assert_that(
        assertthat::has_name(anno_data, "transcript_first_last_df")
    )

    # Prepare transcript data
    tx_trunc_df <- tx_df %>%
        dplyr::filter(.data$splicing_support_level %in% c("PC", "EC", "NC", "DN")) %>%
        dplyr::filter(!is.na(.data$gene_id)) %>%
        dplyr::select("hash_id", "gene_id")
    tx_splice_first_last_df <- get_first_last_grange(
        tx_intron_granges_list[tx_trunc_df$hash_id]
    ) %>%
        dplyr::transmute(
            hash_id = .data$feature_id,
            strand_tx = as.character(BiocGenerics::strand(
                methods::as(.data$first_grange, "GRanges")
            )),
            first_splice_5p_tx = ifelse(.data$strand_tx == "+",
                                        gsub("-\\d+", "", .data$first_grange),
                                        gsub("\\d+-", "", .data$first_grange)),
            last_splice_3p_tx = ifelse(.data$strand_tx == "+",
                                       gsub("\\d+-", "", .data$last_grange),
                                       gsub("-\\d+", "", .data$last_grange))
        )
    tx_trunc_df <- tx_trunc_df %>%
        dplyr::left_join(tx_splice_first_last_df)

    # Prepare reference transcript data
    anno_splice_first_last_df <- anno_data$transcript_first_last_df %>%
        dplyr::transmute(
            gene_id = .data$gene_id,
            strand_ref = as.character(BiocGenerics::strand(
                methods::as(.data$first_intron_ref, "GRanges")
            )),
            first_splice_5p_ref = ifelse(.data$strand_ref == "+",
                                         gsub("-\\d+", "", .data$first_intron_ref),
                                         gsub("\\d+-", "", .data$first_intron_ref)),
            last_splice_3p_ref = ifelse(.data$strand_ref == "+",
                                        gsub("\\d+-", "", .data$last_intron_ref),
                                        gsub("-\\d+", "", .data$last_intron_ref))
        )
    tx_trunc_df <- tx_trunc_df %>%
        dplyr::inner_join(anno_splice_first_last_df)


    # Calculate 5' and 3' support
    tx_trunc_df$has_5p_support <-
        tx_trunc_df$first_splice_5p_tx == tx_trunc_df$first_splice_5p_ref
    tx_trunc_df$has_3p_support <-
        tx_trunc_df$last_splice_3p_tx == tx_trunc_df$last_splice_3p_ref

    # Calculate fivethree support level
    tx_trunc_df <- tx_trunc_df %>%
        dplyr::group_by(.data$hash_id) %>%
        dplyr::summarise(
            has_5p_support = any(.data$has_5p_support),
            has_3p_support = any(.data$has_3p_support)
        ) %>%
        dplyr::transmute(
            hash_id = .data$hash_id,
            fivethree_support_level = ifelse(.data$has_5p_support &
                                                 .data$has_3p_support, "FL",
                                             ifelse(.data$has_5p_support, "3T",
                                                    ifelse(.data$has_3p_support, "5T",
                                                           "FT")))
        )
    tx_fivethree_support_level <- stats::setNames(
        tx_trunc_df$fivethree_support_level, tx_trunc_df$hash_id
    )
    tx_fivethree_support_level <- tx_fivethree_support_level[tx_df$hash_id]
    tx_fivethree_support_level <- unname(tx_fivethree_support_level)
    tx_fivethree_support_level[is.na(tx_fivethree_support_level)] <- "NA"

    return(tx_fivethree_support_level)
}
