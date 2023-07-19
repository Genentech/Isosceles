#' Check splicing status of transcripts detected from the BAM file(s)
#' @noRd
check_splicing_status <- function(tx_df,
                                  anno_data,
                                  bam_parsed,
                                  nr_intron_positions,
                                  tx_intron_idx,
                                  tx_anno_intron_idx,
                                  known_intron_granges = NULL,
                                  min_bam_splice_read_count = 2,
                                  min_bam_splice_fraction = 0.1) {

    # Check arguments
    assertthat::assert_that(is.data.frame(tx_df))
    assertthat::assert_that(is.list(anno_data))
    assertthat::assert_that(assertthat::has_name(anno_data, "intron_df"))
    assertthat::assert_that(is.data.frame(anno_data$intron_df))
    assertthat::assert_that(is.data.frame(bam_parsed))
    assertthat::assert_that(identical(
        colnames(bam_parsed), c("intron_positions", "read_count", "chromosome",
                                "start_positions","end_positions")
    ))
    assertthat::assert_that(is.character(nr_intron_positions))
    assertthat::assert_that(is.list(tx_intron_idx))
    assertthat::assert_that(is.list(tx_anno_intron_idx))
    if (!is.null(known_intron_granges)) {
        assertthat::assert_that(methods::is(known_intron_granges, "GRanges"))
    }
    assertthat::assert_that(assertthat::is.count(min_bam_splice_read_count))
    assertthat::assert_that(is.numeric(min_bam_splice_fraction))
    assertthat::assert_that(length(min_bam_splice_fraction) == 1)
    assertthat::assert_that(min_bam_splice_fraction >= 0)
    assertthat::assert_that(min_bam_splice_fraction <= 1)

    # Prepare a non-redundant known intron set
    nr_known_intron_positions <- unique(c(anno_data$intron_df$position,
                                          as.character(known_intron_granges)))
    nr_known_splice_5p_positions <-
        unique(gsub("-\\d+", "", nr_known_intron_positions))
    nr_known_splice_3p_positions <-
        unique(gsub("\\d+-", "", nr_known_intron_positions))

    # Prepare splice sites confirmed by aligned reads
    bam_splice_df <- data.frame(read_count = bam_parsed$read_count)
    bam_splice_df$intron_idx <- tx_intron_idx
    bam_splice_df <- tidyr::unchop(bam_splice_df, "intron_idx")
    bam_splice_df$splice_5p_position <-
        gsub("-\\d+", "", nr_intron_positions[bam_splice_df$intron_idx])
    bam_splice_df$splice_3p_position <-
        gsub("\\d+-", "", nr_intron_positions[bam_splice_df$intron_idx])
    bam_splice_df$is_splice_5p_known <-
        bam_splice_df$splice_5p_position %in% nr_known_splice_5p_positions
    bam_splice_df$is_splice_3p_known <-
        bam_splice_df$splice_3p_position %in% nr_known_splice_3p_positions
    bam_splice_df <- bam_splice_df %>%
        dplyr::group_by(.data$splice_5p_position) %>%
        dplyr::mutate(
            splice_5p_read_count = sum(.data$read_count)
        ) %>%
        dplyr::ungroup()
    bam_splice_df <- bam_splice_df %>%
        dplyr::group_by(.data$splice_3p_position) %>%
        dplyr::mutate(
            splice_3p_read_count = sum(.data$read_count)
        ) %>%
        dplyr::ungroup()

    bam_splice_5p_df <- bam_splice_df %>%
        dplyr::filter(.data$is_splice_3p_known) %>%
        dplyr::group_by(.data$splice_5p_position, .data$splice_3p_position) %>%
        dplyr::summarise(
            fraction = sum(.data$read_count) /
                unique(.data$splice_3p_read_count)
        ) %>%
        dplyr::ungroup()
    bam_splice_5p_positions <- bam_splice_5p_df %>%
        dplyr::filter(.data$fraction >= min_bam_splice_fraction) %>%
        dplyr::pull(.data$splice_5p_position) %>%
        unique() %>%
        setdiff(nr_known_splice_5p_positions)
    bam_splice_5p_positions <- bam_splice_df %>%
        dplyr::filter(.data$splice_5p_position %in% bam_splice_5p_positions) %>%
        dplyr::group_by(.data$splice_5p_position) %>%
        dplyr::summarise(read_count = sum(.data$read_count)) %>%
        dplyr::filter(.data$read_count >= min_bam_splice_read_count) %>%
        dplyr::pull(.data$splice_5p_position)
    bam_splice_5p_positions <- c(bam_splice_5p_positions,
                                 nr_known_splice_5p_positions)

    bam_splice_3p_df <- bam_splice_df %>%
        dplyr::filter(.data$is_splice_5p_known) %>%
        dplyr::group_by(.data$splice_5p_position, .data$splice_3p_position) %>%
        dplyr::summarise(
            fraction = sum(.data$read_count) /
                unique(.data$splice_5p_read_count)
        ) %>%
        dplyr::ungroup()
    bam_splice_3p_positions <- bam_splice_3p_df %>%
        dplyr::filter(.data$fraction >= min_bam_splice_fraction) %>%
        dplyr::pull(.data$splice_3p_position) %>%
        unique() %>%
        setdiff(nr_known_splice_3p_positions)
    bam_splice_3p_positions <- bam_splice_df %>%
        dplyr::filter(.data$splice_3p_position %in% bam_splice_3p_positions) %>%
        dplyr::group_by(.data$splice_3p_position) %>%
        dplyr::summarise(read_count = sum(.data$read_count)) %>%
        dplyr::filter(.data$read_count >= min_bam_splice_read_count) %>%
        dplyr::pull(.data$splice_3p_position)
    bam_splice_3p_positions <- c(bam_splice_3p_positions,
                                 nr_known_splice_3p_positions)

    # Calculate splicing support level
    tx_strand <- methods::as(tx_df$position, "GRanges") %>%
        BiocGenerics::strand() %>%
        as.character()
    has_correct_strand <- tx_strand != "*"
    has_annotated_gene <- !is.na(tx_df$gene_id)
    is_fusion_transcript <- tx_df$compatible_gene_count > 1
    tx_anno_compatibility_df <- check_splicing_compatibility(tx_intron_idx,
                                                             tx_anno_intron_idx)
    is_anno_compatible <- rep(FALSE, nrow(tx_df))
    is_anno_compatible[unique(tx_anno_compatibility_df$subject_idx)] <- TRUE
    are_all_intron_known <- sapply(tx_intron_idx, function(x) {
        all(fastmatch::`%fin%`(nr_intron_positions[x], nr_known_intron_positions))
    })
    nr_splice_5p_positions <- gsub("-\\d+", "", nr_intron_positions)
    nr_splice_3p_positions <- gsub("\\d+-", "", nr_intron_positions)
    are_all_splice_known <- sapply(tx_intron_idx, function(x) {
        all(fastmatch::`%fin%`(nr_splice_5p_positions[x],
                               nr_known_splice_5p_positions)) &
            all(fastmatch::`%fin%`(nr_splice_3p_positions[x],
                                   nr_known_splice_3p_positions))
    })
    are_all_splice_bam <- sapply(tx_intron_idx, function(x) {
        all(fastmatch::`%fin%`(nr_splice_5p_positions[x],
                               bam_splice_5p_positions)) &
            all(fastmatch::`%fin%`(nr_splice_3p_positions[x],
                                   bam_splice_3p_positions))
    })
    splicing_support_level <-
        ifelse(!has_correct_strand, "AX",
               ifelse(is_fusion_transcript, "AF",
                      ifelse(is_anno_compatible, "PC",
                             ifelse(are_all_intron_known, "EC",
                                    ifelse(are_all_splice_known, "NC",
                                           ifelse(are_all_splice_bam, "DN",
                                                  ifelse(has_annotated_gene, "AS",
                                                         "AX")))))))

    return(splicing_support_level)
}
