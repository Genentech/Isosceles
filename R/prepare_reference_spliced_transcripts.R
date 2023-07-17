#' Prepare reference spliced transcript data
#'
#' Prepares reference spliced transcript data for further processing
#'
#' @param anno_data A list containing genome annotation data returned by
#' the \code{\link{prepare_reference_annotations}} function
#' @param bin_size An integer scalar specifying the bin size for transcript
#' start and end position binning
#' @return A named list containing following elements:
#' \describe{
#'   \item{tx_df}{a data frame storing reference spliced transcript data}
#'   \item{tx_granges}{a GRanges object storing genomic positions of reference spliced transcript}
#'   \item{tx_exon_granges_list}{a GRangesList object storing exon genomic positions of reference spliced transcript}
#'   \item{tx_intron_granges_list}{a GRangesList object storing intron genomic positions of reference spliced transcript}
#' }
#' @keywords internal
prepare_reference_spliced_transcripts <- function(anno_data,
                                                  bin_size = 50) {

    # Check arguments
    assertthat::assert_that(is.list(anno_data))
    assertthat::assert_that(assertthat::has_name(anno_data, "transcript_df"))
    assertthat::assert_that(is.data.frame(anno_data$transcript_df))
    assertthat::assert_that(assertthat::is.count(bin_size))

    # Helper functions
    make_hash <- function(x) {
        x %>%
            sapply(digest::digest) %>%
            substring(0, 16) %>%
            `names<-`(NULL)
    }

    # Prepare spliced transcript data
    transcript_df <- anno_data$transcript_df %>%
        dplyr::filter(.data$is_spliced)

    # Prepare transcript position GRanges object
    tx_granges <- methods::as(transcript_df$position, "GRanges")
    tx_chromosome <- as.character(GenomeInfoDb::seqnames(tx_granges))
    tx_strand <- as.character(BiocGenerics::strand(tx_granges))
    tx_start <- BiocGenerics::start(tx_granges)
    tx_end <- BiocGenerics::end(tx_granges)
    tx_start_bin <- floor(tx_start / bin_size) * bin_size
    tx_end_bin <- (floor(tx_end / bin_size) * bin_size) + bin_size

    # Prepare basic transcript information
    tx_df <- transcript_df %>%
        dplyr::transmute(
            hash_id = make_hash(.data$intron_positions),
            tx_start_bin = tx_start_bin,
            tx_end_bin = tx_end_bin,
            tx_chromosome = tx_chromosome,
            tx_start = tx_start,
            tx_end = tx_end,
            tx_strand = tx_strand,
            intron_positions = .data$intron_positions,
            gene_id = .data$gene_id,
            gene_name = .data$gene_name,
            transcript_id = .data$transcript_id
        ) %>%
        dplyr::group_by(.data$hash_id, .data$tx_start_bin, .data$tx_end_bin) %>%
        dplyr::summarise(
            chromosome = unique(.data$tx_chromosome),
            start = min(.data$tx_start),
            end = max(.data$tx_end),
            strand = unique(.data$tx_strand),
            intron_positions = unique(.data$intron_positions),
            gene_id = unique(.data$gene_id),
            gene_name = unique(.data$gene_name),
            compatible_tx = paste0(unique(.data$transcript_id), collapse = ","),
            .groups = "drop"
        )  %>%
        dplyr::transmute(
            hash_id = .data$hash_id,
            position = paste0(.data$chromosome, ":", .data$start, "-", .data$end,
                              ":", .data$strand),
            intron_positions = .data$intron_positions,
            gene_id = .data$gene_id,
            gene_name = .data$gene_name,
            compatible_gene_ids = .data$gene_id,
            compatible_gene_names = .data$gene_name,
            compatible_tx = .data$compatible_tx,
            splicing_support_level = "AP",
            fivethree_support_level = "FL",
            read_count = as.integer(NA),
            relative_expression = as.numeric(NA)
        ) %>%
        as.data.frame()

    # Prepare transcript position GRanges object
    tx_granges <- methods::as(tx_df$position, "GRanges")
    names(tx_granges) <- tx_df$hash_id

    # Prepare transcript intron & exon GRangesList object
    tx_intron_granges_df <- tx_df$intron_positions %>%
        strsplit(",") %>%
        tibble::enframe() %>%
        dplyr::rename(tx_idx = "name", intron_position = "value") %>%
        tidyr::unchop("intron_position")
    tx_intron_granges <- methods::as(
        tx_intron_granges_df$intron_position, "GRanges"
    )
    tx_intron_granges_list <- GenomicRanges::split(
        tx_intron_granges, tx_intron_granges_df$tx_idx
    )
    names(tx_intron_granges_list) <- tx_df$hash_id
    tx_exon_granges_list <- IRanges::psetdiff(tx_granges, tx_intron_granges_list)

    return(list(
        tx_df = tx_df,
        tx_granges = tx_granges,
        tx_exon_granges_list = tx_exon_granges_list,
        tx_intron_granges_list = tx_intron_granges_list
    ))
}
