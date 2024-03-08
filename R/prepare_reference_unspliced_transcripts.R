#' Prepare reference unspliced transcript data
#'
#' Prepares reference unspliced transcript data for further processing.
#'
#' @param anno_data A list containing genome annotation data returned by
#' the \code{\link{prepare_reference_annotations}} function.
#' @param bin_size An integer scalar specifying the bin size for transcript
#' start and end position binning.
#' @param use_full_hash A logical scalar specifying if full value of the MD5
#' hash (32 characters) should be used for the stable hash identifier rather
#' than its 16-character substring. This option should not be used unless you
#' encounter a hashing collision error (extremely unlikely).
#' @return A named list containing following elements:
#' \describe{
#'   \item{tx_df}{a data frame storing reference unspliced transcript data}
#'   \item{tx_granges}{a GRanges object storing genomic positions of reference unspliced transcript}
#'   \item{tx_exon_granges_list}{a GRangesList object storing exon genomic positions of reference unspliced transcript}
#'   \item{tx_intron_granges_list}{a GRangesList object storing intron genomic positions of reference unspliced transcript}
#' }
#' @keywords internal
prepare_reference_unspliced_transcripts <- function(anno_data,
                                                    bin_size = 50,
                                                    use_full_hash = FALSE) {

    # Check arguments
    assertthat::assert_that(is.list(anno_data))
    assertthat::assert_that(assertthat::has_name(anno_data, "transcript_df"))
    assertthat::assert_that(is.data.frame(anno_data$transcript_df))
    assertthat::assert_that(assertthat::is.count(bin_size))
    assertthat::assert_that(assertthat::is.flag(use_full_hash))

    # Prepare unspliced transcript data
    transcript_df <- anno_data$transcript_df %>%
        dplyr::filter(!.data$is_spliced)

    # Prepare transcript position GRanges object & hash IDs
    tx_granges <- methods::as(transcript_df$position, "GRanges")
    tx_chromosome <- as.character(GenomeInfoDb::seqnames(tx_granges))
    tx_strand <- as.character(BiocGenerics::strand(tx_granges))
    hash_length <- ifelse(use_full_hash, 32, 16)
    hash_id <- paste0(tx_chromosome, tx_strand) %>%
        sapply(digest::digest, algo = "md5") %>%
        substring(0, hash_length - 12) %>%
        unname()
    hash_id <- paste0("000000000000", hash_id)
    names(tx_granges) <- hash_id
    tx_start <- BiocGenerics::start(tx_granges)
    tx_end <- BiocGenerics::end(tx_granges)
    tx_start_bin <- floor(tx_start / bin_size) * bin_size
    tx_end_bin <- (floor(tx_end / bin_size) * bin_size) + bin_size

    # Prepare basic transcript information
    tx_df <- transcript_df %>%
        dplyr::transmute(
            hash_id = hash_id,
            tx_start_bin = tx_start_bin,
            tx_end_bin = tx_end_bin,
            tx_chromosome = tx_chromosome,
            tx_start = tx_start,
            tx_end = tx_end,
            tx_strand = tx_strand,
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
            gene_id = paste0(unique(.data$gene_id), collapse = ","),
            gene_name = paste0(unique(.data$gene_name), collapse = ","),
            compatible_tx = paste0(unique(.data$transcript_id), collapse = ","),
            .groups = "drop"
        )  %>%
        dplyr::transmute(
            hash_id = .data$hash_id,
            position = paste0(.data$chromosome, ":", .data$start, "-", .data$end,
                              ":", .data$strand),
            intron_positions = as.character(NA),
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

    # Prepare transcript intron & exon GRangesList object
    tx_exon_granges_list <- GenomicRanges::split(S4Vectors::unname(tx_granges),
                                                 seq_along(tx_granges))
    names(tx_exon_granges_list) <- tx_df$hash_id
    tx_intron_granges_list <- IRanges::psetdiff(tx_granges, tx_exon_granges_list)

    return(list(
        tx_df = tx_df,
        tx_granges = tx_granges,
        tx_exon_granges_list = tx_exon_granges_list,
        tx_intron_granges_list = tx_intron_granges_list
    ))
}
