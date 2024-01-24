#' Prepare transcript data extracted from BAM file(s)
#'
#' Prepares transcript data extracted from BAM file(s) for further processing.
#'
#' @param bam_parsed A data frame containing non-redundant read structure data
#' returned by the \code{\link{bam_to_read_structures}} function.
#' @param anno_data A list containing genome annotation data returned by
#' the \code{\link{prepare_reference_annotations}} function.
#' @param genome_fasta_file A string containing a genome FASTA file path.
#' @param min_intron_length An integer scalar specifying the minimal length
#' of introns to assign strand to.
#' @param max_intron_length An integer scalar specifying the maximum length
#' of introns to assign strand to.
#' @param known_intron_motifs A character vector specifying the known intron
#' motifs.
#' @param rescue_annotated_introns A logical scalar specifying if introns
#' found in genome annotations should be kept even if they don't have known
#' intron motifs.
#' @param known_intron_granges A GRanges object storing known intron positions
#' (e.g. from short read data) used for transcript classification. If set to
#' NULL, only introns from reference annotations are used.
#' @param min_bam_splice_read_count An integer scalar specifying the read count
#' threshold for splice sites confirmed by aligned reads.
#' @param min_bam_splice_fraction A numeric scalar specifying the minimum
#' connectivity fraction to a known splice site for splice sites confirmed by
#' aligned reads.
#' @return A named list containing following elements:
#' \describe{
#'   \item{tx_df}{a data frame storing extracted transcript data}
#'   \item{tx_granges}{a GRanges object storing genomic positions of extracted transcript}
#'   \item{tx_exon_granges_list}{a GRangesList object storing exon genomic positions of extracted transcript}
#'   \item{tx_intron_granges_list}{a GRangesList object storing intron genomic positions of extracted transcript}
#' }
#' @keywords internal
prepare_bam_transcripts <- function(bam_parsed,
                                    anno_data,
                                    genome_fasta_file,
                                    min_intron_length = 30,
                                    max_intron_length = 5e6,
                                    known_intron_motifs = c("GT-AG"),
                                    rescue_annotated_introns = FALSE,
                                    known_intron_granges = NULL,
                                    min_bam_splice_read_count = 2,
                                    min_bam_splice_fraction = 0.1) {

    # Check arguments
    assertthat::assert_that(is.data.frame(bam_parsed))
    assertthat::assert_that(is.list(anno_data))
    assertthat::assert_that(assertthat::has_name(anno_data, "gene_df"))
    assertthat::assert_that(is.data.frame(anno_data$gene_df))
    assertthat::assert_that(assertthat::has_name(anno_data, "intron_df"))
    assertthat::assert_that(is.data.frame(anno_data$intron_df))
    assertthat::assert_that(assertthat::has_name(anno_data, "splicing_df"))
    assertthat::assert_that(is.data.frame(anno_data$splicing_df))
    assertthat::assert_that(assertthat::has_name(anno_data, "transcript_first_last_df"))
    assertthat::assert_that(is.data.frame(anno_data$transcript_first_last_df))
    assertthat::assert_that(assertthat::is.string(genome_fasta_file))
    assertthat::assert_that(file.exists(genome_fasta_file))
    assertthat::assert_that(assertthat::is.count(min_intron_length))
    assertthat::assert_that(assertthat::is.count(max_intron_length))
    assertthat::assert_that(is.character(known_intron_motifs))
    assertthat::assert_that(assertthat::is.flag(rescue_annotated_introns))
    if (!is.null(known_intron_granges)) {
        assertthat::assert_that(methods::is(known_intron_granges, "GRanges"))
    }
    assertthat::assert_that(assertthat::is.count(min_bam_splice_read_count))
    assertthat::assert_that(is.numeric(min_bam_splice_fraction))
    assertthat::assert_that(length(min_bam_splice_fraction) == 1)
    assertthat::assert_that(min_bam_splice_fraction >= 0)
    assertthat::assert_that(min_bam_splice_fraction <= 1)

    # Helper functions
    make_hash <- function(x) {
        x %>%
            sapply(digest::digest) %>%
            substring(0, 16) %>%
            `names<-`(NULL)
    }
    median_position <- function(x) {
        x %>%
            sapply(stats::median) %>%
            as.integer()
    }

    # Process intron data
    intron_data <- process_intron_data(bam_parsed = bam_parsed,
                                       anno_data = anno_data,
                                       genome_fasta_file = genome_fasta_file,
                                       min_intron_length = min_intron_length,
                                       max_intron_length = max_intron_length,
                                       known_intron_motifs = known_intron_motifs,
                                       rescue_annotated_introns = rescue_annotated_introns)
    nr_intron_positions <- intron_data$nr_intron_positions
    nr_intron_strand <- methods::as(nr_intron_positions, "GRanges") %>%
        BiocGenerics::strand() %>%
        as.character()
    tx_intron_idx <- intron_data$tx_intron_idx
    tx_anno_intron_idx <- intron_data$tx_anno_intron_idx
    tx_intron_positions <- lapply(tx_intron_idx, function(x) {
        nr_intron_positions[x]
    })
    tx_intron_positions <- sapply(tx_intron_positions, paste0, collapse = ",")

    # Prepare basic transcript information
    tx_start <- median_position(bam_parsed$start_positions)
    tx_end <- median_position(bam_parsed$end_positions)
    tx_intron_strand <- unname(lapply(tx_intron_idx, function(x) {
        nr_intron_strand[x]
    }))
    tx_strand <- sapply(tx_intron_strand, function(x) {
        paste0(unique(x), collapse = ",")
    })
    tx_strand <- ifelse(tx_strand %in% c("+", "-"), tx_strand, "*")
    tx_position <- glue::glue("{bam_parsed$chromosome}:{tx_start}-{tx_end}:{tx_strand}")
    tx_read_counts <- bam_parsed$read_count
    tx_df <- data.frame(
        hash_id = make_hash(tx_intron_positions),
        read_count = tx_read_counts,
        position = as.character(tx_position),
        intron_positions = tx_intron_positions
    )
    assertthat::assert_that(
        length(tx_df$hash_id) == length(unique(tx_df$hash_id))
    )

    # Assign transcripts to genes
    tx_gene_df <- assign_transcript_gene(tx_df, anno_data)
    tx_df <- cbind(tx_df, tx_gene_df[, -1])

    # Assign matching annotated transcript IDs
    tx_splicing_df <- anno_data$splicing_df %>%
        dplyr::transmute(
            hash_id = make_hash(.data$intron_positions),
            compatible_tx = .data$compatible_tx
        )
    tx_df <- tx_df %>%
        dplyr::left_join(tx_splicing_df)

    # Calculate splicing support level
    tx_df$splicing_support_level <- check_splicing_status(
        tx_df = tx_df, anno_data = anno_data, bam_parsed = bam_parsed,
        nr_intron_positions = nr_intron_positions, tx_intron_idx = tx_intron_idx,
        tx_anno_intron_idx = tx_anno_intron_idx,
        known_intron_granges = known_intron_granges,
        min_bam_splice_read_count = min_bam_splice_read_count,
        min_bam_splice_fraction = min_bam_splice_fraction
    )

    # Prepare transcript position GRanges object
    tx_granges <- methods::as(tx_df$position, "GRanges")
    names(tx_granges) <- tx_df$hash_id

    # Prepare transcript intron & exon GRangesList object
    tx_intron_granges_df <- tx_intron_idx %>%
        tibble::enframe() %>%
        dplyr::rename(tx_idx = "name", intron_idx = "value") %>%
        tidyr::unchop("intron_idx")
    tx_intron_granges <- methods::as(
        nr_intron_positions[tx_intron_granges_df$intron_idx], "GRanges"
    )
    tx_intron_granges_list <- GenomicRanges::split(
        tx_intron_granges, tx_intron_granges_df$tx_idx
    )
    names(tx_intron_granges_list) <- tx_df$hash_id
    tx_exon_granges_list <- IRanges::psetdiff(tx_granges, tx_intron_granges_list)

    # Calculate fivethree support level
    tx_df$fivethree_support_level <- check_truncation_status(
        tx_df,tx_exon_granges_list, tx_intron_granges_list, anno_data
    )

    # Calculate transcript relative expression
    tx_df$relative_expression <- tx_df %>%
        dplyr::group_by(.data$gene_id) %>%
        dplyr::mutate(rel_expr = .data$read_count / sum(.data$read_count)) %>%
        dplyr::pull(.data$rel_expr)
    tx_df$relative_expression[is.na(tx_df$gene_id)] <- NA

    # Prepare the output data
    tx_df <- tx_df[, c("hash_id", "position", "intron_positions", "gene_id",
                       "gene_name", "compatible_gene_ids",
                       "compatible_gene_names", "compatible_tx",
                       "splicing_support_level", "fivethree_support_level",
                       "read_count", "relative_expression")]

    return(list(
        tx_df = tx_df,
        tx_granges = tx_granges,
        tx_exon_granges_list = tx_exon_granges_list,
        tx_intron_granges_list = tx_intron_granges_list
    ))
}
