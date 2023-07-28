#' Transcript data preparation
#'
#' Prepare transcript data (reference and extracted from the BAM files) for
#' further analysis
#'
#' @param gtf_file A string containing a GTF file path
#' @param genome_fasta_file A string containing a genome FASTA file path
#' @param bam_parsed A data frame containing non-redundant read structure data
#' returned by the \code{\link{extract_read_structures}} function. If NULL,
#' only reference transcripts are used
#' @param min_intron_length An integer scalar specifying the minimal length
#' of introns to assign strand to
#' @param known_intron_motifs A character vector specifying the known intron
#' motifs
#' @param rescue_annotated_introns A logical scalar specifying if introns
#' found in genome annotations should be kept even if they don't have known
#' intron motifs
#' @param known_intron_granges A GRanges object storing known intron positions
#' (e.g. from short read data) used for transcript classification. If set to
#' NULL, only introns from reference annotations are used
#' @param min_bam_splice_read_count An integer scalar specifying the read count
#' threshold for splice sites confirmed by aligned reads
#' @param min_bam_splice_fraction A numeric scalar specifying the minimum
#' connectivity fraction to a known splice site for splice sites confirmed by
#' aligned reads
#' @param bin_size An integer scalar specifying the bin size for transcript
#' start and end position binning
#' @return A named list containing following elements:
#' \describe{
#'   \item{tx_df}{a data frame storing extracted transcript data}
#'   \item{tx_granges}{a GRanges object storing genomic positions of extracted transcript}
#'   \item{tx_exon_granges_list}{a GRangesList object storing exon genomic positions of extracted transcript}
#'   \item{tx_intron_granges_list}{a GRangesList object storing intron genomic positions of extracted transcript}
#' }
#' @export
prepare_transcripts <- function(gtf_file,
                                genome_fasta_file,
                                bam_parsed,
                                min_intron_length = 30,
                                known_intron_motifs = c("GT-AG"),
                                rescue_annotated_introns = FALSE,
                                known_intron_granges = NULL,
                                min_bam_splice_read_count = 2,
                                min_bam_splice_fraction = 0.1,
                                bin_size = 50) {

    # Check arguments
    assertthat::assert_that(assertthat::is.string(gtf_file))
    assertthat::assert_that(file.exists(gtf_file))
    assertthat::assert_that(assertthat::is.string(genome_fasta_file))
    assertthat::assert_that(file.exists(genome_fasta_file))
    if (!is.null(bam_parsed)) {
        assertthat::assert_that(is.data.frame(bam_parsed))
    }
    assertthat::assert_that(assertthat::is.count(min_intron_length))
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
    assertthat::assert_that(assertthat::is.count(bin_size))

    # Prepare reference annotation data
    anno_data <- prepare_reference_annotations(gtf_file)

    # Prepare spliced reference transcripts
    tx_list_spliced <- prepare_reference_spliced_transcripts(
        anno_data, bin_size = bin_size
    )
    tx_df <- tx_list_spliced$tx_df
    tx_granges <- tx_list_spliced$tx_granges
    tx_exon_granges_list <- tx_list_spliced$tx_exon_granges_list
    tx_intron_granges_list <- tx_list_spliced$tx_intron_granges_list

    # Prepare unspliced reference transcripts
    if (any(!anno_data$transcript_df$is_spliced)) {
        tx_list_unspliced <- prepare_reference_unspliced_transcripts(
            anno_data, bin_size = bin_size
        )
        tx_df <- rbind(tx_df, tx_list_unspliced$tx_df)
        tx_granges <- suppressWarnings(c(tx_granges, tx_list_unspliced$tx_granges))
        tx_exon_granges_list <- suppressWarnings(
            c(tx_exon_granges_list, tx_list_unspliced$tx_exon_granges_list)
        )
        tx_intron_granges_list <- suppressWarnings(
            c(tx_intron_granges_list, tx_list_unspliced$tx_intron_granges_list)
        )
    }

    # Prepare transcripts from the BAM files
    if (!is.null(bam_parsed)) {
        tx_list_bam <- prepare_bam_transcripts(
            bam_parsed = bam_parsed, anno_data = anno_data,
            genome_fasta_file = genome_fasta_file, min_intron_length = min_intron_length,
            known_intron_motifs = known_intron_motifs,
            rescue_annotated_introns = rescue_annotated_introns,
            known_intron_granges = known_intron_granges,
            min_bam_splice_read_count = min_bam_splice_read_count,
            min_bam_splice_fraction = min_bam_splice_fraction
        )
        tx_df <- rbind(tx_df, tx_list_bam$tx_df)
        tx_granges <- suppressWarnings(c(tx_granges, tx_list_bam$tx_granges))
        tx_exon_granges_list <- suppressWarnings(
            c(tx_exon_granges_list, tx_list_bam$tx_exon_granges_list)
        )
        tx_intron_granges_list <- suppressWarnings(
            c(tx_intron_granges_list, tx_list_bam$tx_intron_granges_list)
        )
    }

    # Reorder transcripts by gene ID
    tx_order_idx <- order(tx_df$gene_id)
    tx_df <- tx_df[tx_order_idx,]
    rownames(tx_df) <- NULL
    tx_granges <- tx_granges[tx_order_idx]
    tx_exon_granges_list <- tx_exon_granges_list[tx_order_idx]
    tx_intron_granges_list <- tx_intron_granges_list[tx_order_idx]

    # Prepare transcript IDs
    tx_start_bin <- floor(BiocGenerics::start(tx_granges) / bin_size) * bin_size
    tx_end_bin <- (floor(BiocGenerics::end(tx_granges) / bin_size) * bin_size) + bin_size
    transcript_id <- paste0(
        "ISOT", gsub("(.{4})", "-\\1", tx_df$hash_id), ":s", tx_start_bin,
        ":e", tx_end_bin, ":", tx_df$splicing_support_level, ":",
        tx_df$fivethree_support_level
    )
    assertthat::assert_that(
        length(transcript_id) == length(unique(transcript_id))
    )
    tx_df <- dplyr::rename(tx_df, transcript_id = "hash_id")
    tx_df$transcript_id <- transcript_id
    names(tx_granges) <- transcript_id
    names(tx_exon_granges_list) <- transcript_id
    names(tx_intron_granges_list) <- transcript_id

    return(list(
        tx_df = tx_df,
        tx_granges = tx_granges,
        tx_exon_granges_list = tx_exon_granges_list,
        tx_intron_granges_list = tx_intron_granges_list
    ))
}
