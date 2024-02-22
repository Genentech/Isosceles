#' Create an intron BED file from a GTF file
#'
#' Creates an intron BED file from a GTF file.
#'
#' @param gtf_file A string containing a GTF file path.
#' @param file A string specifying the output file path.
#' @return Nothing is returned.
#' @export
gtf_to_intron_bed <- function(gtf_file,
                              file) {

    # Check arguments
    assertthat::assert_that(assertthat::is.string(gtf_file))
    assertthat::assert_that(file.exists(gtf_file))
    assertthat::assert_that(assertthat::is.string(file))

    # Read annotation data from the GTF file
    txdb <- suppressWarnings(
        GenomicFeatures::makeTxDbFromGFF(gtf_file, format = "gtf")
    )

    # Prepare intron GRanges data
    transcript_granges <- GenomicFeatures::transcripts(
        txdb, columns = c("tx_name", "gene_id")
    )
    exon_granges_list <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
    exon_granges_list <- exon_granges_list[S4Vectors::mcols(transcript_granges)$tx_name]
    tx_exon_strands <- sapply(exon_granges_list, function(granges) {
        unique(BiocGenerics::strand(granges))
    })
    original_tx_strands <- BiocGenerics::strand(transcript_granges)
    BiocGenerics::strand(transcript_granges) <- tx_exon_strands
    intron_granges_list <- IRanges::psetdiff(transcript_granges, exon_granges_list)
    intron_granges <- BiocGenerics::unique(unlist(intron_granges_list))
    names(intron_granges) <- NULL

    # Give a warning if there are any issues with the input GTF file
    wrong_strand_tx_ids <- names(tx_exon_strands)[
        as.vector(tx_exon_strands != original_tx_strands)
    ]
    if (length(wrong_strand_tx_ids > 0)) {
        warning(paste0("WARNING: the following transcripts have different ",
                       "strand than their exons: ",
                       paste(wrong_strand_tx_ids, collapse = ", "),
                       "\n",
                       "The transcript strand values will be corrected."))
    }

    # Save the BED file
    rtracklayer::export(intron_granges, file, "bed")

    # Return nothing
    invisible()
}
