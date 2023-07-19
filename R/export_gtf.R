#' Data export to a GTF file
#'
#' Export transcripts from a SummarizedExperiment to a GTF file
#'
#' @param se A transcript-level SummarizedExperiment object returned by the
#' \code{\link{prepare_transcript_se}} function
#' @param file A string specifying the output file path
#' @return Nothing is returned
#' @export
export_gtf <- function(se,
                       file) {

    # Check arguments
    assertthat::assert_that(methods::is(se, "SummarizedExperiment"))
    assertthat::assert_that(methods::is(
        SummarizedExperiment::rowRanges(se), "GRangesList"
    ))
    assertthat::assert_that(assertthat::is.string(file))

    # Prepare GTF file exon GRanges
    gtf_data <- SummarizedExperiment::rowData(se)
    gtf_ranges <- SummarizedExperiment::rowRanges(se)
    gtf_granges_nrow <- S4Vectors::elementNROWS(gtf_ranges)
    gtf_ranges <- unlist(gtf_ranges)
    S4Vectors::mcols(gtf_ranges)$type <- "exon"
    S4Vectors::mcols(gtf_ranges)$gene_id <-
        rep(gtf_data$gene_id, gtf_granges_nrow)
    S4Vectors::mcols(gtf_ranges)$transcript_id <-
        rep(gtf_data$transcript_id, gtf_granges_nrow)
    S4Vectors::mcols(gtf_ranges)$gene_name <-
        rep(gtf_data$gene_name, gtf_granges_nrow)
    S4Vectors::mcols(gtf_ranges)$compatible_gene_ids <-
        rep(gtf_data$compatible_gene_ids, gtf_granges_nrow)
    S4Vectors::mcols(gtf_ranges)$compatible_gene_names <-
        rep(gtf_data$compatible_gene_names, gtf_granges_nrow)
    S4Vectors::mcols(gtf_ranges)$compatible_tx <-
        rep(gtf_data$compatible_tx, gtf_granges_nrow)
    S4Vectors::mcols(gtf_ranges)$splicing_support_level <-
        rep(gtf_data$splicing_support_level, gtf_granges_nrow)
    S4Vectors::mcols(gtf_ranges)$fivethree_support_level <-
        rep(gtf_data$fivethree_support_level, gtf_granges_nrow)
    names(gtf_ranges) <- NULL

    # Save the GTF file
    unlink(file)
    rtracklayer::export(gtf_ranges, file, "gtf", append = TRUE)

    # Return nothing
    invisible()
}
