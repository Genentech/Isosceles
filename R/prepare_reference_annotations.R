#' GTF reference annotation data processing
#'
#' Extract required reference annotation from a GTF file
#'
#' @param gtf_file A string containing a GTF file path
#' @return A named list containing following elements:
#' \describe{
#'   \item{gene_df}{a data frame storing annotated gene data}
#'   \item{transcript_df}{a data frame storing annotated transcript data}
#'   \item{splicing_df}{a data frame storing non-redundant annotated transcript intron structure data}
#'   \item{intron_df}{a data frame storing annotated intron data}
#'   \item{transcript_first_last_df}{a data frame storing confirmed transcript first/last exon/intron positions}
#' }
#' @keywords internal
prepare_reference_annotations <- function(gtf_file) {

    # Check arguments
    assertthat::assert_that(assertthat::is.string(gtf_file))
    assertthat::assert_that(file.exists(gtf_file))

    # Read annotation data from the GTF file
    txdb <- suppressWarnings(
        GenomicFeatures::makeTxDbFromGFF(gtf_file, format = "gtf")
    )
    gtf_granges <- rtracklayer::import(gtf_file, format = "gtf")
    exon_granges <- gtf_granges[S4Vectors::mcols(gtf_granges)$type == "exon"]
    if (is.null(S4Vectors::mcols(exon_granges)$gene_name)) {
        S4Vectors::mcols(exon_granges)$gene_name <-
            S4Vectors::mcols(exon_granges)$gene_id
    }

    # Prepare gene data
    gene_df <- S4Vectors::mcols(exon_granges) %>%
        as.data.frame() %>%
        dplyr::select("gene_id", "gene_name") %>%
        dplyr::distinct()

    # Prepare intron data
    intron_granges <- GenomicFeatures::tidyIntrons(txdb)
    intron_df <- data.frame(
        position = as.character(intron_granges),
        transcript_id = S4Vectors::mcols(intron_granges)$tx_name,
        gene_id = S4Vectors::mcols(intron_granges)$gene_id
    )

    # Prepare transcript data
    transcript_granges <- GenomicFeatures::transcripts(
        txdb, columns = c("tx_name", "gene_id")
    )
    transcript_position_df <- data.frame(
        transcript_id = S4Vectors::mcols(transcript_granges)$tx_name,
        tx_chromosome = as.character(GenomeInfoDb::seqnames(transcript_granges)),
        tx_start = BiocGenerics::start(transcript_granges),
        tx_end = BiocGenerics::end(transcript_granges),
        tx_strand = as.character(BiocGenerics::strand(transcript_granges)),
        gene_id = unlist(S4Vectors::mcols(transcript_granges)$gene_id)
    ) %>%
        dplyr::left_join(gene_df)
    transcript_splicing_df <- intron_df %>%
        dplyr::group_by(.data$transcript_id) %>%
        dplyr::summarise(
            intron_positions = paste0(.data$position, collapse = ",")
        )
    transcript_df <- data.frame(
        transcript_id = S4Vectors::mcols(transcript_granges)$tx_name,
        position = as.character(transcript_granges),
        is_spliced = S4Vectors::mcols(transcript_granges)$tx_name %in%
            intron_df$transcript_id
    ) %>%
        dplyr::left_join(transcript_splicing_df) %>%
        tidyr::replace_na(list(intron_positions = ""))
    transcript_df$gene_id <- unlist(S4Vectors::mcols(transcript_granges)$gene_id)
    transcript_df <- transcript_df %>%
        dplyr::left_join(gene_df)

    # Merge gene IDs when genes share introns
    merged_gene_df <- merge_annotated_genes(gene_df, intron_df)
    merged_gene_ids <- stats::setNames(merged_gene_df$merged_gene_id,
                                       merged_gene_df$gene_id)
    merged_gene_names <- stats::setNames(merged_gene_df$merged_gene_name,
                                         merged_gene_df$gene_id)
    gene_df$gene_name <- merged_gene_names[gene_df$gene_id]
    gene_df$gene_id <- merged_gene_ids[gene_df$gene_id]
    gene_df <- dplyr::distinct(gene_df)
    intron_df$gene_id <- merged_gene_ids[intron_df$gene_id]
    transcript_position_df$gene_name <- merged_gene_names[transcript_position_df$gene_id]
    transcript_position_df$gene_id <- merged_gene_ids[transcript_position_df$gene_id]
    transcript_df$gene_name <- merged_gene_names[transcript_df$gene_id]
    transcript_df$gene_id <- merged_gene_ids[transcript_df$gene_id]

    # Prepare splicing data
    splicing_df <- transcript_splicing_df %>%
        dplyr::left_join(transcript_position_df) %>%
        dplyr::group_by(.data$intron_positions) %>%
        dplyr::summarise(
            chromosome = unique(.data$tx_chromosome),
            start = min(.data$tx_start),
            end = max(.data$tx_end),
            strand = unique(.data$tx_strand),
            gene_id = unique(.data$gene_id),
            gene_name = unique(.data$gene_name),
            compatible_tx = paste0(unique(.data$transcript_id), collapse = ",")
        ) %>%
        dplyr::transmute(
            intron_positions = .data$intron_positions,
            position = paste0(.data$chromosome, ":", .data$start, "-", .data$end,
                              ":", .data$strand),
            gene_id = .data$gene_id,
            gene_name = .data$gene_name,
            compatible_tx = .data$compatible_tx
        ) %>%
        as.data.frame()

    # Prepare transcript first/last exon/intron data
    transcript_exon_first_last_df <-
        GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE) %>%
        get_first_last_grange() %>%
        dplyr::rename(transcript_id = "feature_id",
                      first_exon_ref = "first_grange",
                      last_exon_ref = "last_grange")
    transcript_intron_first_last_df <-
        GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE) %>%
        get_first_last_grange() %>%
        dplyr::rename(transcript_id = "feature_id",
                      first_intron_ref = "first_grange",
                      last_intron_ref = "last_grange")
    transcript_first_last_df <- transcript_df %>%
        dplyr::select("transcript_id", "gene_id") %>%
        dplyr::inner_join(transcript_exon_first_last_df) %>%
        dplyr::inner_join(transcript_intron_first_last_df)

    return(list(
        gene_df = gene_df,
        transcript_df = transcript_df,
        splicing_df = splicing_df,
        intron_df = intron_df,
        transcript_first_last_df = transcript_first_last_df
    ))
}
