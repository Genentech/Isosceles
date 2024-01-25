#' Prepare intron data extracted from BAM file(s) and reference annotations
#' for further processing
#' @noRd
process_intron_data <- function(bam_parsed,
                                anno_data,
                                genome_fasta_file,
                                min_intron_length = 30,
                                max_intron_length = 5e6,
                                known_intron_motifs = c("GT-AG"),
                                rescue_annotated_introns = FALSE) {

    # Check arguments
    assertthat::assert_that(is.data.frame(bam_parsed))
    assertthat::assert_that(is.list(anno_data))
    assertthat::assert_that(assertthat::has_name(anno_data, "intron_df"))
    assertthat::assert_that(is.data.frame(anno_data$intron_df))
    assertthat::assert_that(assertthat::has_name(anno_data, "splicing_df"))
    assertthat::assert_that(is.data.frame(anno_data$splicing_df))
    assertthat::assert_that(assertthat::is.string(genome_fasta_file))
    assertthat::assert_that(file.exists(genome_fasta_file))
    assertthat::assert_that(assertthat::is.count(min_intron_length))
    assertthat::assert_that(assertthat::is.count(max_intron_length))
    assertthat::assert_that(is.character(known_intron_motifs))
    assertthat::assert_that(assertthat::is.flag(rescue_annotated_introns))

    # Prepare a non-redundant annotated intron set
    nr_anno_intron_positions <- unique(anno_data$intron_df$position)
    nr_anno_intron_granges <- methods::as(nr_anno_intron_positions, "GRanges")

    # Prepare a non-redundant intron set from extracted read data
    nr_unstranded_intron_positions <- bam_parsed$intron_positions %>%
        strsplit(",") %>%
        unlist() %>%
        unique()
    nr_intron_granges <- methods::as(nr_unstranded_intron_positions, "GRanges") %>%
        GenomicRanges::sort()
    nr_unstranded_intron_positions <- as.character(nr_intron_granges)
    nr_intron_granges <- assign_intron_strand(nr_intron_granges,
                                              anno_data = anno_data,
                                              genome_fasta_file = genome_fasta_file,
                                              min_intron_length = min_intron_length,
                                              max_intron_length = max_intron_length,
                                              known_intron_motifs = known_intron_motifs,
                                              rescue_annotated_introns = rescue_annotated_introns)
    nr_intron_positions <- as.character(nr_intron_granges)

    # Process transcript intron structures
    # (convert to integer indices & assign strand to introns)
    tx_intron_idx <- lapply(strsplit(bam_parsed$intron_positions, ","),
                            fastmatch::fmatch, nr_unstranded_intron_positions)
    tx_intron_positions <- lapply(tx_intron_idx, function(x) {
        nr_intron_positions[x]
    })

    # Add annotated introns to the non-redundant intron set
    nr_intron_granges <- c(nr_intron_granges, nr_anno_intron_granges) %>%
        BiocGenerics::unique() %>%
        GenomicRanges::sort()
    nr_intron_positions <- as.character(nr_intron_granges)

    # Recalculate integer indices of transcript intron structures
    tx_intron_idx <- lapply(tx_intron_positions,
                            fastmatch::fmatch, nr_intron_positions)
    tx_intron_positions <- sapply(tx_intron_positions, paste0, collapse = ",")

    # Calculate integer indices of annotated transcript intron structures
    tx_anno_intron_idx <- lapply(strsplit(anno_data$splicing_df$intron_positions, ","),
                                 fastmatch::fmatch, nr_intron_positions)

    return(list(
        nr_intron_positions = as.character(nr_intron_positions),
        tx_intron_idx = tx_intron_idx,
        tx_anno_intron_idx = tx_anno_intron_idx
    ))
}
