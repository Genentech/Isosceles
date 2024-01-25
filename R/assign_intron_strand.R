#' Assign strand to introns
#'
#' Assigns strand to the intron GRanges object based on known intron motifs.
#'
#' @param intron_granges A GRanges object storing intron positions.
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
#' @return A copy of the intron_granges GRanges object with intron strands
#' assigned.
#' @keywords internal
assign_intron_strand <- function(intron_granges,
                                 anno_data,
                                 genome_fasta_file,
                                 min_intron_length = 30,
                                 max_intron_length = 5e6,
                                 known_intron_motifs = c("GT-AG"),
                                 rescue_annotated_introns = FALSE) {

    # Check arguments
    assertthat::assert_that(methods::is(intron_granges, "GRanges"))
    assertthat::assert_that(is.list(anno_data))
    assertthat::assert_that(assertthat::has_name(anno_data, "intron_df"))
    assertthat::assert_that(is.data.frame(anno_data$intron_df))
    assertthat::assert_that(assertthat::is.string(genome_fasta_file))
    assertthat::assert_that(file.exists(genome_fasta_file))
    assertthat::assert_that(assertthat::is.count(min_intron_length))
    assertthat::assert_that(assertthat::is.count(max_intron_length))
    assertthat::assert_that(is.character(known_intron_motifs))
    assertthat::assert_that(assertthat::is.flag(rescue_annotated_introns))

    # Read the genome sequence from the FASTA file
    genome_seq <- Biostrings::readDNAStringSet(genome_fasta_file,
                                               format = "fasta")
    names(genome_seq) <- sapply(strsplit(names(genome_seq), "\\s+"), "[", 1)

    # Prepare known intron data
    known_intron_positions <- anno_data$intron_df %>%
        dplyr::pull(.data$position) %>%
        unique()
    known_intron_positions_plus <-
        known_intron_positions[grepl(":\\+$", known_intron_positions)]
    known_intron_positions_plus <-
        unique(gsub(":\\+", "", known_intron_positions_plus))
    known_intron_positions_minus <-
        known_intron_positions[grepl(":-$", known_intron_positions)]
    known_intron_positions_minus <-
        unique(gsub(":-", "", known_intron_positions_minus))

    # Extract intron motifs (+ strand)
    intron_granges_plus <- intron_granges
    BiocGenerics::strand(intron_granges_plus) <- "+"
    intron_seq_plus <- BSgenome::getSeq(genome_seq, intron_granges_plus)
    intron_lenght_plus <- BiocGenerics::width(intron_seq_plus)
    intron_seq_plus[intron_lenght_plus < min_intron_length] <-
        strrep("N", 4)
    intron_seq_plus[intron_lenght_plus > max_intron_length] <-
        strrep("N", 4)
    intron_motif_plus <- glue::glue(
        "{XVector::subseq(intron_seq_plus, start = 1, width = 2)}", "-",
        "{XVector::subseq(intron_seq_plus, end = BiocGenerics::width(intron_seq_plus), width = 2)}"
    )

    # Extract intron motifs (- strand)
    intron_granges_minus <- intron_granges
    BiocGenerics::strand(intron_granges_minus) <- "-"
    intron_seq_minus <- BSgenome::getSeq(genome_seq, intron_granges_minus)
    intron_lenght_minus <- BiocGenerics::width(intron_seq_minus)
    intron_seq_minus[intron_lenght_minus < min_intron_length] <-
        strrep("N", 4)
    intron_seq_minus[intron_lenght_minus > max_intron_length] <-
        strrep("N", 4)
    intron_motif_minus <- glue::glue(
        "{XVector::subseq(intron_seq_minus, start = 1, width = 2)}", "-",
        "{XVector::subseq(intron_seq_minus, end = BiocGenerics::width(intron_seq_minus), width = 2)}"
    )

    # Infer intron strands
    intron_positions <- as.character(BiocGenerics::unstrand(intron_granges))
    is_strand_plus <- intron_motif_plus %in% known_intron_motifs
    is_strand_minus <- intron_motif_minus %in% known_intron_motifs
    if (rescue_annotated_introns) {
        is_rescue_plus <- intron_positions %in% known_intron_positions_plus
        is_rescue_minus <- intron_positions %in% known_intron_positions_minus
        is_rescue_both <- is_rescue_plus & is_rescue_minus
        is_rescue_plus[is_rescue_both] <- FALSE
        is_rescue_minus[is_rescue_both] <- FALSE
        is_strand_plus <- is_strand_plus | is_rescue_plus
        is_strand_minus <- is_strand_minus | is_rescue_minus
    }
    intron_strand <- rep("*", times = length(intron_granges))
    intron_strand[is_strand_plus] <- "+"
    intron_strand[is_strand_minus] <- "-"
    BiocGenerics::strand(intron_granges) <- intron_strand

    return(intron_granges)
}
