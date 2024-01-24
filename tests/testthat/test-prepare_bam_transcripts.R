test_that("prepare_bam_transcripts works as expected", {

    # Preparing test data
    bam_file <- system.file(
        "extdata", "bulk_rnaseq.bam",
        package = "Isosceles"
    )
    gtf_file <- system.file(
        "extdata", "bulk_rnaseq.gtf",
        package = "Isosceles"
    )
    genome_fasta_file <- system.file(
        "extdata", "bulk_rnaseq.fa",
        package = "Isosceles"
    )
    anno_data <- prepare_reference_annotations(gtf_file)
    bam_parsed <- bam_to_read_structures(bam_file)

    # Testing if function throws the expected errors
    expect_error(prepare_bam_transcripts(bam_parsed = NULL),
                 regexp = "bam_parsed is not a data frame",
                 fixed = TRUE)
    expect_error(prepare_bam_transcripts(bam_parsed = bam_parsed,
                                         anno_data = NULL),
                 regexp = "anno_data is not a list",
                 fixed = TRUE)
    expect_error(prepare_bam_transcripts(bam_parsed = bam_parsed,
                                         anno_data = list()),
                 regexp = "anno_data does not have all of these name(s): 'gene_df'",
                 fixed = TRUE)
    expect_error(prepare_bam_transcripts(bam_parsed = bam_parsed,
                                         anno_data = list(
                                             gene_df = 42
                                         )),
                 regexp = "anno_data$gene_df is not a data frame",
                 fixed = TRUE)
    expect_error(prepare_bam_transcripts(bam_parsed = bam_parsed,
                                         anno_data = list(
                                             gene_df = anno_data$gene_df
                                         )),
                 regexp = "anno_data does not have all of these name(s): 'intron_df'",
                 fixed = TRUE)
    expect_error(prepare_bam_transcripts(bam_parsed = bam_parsed,
                                         anno_data = list(
                                             gene_df = anno_data$gene_df,
                                             intron_df = 42
                                         )),
                 regexp = "anno_data$intron_df is not a data frame",
                 fixed = TRUE)
    expect_error(prepare_bam_transcripts(bam_parsed = bam_parsed,
                                         anno_data = list(
                                             gene_df = anno_data$gene_df,
                                             intron_df = anno_data$intron_df
                                         )),
                 regexp = "anno_data does not have all of these name(s): 'splicing_df'",
                 fixed = TRUE)
    expect_error(prepare_bam_transcripts(bam_parsed = bam_parsed,
                                         anno_data = list(
                                             gene_df = anno_data$gene_df,
                                             intron_df = anno_data$intron_df,
                                             splicing_df = 42
                                         )),
                 regexp = "anno_data$splicing_df is not a data frame",
                 fixed = TRUE)
    expect_error(prepare_bam_transcripts(bam_parsed = bam_parsed,
                                         anno_data = list(
                                             gene_df = anno_data$gene_df,
                                             intron_df = anno_data$intron_df,
                                             splicing_df = anno_data$splicing_df
                                         )),
                 regexp = "anno_data does not have all of these name(s): 'transcript_first_last_df'",
                 fixed = TRUE)
    expect_error(prepare_bam_transcripts(bam_parsed = bam_parsed,
                                         anno_data = list(
                                             gene_df = anno_data$gene_df,
                                             intron_df = anno_data$intron_df,
                                             splicing_df = anno_data$splicing_df,
                                             transcript_first_last_df = 42
                                         )),
                 regexp = "anno_data$transcript_first_last_df is not a data frame",
                 fixed = TRUE)
    expect_error(prepare_bam_transcripts(bam_parsed = bam_parsed,
                                         anno_data = anno_data,
                                         genome_fasta_file = NULL),
                 regexp = "genome_fasta_file is not a string",
                 fixed = TRUE)
    expect_error(prepare_bam_transcripts(bam_parsed = bam_parsed,
                                         anno_data = anno_data,
                                         genome_fasta_file = "jabberwocky"),
                 regexp = "Path 'jabberwocky' does not exist",
                 fixed = TRUE)
    expect_error(prepare_bam_transcripts(bam_parsed = bam_parsed,
                                         anno_data = anno_data,
                                         genome_fasta_file = genome_fasta_file,
                                         min_intron_length = NULL),
                 regexp = "min_intron_length is not a count",
                 fixed = TRUE)
    expect_error(prepare_bam_transcripts(bam_parsed = bam_parsed,
                                         anno_data = anno_data,
                                         genome_fasta_file = genome_fasta_file,
                                         max_intron_length = NULL),
                 regexp = "max_intron_length is not a count",
                 fixed = TRUE)
    expect_error(prepare_bam_transcripts(bam_parsed = bam_parsed,
                                         anno_data = anno_data,
                                         genome_fasta_file = genome_fasta_file,
                                         known_intron_motifs = NULL),
                 regexp = "known_intron_motifs is not a character vector",
                 fixed = TRUE)
    expect_error(prepare_bam_transcripts(bam_parsed = bam_parsed,
                                         anno_data = anno_data,
                                         genome_fasta_file = genome_fasta_file,
                                         rescue_annotated_introns = NULL),
                 regexp = "rescue_annotated_introns is not a flag",
                 fixed = TRUE)
    expect_error(prepare_bam_transcripts(bam_parsed = bam_parsed,
                                         anno_data = anno_data,
                                         genome_fasta_file = genome_fasta_file,
                                         known_intron_granges = 42),
                 regexp = "methods::is(object = known_intron_granges, class2 = ",
                 fixed = TRUE)
    expect_error(prepare_bam_transcripts(bam_parsed = bam_parsed,
                                         anno_data = anno_data,
                                         genome_fasta_file = genome_fasta_file,
                                         min_bam_splice_read_count = NULL),
                 regexp = "min_bam_splice_read_count is not a count",
                 fixed = TRUE)
    expect_error(prepare_bam_transcripts(bam_parsed = bam_parsed,
                                         anno_data = anno_data,
                                         genome_fasta_file = genome_fasta_file,
                                         min_bam_splice_fraction = NULL),
                 regexp = "min_bam_splice_fraction is not a numeric or integer vector",
                 fixed = TRUE)
    expect_error(prepare_bam_transcripts(bam_parsed = bam_parsed,
                                         anno_data = anno_data,
                                         genome_fasta_file = genome_fasta_file,
                                         min_bam_splice_fraction = 1:10),
                 regexp = "length(min_bam_splice_fraction) not equal to 1",
                 fixed = TRUE)
    expect_error(prepare_bam_transcripts(bam_parsed = bam_parsed,
                                         anno_data = anno_data,
                                         genome_fasta_file = genome_fasta_file,
                                         min_bam_splice_fraction = -42),
                 regexp = "min_bam_splice_fraction not greater than or equal to 0",
                 fixed = TRUE)
    expect_error(prepare_bam_transcripts(bam_parsed = bam_parsed,
                                         anno_data = anno_data,
                                         genome_fasta_file = genome_fasta_file,
                                         min_bam_splice_fraction = 42),
                 regexp = "min_bam_splice_fraction not less than or equal to 1",
                 fixed = TRUE)

    # Testing if function returns the expected output
    expect_message(
        tx_list <- prepare_bam_transcripts(
            bam_parsed = bam_parsed, anno_data = anno_data,
            genome_fasta_file = genome_fasta_file
        ),
        regexp = "hash_id",
        fixed = TRUE
    )
    expect_true(is.list(tx_list))
    expect_identical(length(tx_list), 4L)
    expect_identical(names(tx_list),
                     c("tx_df", "tx_granges","tx_exon_granges_list",
                       "tx_intron_granges_list"))
    expect_identical(dim(tx_list$tx_df), c(74L, 12L))
    expect_identical(colnames(tx_list$tx_df),
                     c("hash_id", "position", "intron_positions",
                       "gene_id", "gene_name", "compatible_gene_ids",
                       "compatible_gene_names", "compatible_tx",
                       "splicing_support_level", "fivethree_support_level",
                       "read_count", "relative_expression"))
    expect_identical(sum(is.na(tx_list$tx_df$hash_id)), 0L)
    expect_identical(sum(is.na(tx_list$tx_df$position)), 0L)
    expect_identical(sum(is.na(tx_list$tx_df$intron_positions)), 0L)
    expect_identical(sum(is.na(tx_list$tx_df$gene_id)), 11L)
    expect_identical(sum(is.na(tx_list$tx_df$gene_name)), 11L)
    expect_identical(sum(is.na(tx_list$tx_df$compatible_gene_ids)), 1L)
    expect_identical(sum(is.na(tx_list$tx_df$compatible_gene_names)), 1L)
    expect_identical(sum(is.na(tx_list$tx_df$compatible_tx)), 69L)
    expect_identical(names(table(tx_list$tx_df$splicing_support_level)),
                     c("AF", "AS", "AX", "DN", "EC", "NC", "PC"))
    expect_identical(as.numeric(table(tx_list$tx_df$splicing_support_level)),
                     c(10, 24, 4, 15, 1, 1, 19))
    expect_identical(names(table(tx_list$tx_df$fivethree_support_level)),
                     c("3T", "5T", "FL", "FT", "NA"))
    expect_identical(as.numeric(table(tx_list$tx_df$fivethree_support_level)),
                     c(1, 17, 13, 4, 39))
    expect_identical(sum(tx_list$tx_df$read_count), 114L)
    expect_identical(sum(is.na(tx_list$tx_df$relative_expression)), 11L)
    expect_identical(sum(tx_list$tx_df$relative_expression, na.rm = TRUE), 4)
    expect_true(class(tx_list$tx_granges) == "GRanges")
    expect_identical(length(tx_list$tx_granges), 74L)
    expect_identical(
        names(table(BiocGenerics::strand(tx_list$tx_granges))),
        c("+", "-", "*")
    )
    expect_identical(
        as.numeric(table(BiocGenerics::strand(tx_list$tx_granges))),
        c(30, 40, 4)
    )
    expect_true(grepl("GRangesList", class(tx_list$tx_exon_granges_list)))
    expect_identical(length(tx_list$tx_exon_granges_list), 74L)
    expect_identical(length(unlist(tx_list$tx_exon_granges_list)), 521L)
    expect_true(grepl("GRangesList", class(tx_list$tx_intron_granges_list)))
    expect_identical(length(tx_list$tx_intron_granges_list), 74L)
    expect_identical(length(unlist(tx_list$tx_intron_granges_list)), 447L)
})
