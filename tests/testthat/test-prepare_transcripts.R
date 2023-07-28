test_that("prepare_transcripts works as expected", {

    # Preparing test data
    bam_file <- system.file(
        "extdata", "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam",
        package = "Isosceles"
    )
    gtf_file <- system.file(
        "extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf",
        package = "Isosceles"
    )
    genome_fasta_file <- system.file(
        "extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa",
        package = "Isosceles"
    )
    bam_parsed <- extract_read_structures(bam_file)

    # Testing if function throws the expected errors
    expect_error(prepare_transcripts(gtf_file = NULL),
                 regexp = "gtf_file is not a string",
                 fixed = TRUE)
    expect_error(prepare_transcripts(gtf_file = "jabberwocky"),
                 regexp = "Path 'jabberwocky' does not exist",
                 fixed = TRUE)
    expect_error(prepare_transcripts(gtf_file = gtf_file,
                                     genome_fasta_file = NULL),
                 regexp = "genome_fasta_file is not a string",
                 fixed = TRUE)
    expect_error(prepare_transcripts(gtf_file = gtf_file,
                                     genome_fasta_file = "jabberwocky"),
                 regexp = "Path 'jabberwocky' does not exist",
                 fixed = TRUE)
    expect_error(prepare_transcripts(gtf_file = gtf_file,
                                     genome_fasta_file = genome_fasta_file,
                                     bam_parsed = 42),
                 regexp = "bam_parsed is not a data frame",
                 fixed = TRUE)
    expect_error(prepare_transcripts(gtf_file = gtf_file,
                                     genome_fasta_file = genome_fasta_file,
                                     bam_parsed = bam_parsed,
                                     min_intron_length = NULL),
                 regexp = "min_intron_length is not a count",
                 fixed = TRUE)
    expect_error(prepare_transcripts(gtf_file = gtf_file,
                                     genome_fasta_file = genome_fasta_file,
                                     bam_parsed = bam_parsed,
                                     known_intron_motifs = NULL),
                 regexp = "known_intron_motifs is not a character vector",
                 fixed = TRUE)
    expect_error(prepare_transcripts(gtf_file = gtf_file,
                                     genome_fasta_file = genome_fasta_file,
                                     bam_parsed = bam_parsed,
                                     rescue_annotated_introns = NULL),
                 regexp = "rescue_annotated_introns is not a flag",
                 fixed = TRUE)
    expect_error(prepare_transcripts(gtf_file = gtf_file,
                                     genome_fasta_file = genome_fasta_file,
                                     bam_parsed = bam_parsed,
                                     known_intron_granges = 42),
                 regexp = "methods::is(object = known_intron_granges, class2 =",
                 fixed = TRUE)
    expect_error(prepare_transcripts(gtf_file = gtf_file,
                                     genome_fasta_file = genome_fasta_file,
                                     bam_parsed = bam_parsed,
                                     min_bam_splice_read_count = NULL),
                 regexp = "min_bam_splice_read_count is not a count",
                 fixed = TRUE)
    expect_error(prepare_transcripts(gtf_file = gtf_file,
                                     genome_fasta_file = genome_fasta_file,
                                     bam_parsed = bam_parsed,
                                     min_bam_splice_fraction = NULL),
                 regexp = "min_bam_splice_fraction is not a numeric or integer vector",
                 fixed = TRUE)
    expect_error(prepare_transcripts(gtf_file = gtf_file,
                                     genome_fasta_file = genome_fasta_file,
                                     bam_parsed = bam_parsed,
                                     min_bam_splice_fraction = 1:10),
                 regexp = "length(min_bam_splice_fraction) not equal to 1",
                 fixed = TRUE)
    expect_error(prepare_transcripts(gtf_file = gtf_file,
                                     genome_fasta_file = genome_fasta_file,
                                     bam_parsed = bam_parsed,
                                     min_bam_splice_fraction = -42),
                 regexp = "min_bam_splice_fraction not greater than or equal to 0",
                 fixed = TRUE)
    expect_error(prepare_transcripts(gtf_file = gtf_file,
                                     genome_fasta_file = genome_fasta_file,
                                     bam_parsed = bam_parsed,
                                     min_bam_splice_fraction = 42),
                 regexp = "min_bam_splice_fraction not less than or equal to 1",
                 fixed = TRUE)
    expect_error(prepare_transcripts(gtf_file = gtf_file,
                                     genome_fasta_file = genome_fasta_file,
                                     bam_parsed = bam_parsed,
                                     bin_size = NULL),
                 regexp = "bin_size is not a count",
                 fixed = TRUE)

    # Testing if function returns the expected output (reference and BAM data)
    expect_message(
        tx_list <- prepare_transcripts(
            gtf_file = gtf_file, genome_fasta_file = genome_fasta_file,
            bam_parsed = bam_parsed
        ),
        regexp = "hash_id",
        fixed = TRUE
    )
    expect_true(is.list(tx_list))
    expect_identical(length(tx_list), 4L)
    expect_identical(names(tx_list),
                     c("tx_df", "tx_granges","tx_exon_granges_list",
                       "tx_intron_granges_list"))
    expect_identical(dim(tx_list$tx_df), c(249L, 12L))
    expect_identical(colnames(tx_list$tx_df),
                     c("transcript_id", "position", "intron_positions",
                       "gene_id", "gene_name", "compatible_gene_ids",
                       "compatible_gene_names", "compatible_tx",
                       "splicing_support_level", "fivethree_support_level",
                       "read_count", "relative_expression"))
    expect_identical(sum(is.na(tx_list$tx_df$hash_id)), 0L)
    expect_identical(sum(is.na(tx_list$tx_df$position)), 0L)
    expect_identical(sum(is.na(tx_list$tx_df$intron_positions)), 8L)
    expect_identical(sum(is.na(tx_list$tx_df$gene_id)), 12L)
    expect_identical(sum(is.na(tx_list$tx_df$gene_name)), 12L)
    expect_identical(sum(is.na(tx_list$tx_df$compatible_gene_ids)), 2L)
    expect_identical(sum(is.na(tx_list$tx_df$compatible_gene_names)), 2L)
    expect_identical(sum(is.na(tx_list$tx_df$compatible_tx)), 138L)
    expect_identical(names(table(tx_list$tx_df$splicing_support_level)),
                     c("AF", "AP", "AS", "AX", "DN", "EC", "NC", "PC"))
    expect_identical(as.numeric(table(tx_list$tx_df$splicing_support_level)),
                     c(10, 105, 53, 14, 26, 6, 5, 30))
    expect_identical(names(table(tx_list$tx_df$fivethree_support_level)),
                     c("3T", "5T", "FL", "FT", "NA"))
    expect_identical(as.numeric(table(tx_list$tx_df$fivethree_support_level)),
                     c(3, 36, 123, 9, 78))
    expect_identical(sum(is.na(tx_list$tx_df$read_count)), 105L)
    expect_identical(sum(tx_list$tx_df$read_count, na.rm = TRUE), 399L)
    expect_identical(sum(is.na(tx_list$tx_df$relative_expression)), 117L)
    expect_identical(sum(tx_list$tx_df$relative_expression, na.rm = TRUE), 5)
    expect_true(class(tx_list$tx_granges) == "GRanges")
    expect_identical(length(tx_list$tx_granges), 249L)
    expect_identical(
        names(table(BiocGenerics::strand(tx_list$tx_granges))),
        c("+", "-", "*")
    )
    expect_identical(
        as.numeric(table(BiocGenerics::strand(tx_list$tx_granges))),
        c(82, 153, 14)
    )
    expect_true(grepl("GRangesList", class(tx_list$tx_exon_granges_list)))
    expect_identical(length(tx_list$tx_exon_granges_list), 249L)
    expect_identical(length(unlist(tx_list$tx_exon_granges_list)), 1991L)
    expect_true(grepl("GRangesList", class(tx_list$tx_intron_granges_list)))
    expect_identical(length(tx_list$tx_intron_granges_list), 249L)
    expect_identical(length(unlist(tx_list$tx_intron_granges_list)), 1742L)

    # Testing if function returns the expected output (reference data only)
    expect_message(
        tx_list <- prepare_transcripts(
            gtf_file = gtf_file, genome_fasta_file = genome_fasta_file,
            bam_parsed = NULL
        ),
        regexp = "feature_id",
        fixed = TRUE
    )
    expect_true(is.list(tx_list))
    expect_identical(length(tx_list), 4L)
    expect_identical(names(tx_list),
                     c("tx_df", "tx_granges","tx_exon_granges_list",
                       "tx_intron_granges_list"))
    expect_identical(dim(tx_list$tx_df), c(105L, 12L))
    expect_identical(colnames(tx_list$tx_df),
                     c("transcript_id", "position", "intron_positions",
                       "gene_id", "gene_name", "compatible_gene_ids",
                       "compatible_gene_names", "compatible_tx",
                       "splicing_support_level", "fivethree_support_level",
                       "read_count", "relative_expression"))
    expect_identical(sum(is.na(tx_list$tx_df$hash_id)), 0L)
    expect_identical(sum(is.na(tx_list$tx_df$position)), 0L)
    expect_identical(sum(is.na(tx_list$tx_df$intron_positions)), 8L)
    expect_identical(sum(is.na(tx_list$tx_df$gene_id)), 0L)
    expect_identical(sum(is.na(tx_list$tx_df$gene_name)), 0L)
    expect_identical(sum(is.na(tx_list$tx_df$compatible_gene_ids)), 0L)
    expect_identical(sum(is.na(tx_list$tx_df$compatible_gene_names)), 0L)
    expect_identical(sum(is.na(tx_list$tx_df$compatible_tx)), 0L)
    expect_identical(names(table(tx_list$tx_df$splicing_support_level)),
                     "AP")
    expect_identical(as.numeric(table(tx_list$tx_df$splicing_support_level)),
                     105)
    expect_identical(names(table(tx_list$tx_df$fivethree_support_level)),
                     "FL")
    expect_identical(as.numeric(table(tx_list$tx_df$fivethree_support_level)),
                     105)
    expect_true(all(is.na(tx_list$tx_df$read_count)))
    expect_true(all(is.na(tx_list$tx_df$relative_expression)))
    expect_true(class(tx_list$tx_granges) == "GRanges")
    expect_identical(length(tx_list$tx_granges), 105L)
    expect_identical(
        names(table(BiocGenerics::strand(tx_list$tx_granges))),
        c("+", "-")
    )
    expect_identical(
        as.numeric(table(BiocGenerics::strand(tx_list$tx_granges))),
        c(49, 56)
    )
    expect_true(grepl("GRangesList", class(tx_list$tx_exon_granges_list)))
    expect_identical(length(tx_list$tx_exon_granges_list), 105L)
    expect_identical(length(unlist(tx_list$tx_exon_granges_list)), 831L)
    expect_true(grepl("GRangesList", class(tx_list$tx_intron_granges_list)))
    expect_identical(length(tx_list$tx_intron_granges_list), 105L)
    expect_identical(length(unlist(tx_list$tx_intron_granges_list)), 726L)
})
