test_that("plot_splicing_support_levels works as expected", {

    # Preparing test data
    bam_file <- system.file(
        "extdata", "bulk_rnaseq.bam",
        package = "Isosceles"
    )
    bam_files <- c(Sample = bam_file)
    gtf_file <- system.file(
        "extdata", "bulk_rnaseq.gtf",
        package = "Isosceles"
    )
    genome_fasta_file <- system.file(
        "extdata", "bulk_rnaseq.fa",
        package = "Isosceles"
    )
    bam_parsed <- bam_to_read_structures(bam_files)
    transcript_data <- prepare_transcripts(
        gtf_file = gtf_file, genome_fasta_file = genome_fasta_file,
        bam_parsed = bam_parsed, min_bam_splice_read_count = 2,
        min_bam_splice_fraction = 0.01
    )
    transcript_data_ref_only <- prepare_transcripts(
        gtf_file = gtf_file, genome_fasta_file = genome_fasta_file,
        bam_parsed = NULL
    )

    # Testing if function throws the expected errors
    expect_error(plot_splicing_support_levels(transcript_data = NULL),
                 regexp = "transcript_data is not a list",
                 fixed = TRUE)
    expect_error(plot_splicing_support_levels(transcript_data = list()),
                 regexp = "names(transcript_data) not identical to",
                 fixed = TRUE)
    expect_error(
        plot_splicing_support_levels(transcript_data = list(
            tx_df = 42,
            tx_granges = transcript_data$tx_granges,
            tx_exon_granges_list = transcript_data$tx_exon_granges_list,
            tx_intron_granges_list = transcript_data$tx_intron_granges_list
        )),
        regexp = "transcript_data$tx_df is not a data frame",
        fixed = TRUE)
    expect_error(
        plot_splicing_support_levels(transcript_data = transcript_data_ref_only),
        regexp = "transcript_data contains reference transcripts only",
        fixed = TRUE)
    expect_error(
        plot_splicing_support_levels(transcript_data = list(
            tx_df = transcript_data$tx_df,
            tx_granges = 42,
            tx_exon_granges_list = transcript_data$tx_exon_granges_list,
            tx_intron_granges_list = transcript_data$tx_intron_granges_list
        )),
        regexp = 'class(transcript_data$tx_granges) not equal to "GRanges"',
        fixed = TRUE)
    expect_error(
        plot_splicing_support_levels(transcript_data = list(
            tx_df = transcript_data$tx_df,
            tx_granges = transcript_data$tx_granges,
            tx_exon_granges_list = 42,
            tx_intron_granges_list = transcript_data$tx_intron_granges_list
        )),
        regexp = 'grepl(pattern = "GRangesList", x = class',
        fixed = TRUE)
    expect_error(
        plot_splicing_support_levels(transcript_data = list(
            tx_df = transcript_data$tx_df,
            tx_granges = transcript_data$tx_granges,
            tx_exon_granges_list = transcript_data$tx_exon_granges_list,
            tx_intron_granges_list = 42
        )),
        regexp = 'grepl(pattern = "GRangesList", x = class',
        fixed = TRUE)

    # Testing if function returns the expected output
    expect_silent(
        plot_data <- plot_splicing_support_levels(
            transcript_data = transcript_data
        )
    )
    expect_true(methods::is(plot_data, "ggplot"))
})
