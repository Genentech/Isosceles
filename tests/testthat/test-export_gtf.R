test_that("export_gtf works as expected", {

    # Preparing test data
    bam_file <- system.file(
        "extdata", "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam",
        package = "Isosceles"
    )
    bam_files <- c(Sample = bam_file)
    gtf_file <- system.file(
        "extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf",
        package = "Isosceles"
    )
    genome_fasta_file <- system.file(
        "extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa",
        package = "Isosceles"
    )
    bam_parsed <- bam_to_read_structures(bam_files)
    transcript_data <- prepare_transcripts(
        gtf_file = gtf_file, genome_fasta_file = genome_fasta_file,
        bam_parsed = bam_parsed, min_bam_splice_read_count = 2,
        min_bam_splice_fraction = 0.01
    )
    se_tcc <- bam_to_tcc(
        bam_files = bam_files, transcript_data = transcript_data,
        run_mode = "de_novo_loose", min_relative_expression = 0
    )
    se <- tcc_to_transcript(
        se_tcc = se_tcc, use_length_normalization = TRUE
    )
    ref_gtf_file <- system.file(
        "extdata", "transcripts.gtf",
        package = "Isosceles"
    )

    # Testing if function throws the expected errors
    expect_error(export_gtf(se = NULL),
                 regexp = "methods::is(object = se, class2 =",
                 fixed = TRUE)
    se_copy <- se
    SummarizedExperiment::rowRanges(se_copy) <- NULL
    expect_error(export_gtf(se = se_copy),
                 regexp = "SummarizedExperiment::rowRanges(se), class2 =",
                 fixed = TRUE)
    expect_error(export_gtf(se = se,
                            file = NULL),
                 regexp = "file is not a string",
                 fixed = TRUE)

    # Preparing a temporary directory
    temp_dir <- tempfile(pattern = "tmp_")
    unlink(temp_dir, recursive = TRUE)
    dir.create(temp_dir, recursive = TRUE)

    # Testing if function returns the expected output
    expect_silent(
        output <- export_gtf(se, file.path(temp_dir, "transcripts.gtf"))
    )
    expect_true(is.null(output))
    expect_true(file.exists(file.path(temp_dir, "transcripts.gtf")))
    expect_true(tools::md5sum(file.path(temp_dir, "transcripts.gtf"))
                == tools::md5sum(ref_gtf_file))

    # Removing the temporary directory
    unlink(temp_dir, recursive = TRUE)
})
