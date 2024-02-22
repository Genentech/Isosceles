test_that("gtf_to_intron_bed works as expected", {

    # Preparing a temporary directory
    temp_dir <- tempfile(pattern = "tmp_")
    unlink(temp_dir, recursive = TRUE)
    dir.create(temp_dir, recursive = TRUE)

    # Preparing test data
    gtf_file <- system.file(
        "extdata", "bulk_rnaseq.gtf",
        package = "Isosceles"
    )
    ref_bed_file <- system.file(
        "extdata", "introns.bed",
        package = "Isosceles"
    )

    # Testing if function throws the expected errors
    expect_error(gtf_to_intron_bed(gtf_file = NULL),
                 regexp = "gtf_file is not a string",
                 fixed = TRUE)
    expect_error(gtf_to_intron_bed(gtf_file = "jabberwocky"),
                 regexp = "Path 'jabberwocky' does not exist",
                 fixed = TRUE)
    expect_error(gtf_to_intron_bed(gtf_file = gtf_file,
                                   file = NULL),
                 regexp = "file is not a string",
                 fixed = TRUE)

    # Testing if function returns the expected output
    expect_message(
        output <- gtf_to_intron_bed(
            gtf_file = gtf_file,
            file = file.path(temp_dir, "introns.bed")
        ),
        regexp = "TxDb",
        fixed = TRUE
    )
    expect_true(is.null(output))
    expect_true(file.exists(file.path(temp_dir, "introns.bed")))
    expect_true(tools::md5sum(file.path(temp_dir, "introns.bed"))
                == tools::md5sum(ref_bed_file))

    # Preparing test data (GTF file with wrong transcript strand values)
    gtf_file <- system.file(
        "extdata", "wrong_strand.gtf",
        package = "Isosceles"
    )
    ref_bed_file <- system.file(
        "extdata", "introns_wrong_strand.bed",
        package = "Isosceles"
    )

    # Testing if function returns the expected output (GTF file with wrong
    # transcript strand values)
    expect_warning(
        output <- gtf_to_intron_bed(
            gtf_file = gtf_file,
            file = file.path(temp_dir, "introns_wrong_strand.bed")
        ),
        regexp = "different strand than their exons: Transcript_2, Transcript_3",
        fixed = TRUE
    )
    expect_true(is.null(output))
    expect_true(file.exists(file.path(temp_dir, "introns_wrong_strand.bed")))
    expect_true(tools::md5sum(file.path(temp_dir, "introns_wrong_strand.bed"))
                == tools::md5sum(ref_bed_file))

    # Removing the temporary directory
    unlink(temp_dir, recursive = TRUE)
})
