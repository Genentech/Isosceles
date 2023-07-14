test_that("extract_read_structures works as expected", {

    # Preparing test data
    bam_file <- system.file(
        "extdata", "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam",
        package = "Isosceles"
    )

    # Testing if function throws the expected errors
    expect_error(extract_read_structures(bam_files = NULL),
                 regexp = "bam_files is not a character vector",
                 fixed = TRUE)
    expect_error(extract_read_structures(bam_files = character(0)),
                 regexp = "length(bam_files) not greater than 0",
                 fixed = TRUE)
    expect_error(extract_read_structures(bam_files = "jabberwocky"),
                 regexp = "Elements 1 of file.exists(bam_files) are not true",
                 fixed = TRUE)
    expect_error(extract_read_structures(bam_files = bam_file,
                                         chunk_size = NULL),
                 regexp = "chunk_size is not a count",
                 fixed = TRUE)
    expect_error(extract_read_structures(bam_files = bam_file,
                                         ncpu = NULL),
                 regexp = "ncpu is not a count",
                 fixed = TRUE)


    # Testing if function returns the expected output
    expect_message(bam_data <- extract_read_structures(bam_files = bam_file),
                   regexp = 'Joining, by = "read_id"',
                   fixed = TRUE)
    expect_true(is.data.frame(bam_data))
    expect_identical(dim(bam_data), c(144L, 5L))
    expect_identical(colnames(bam_data),
                     c("intron_positions", "read_count", "chromosome",
                       "start_positions", "end_positions"))
})
