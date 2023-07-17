test_that("prepare_reference_annotations works as expected", {

    # Preparing test data
    gtf_file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf",
                            package = "Isosceles")

    # Testing if function throws the expected errors
    expect_error(prepare_reference_annotations(gtf_file = NULL),
                 regexp = "gtf_file is not a string",
                 fixed = TRUE)
    expect_error(prepare_reference_annotations(gtf_file = "jabberwocky"),
                 regexp = "Path 'jabberwocky' does not exist",
                 fixed = TRUE)
    expect_error(prepare_reference_annotations(gtf_file = gtf_file,
                                               is_technical = NULL),
                 regexp = "is_technical is not a flag",
                 fixed = TRUE)

    # Testing if function returns the expected output
    expect_message(
        anno_data <- prepare_reference_annotations(
            gtf_file = gtf_file, is_technical = FALSE
        ),
        regexp = "gene_id",
        fixed = TRUE
    )
    expect_true(is.list(anno_data))
    expect_identical(names(anno_data),
                     c("gene_df", "transcript_df", "splicing_df",
                       "intron_df", "transcript_first_last_df"))
    expect_true(all(sapply(anno_data, is.data.frame)))
    expect_identical(dim(anno_data$gene_df), c(23L, 2L))
    expect_identical(dim(anno_data$transcript_df), c(105L, 6L))
    expect_identical(dim(anno_data$splicing_df), c(96L, 5L))
    expect_identical(dim(anno_data$intron_df), c(726L, 3L))
    expect_identical(dim(anno_data$transcript_first_last_df), c(37L, 6L))
})
