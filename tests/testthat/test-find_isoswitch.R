test_that("find_isoswitch works as expected", {

    # Preparing test data
    sce_transcript <- readRDS(system.file("extdata", "sce_transcript_mouse_e18.rds",
                                          package = "Isosceles"))

    # Testing if function throws the expected errors
    expect_error(find_isoswitch(se = NULL),
                 regexp = "methods::is(object = se, class2 =",
                 fixed = TRUE)
    sce_copy <- sce_transcript
    SummarizedExperiment::assay(sce_copy, "logcounts") <- NULL
    expect_error(find_isoswitch(se = sce_copy),
                 regexp = 'is.element(el = "logcounts",',
                 fixed = TRUE)
    sce_copy <- sce_transcript
    SummarizedExperiment::rowRanges(sce_copy) <- NULL
    expect_error(find_isoswitch(se = sce_copy),
                 regexp = "rowRanges(se))) is not TRUE",
                 fixed = TRUE)
    expect_error(find_isoswitch(se = sce_transcript,
                                cell_labels = NULL),
                 regexp = "length(cell_labels) not equal to ncol(se)",
                 fixed = TRUE)
    expect_error(find_isoswitch(se = sce_transcript,
                                cell_labels = sce_transcript$cluster,
                                min_fdr = NULL),
                 regexp = "min_fdr is not a number",
                 fixed = TRUE)
    expect_error(find_isoswitch(se = sce_transcript,
                                cell_labels = sce_transcript$cluster,
                                min_fdr = 0),
                 regexp = "min_fdr not greater than 0",
                 fixed = TRUE)
    expect_error(find_isoswitch(se = sce_transcript,
                                cell_labels = sce_transcript$cluster,
                                min_fdr = 1),
                 regexp = "min_fdr not less than 1",
                 fixed = TRUE)
    expect_error(find_isoswitch(se = sce_transcript,
                                cell_labels = sce_transcript$cluster,
                                ncpu = NULL),
                 regexp = "ncpu is not a count",
                 fixed = TRUE)

    # Testing if function returns the expected output (min_fdr = 0.05)
    expect_silent(
        isoswitch_df <- find_isoswitch(
            se = sce_transcript, cell_labels = sce_transcript$cluster,
            min_fdr = 0.05
        )
    )
    expect_true(is.data.frame(isoswitch_df))
    expect_identical(dim(isoswitch_df), c(95L, 10L))
    expect_identical(colnames(isoswitch_df),
                     c("transcript_id", "compatible_tx", "gene_id",
                       "gene_name", "pvalue", "fdr", "auc", "group_1",
                       "group_2", "contrast"))

    # Testing if function returns the expected output (min_fdr = 0.01)
    expect_silent(
        isoswitch_df <- find_isoswitch(
            se = sce_transcript, cell_labels = sce_transcript$cluster,
            min_fdr = 0.01
        )
    )
    expect_true(is.data.frame(isoswitch_df))
    expect_identical(dim(isoswitch_df), c(35L, 10L))
    expect_identical(colnames(isoswitch_df),
                     c("transcript_id", "compatible_tx", "gene_id",
                       "gene_name", "pvalue", "fdr", "auc", "group_1",
                       "group_2", "contrast"))
})
