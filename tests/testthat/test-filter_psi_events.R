test_that("filter_psi_events works as expected", {

    # Preparing test data
    sce_psi <- readRDS(system.file(
        "extdata", "sce_psi_mouse_e18.rds",
        package = "Isosceles"
    ))

    # Testing if function throws the expected errors
    expect_error(filter_psi_events(se_psi = NULL),
                 regexp = "methods::is(object = se_psi, class2 =",
                 fixed = TRUE)
    sce_psi_copy <- sce_psi
    SummarizedExperiment::assay(sce_psi_copy, "psi") <- NULL
    expect_error(filter_psi_events(se_psi = sce_psi_copy),
                 regexp = 'is.element(el = "psi",',
                 fixed = TRUE)
    sce_psi_copy <- sce_psi
    SummarizedExperiment::rowData(sce_psi_copy)$gene_id <- NULL
    expect_error(filter_psi_events(se_psi = sce_psi_copy),
                 regexp = 'is.element(el = "gene_id", set = colnames',
                 fixed = TRUE)
    expect_error(filter_psi_events(se_psi = sce_psi,
                                   gene_ids = 1:10),
                 regexp = "gene_ids is not a character vector",
                 fixed = TRUE)
    expect_error(filter_psi_events(se_psi = sce_psi,
                                   exclude_tss = NULL),
                 regexp = "exclude_tss is not a flag",
                 fixed = TRUE)
    expect_error(filter_psi_events(se_psi = sce_psi,
                                   exclude_tes = NULL),
                 regexp = "exclude_tes is not a flag",
                 fixed = TRUE)
    expect_error(filter_psi_events(se_psi = sce_psi,
                                   min_mean_psi = NULL),
                 regexp = "min_mean_psi is not a number",
                 fixed = TRUE)
    expect_error(filter_psi_events(se_psi = sce_psi,
                                   min_mean_psi = -1),
                 regexp = "min_mean_psi not greater than or equal to 0",
                 fixed = TRUE)
    expect_error(filter_psi_events(se_psi = sce_psi,
                                   min_mean_psi = 2),
                 regexp = "min_mean_psi not less than or equal to 1",
                 fixed = TRUE)
    expect_error(filter_psi_events(se_psi = sce_psi,
                                   max_mean_psi = NULL),
                 regexp = "max_mean_psi is not a number",
                 fixed = TRUE)
    expect_error(filter_psi_events(se_psi = sce_psi,
                                   max_mean_psi = -1),
                 regexp = "max_mean_psi not greater than or equal to 0",
                 fixed = TRUE)
    expect_error(filter_psi_events(se_psi = sce_psi,
                                   max_mean_psi = 2),
                 regexp = "max_mean_psi not less than or equal to 1",
                 fixed = TRUE)
    expect_error(filter_psi_events(se_psi = sce_psi,
                                   min_alt_psi_count = NULL),
                 regexp = "min_alt_psi_count is not a count",
                 fixed = TRUE)
    expect_error(filter_psi_events(se_psi = sce_psi,
                                   min_incl_psi = NULL),
                 regexp = "min_incl_psi is not a number",
                 fixed = TRUE)
    expect_error(filter_psi_events(se_psi = sce_psi,
                                   min_incl_psi = -1),
                 regexp = "min_incl_psi not greater than or equal to 0",
                 fixed = TRUE)
    expect_error(filter_psi_events(se_psi = sce_psi,
                                   min_incl_psi = 2),
                 regexp = "min_incl_psi not less than or equal to 1",
                 fixed = TRUE)
    expect_error(filter_psi_events(se_psi = sce_psi,
                                   min_incl_psi_count = NULL),
                 regexp = "min_incl_psi_count is not a count",
                 fixed = TRUE)

    # Testing if function returns the expected output
    expect_silent(
        psi_events <- filter_psi_events(sce_psi)
    )
    expect_true(is.character(psi_events))
    expect_identical(length(psi_events), 30L)
    expect_identical(length(filter_psi_events(sce_psi, gene_ids = "ENSMUSG00000002107")),
                     6L)
    sce_psi_copy <- sce_psi
    rownames(sce_psi_copy)[1] <- "ENSMUSG00000002107:chr2:6664115:-:TSS"
    expect_identical(length(filter_psi_events(sce_psi_copy, exclude_tss = FALSE)),
                     30L)
    expect_identical(length(filter_psi_events(sce_psi_copy, exclude_tss = TRUE)),
                     29L)
    sce_psi_copy <- sce_psi
    rownames(sce_psi_copy)[1] <- "ENSMUSG00000002107:chr2:6664115:-:TES"
    expect_identical(length(filter_psi_events(sce_psi_copy, exclude_tes = FALSE)),
                     30L)
    expect_identical(length(filter_psi_events(sce_psi_copy, exclude_tes = TRUE)),
                     29L)
    expect_identical(length(filter_psi_events(sce_psi, min_mean_psi = 0.1)),
                     23L)
    expect_identical(length(filter_psi_events(sce_psi, max_mean_psi = 0.8)),
                     29L)
    expect_identical(length(filter_psi_events(sce_psi, min_alt_psi_count = 100)),
                     26L)
    expect_identical(length(filter_psi_events(sce_psi, min_incl_psi = 0.3)),
                     24L)
    expect_identical(length(filter_psi_events(sce_psi, min_incl_psi_count = 100)),
                     29L)
})
