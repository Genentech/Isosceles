test_that("calculate_psi_ratio_matrix works as expected", {

    # Preparing test data
    se_tcc <- readRDS(system.file("extdata", "se_tcc_mouse_e18.rds",
                                  package = "Isosceles"))
    pseudotime_matrix <- readRDS(system.file(
        "extdata", "pseudotime_matrix_mouse_e18.rds", package = "Isosceles"
    ))
    psi_events <- c(
        "ENSMUSG00000002107:chr2:6664115-6664197:-:CE",
        "ENSMUSG00000002107:chr2:6560659-6560670:-:A5",
        "ENSMUSG00000002107:chr2:6553965-6553982:-:A3",
        "ENSMUSG00000002107:chr2:6549832-6549975:-:CE",
        "ENSMUSG00000002107:chr2:6546780-6547041:-:RI",
        "ENSMUSG00000002107:chr2:6545676-6546779:-:A3"
    )
    window_sizes <- c(
        glut_1 = 100, glut_2 = 100, gaba = 100,
        rad_glia = 50, cyc_rad_glia = 50, cr = 50
    )
    window_steps <- c(
        glut_1 = 3, glut_2 = 3, gaba = 3,
        rad_glia = 3, cyc_rad_glia = 3, cr = 3
    )

    # Testing if function throws the expected errors
    expect_error(calculate_psi_ratio_matrix(se_tcc = NULL),
                 regexp = "methods::is(object = se_tcc, class2 =",
                 fixed = TRUE)
    se_copy <- se_tcc
    SummarizedExperiment::assay(se_copy, "counts") <- NULL
    expect_error(calculate_psi_ratio_matrix(se_tcc = se_copy),
                 regexp = 'is.element(el = "counts",',
                 fixed = TRUE)
    se_copy <- se_tcc
    SummarizedExperiment::rowData(se_copy)$gene_id <- NULL
    expect_error(calculate_psi_ratio_matrix(se_tcc = se_copy),
                 regexp = 'is.element(el = "gene_id",',
                 fixed = TRUE)
    expect_error(calculate_psi_ratio_matrix(se_tcc = se_tcc,
                                            pseudotime_matrix = NULL),
                 regexp = "pseudotime_matrix is not a matrix",
                 fixed = TRUE)
    expect_error(calculate_psi_ratio_matrix(se_tcc = se_tcc,
                                            pseudotime_matrix = pseudotime_matrix[1:10,]),
                 regexp = "nrow(pseudotime_matrix) not identical to ncol(se_tcc)",
                 fixed = TRUE)
    expect_error(calculate_psi_ratio_matrix(se_tcc = se_tcc,
                                            pseudotime_matrix = pseudotime_matrix,
                                            psi_events = NULL),
                 regexp = "psi_events is not a character vector",
                 fixed = TRUE)
    expect_error(calculate_psi_ratio_matrix(se_tcc = se_tcc,
                                            pseudotime_matrix = pseudotime_matrix,
                                            psi_events = "jabberwocky"),
                 regexp = "length(psi_events) not greater than 1",
                 fixed = TRUE)
    expect_error(calculate_psi_ratio_matrix(se_tcc = se_tcc,
                                            pseudotime_matrix = pseudotime_matrix,
                                            psi_events = psi_events,
                                            window_sizes = NULL),
                 regexp = "window_sizes is not a numeric or integer vector",
                 fixed = TRUE)
    expect_error(calculate_psi_ratio_matrix(se_tcc = se_tcc,
                                            pseudotime_matrix = pseudotime_matrix,
                                            psi_events = psi_events,
                                            window_sizes = window_sizes[1:5]),
                 regexp = "sort(names(window_sizes)) not identical to",
                 fixed = TRUE)
    expect_error(calculate_psi_ratio_matrix(se_tcc = se_tcc,
                                            pseudotime_matrix = pseudotime_matrix,
                                            psi_events = psi_events,
                                            window_sizes = window_sizes,
                                            window_steps = NULL),
                 regexp = "window_steps is not a numeric or integer vector",
                 fixed = TRUE)
    expect_error(calculate_psi_ratio_matrix(se_tcc = se_tcc,
                                            pseudotime_matrix = pseudotime_matrix,
                                            psi_events = psi_events,
                                            window_sizes = window_sizes,
                                            window_steps = window_steps[1:5]),
                 regexp = "sort(names(window_steps)) not identical to",
                 fixed = TRUE)
    expect_error(calculate_psi_ratio_matrix(se_tcc = se_tcc,
                                            pseudotime_matrix = pseudotime_matrix,
                                            psi_events = psi_events,
                                            window_sizes = window_sizes,
                                            window_steps = window_steps,
                                            trim = NULL),
                 regexp = "trim is not a number",
                 fixed = TRUE)
    expect_error(calculate_psi_ratio_matrix(se_tcc = se_tcc,
                                            pseudotime_matrix = pseudotime_matrix,
                                            psi_events = psi_events,
                                            window_sizes = window_sizes,
                                            window_steps = window_steps,
                                            trim = -1),
                 regexp = "trim not greater than or equal to 0",
                 fixed = TRUE)
    expect_error(calculate_psi_ratio_matrix(se_tcc = se_tcc,
                                            pseudotime_matrix = pseudotime_matrix,
                                            psi_events = psi_events,
                                            window_sizes = window_sizes,
                                            window_steps = window_steps,
                                            trim = 1),
                 regexp = "trim not less than 0.5",
                 fixed = TRUE)
    expect_error(calculate_psi_ratio_matrix(se_tcc = se_tcc,
                                            pseudotime_matrix = pseudotime_matrix,
                                            psi_events = psi_events,
                                            window_sizes = window_sizes,
                                            window_steps = window_steps,
                                            n_perm = NULL),
                 regexp = "n_perm is not a count",
                 fixed = TRUE)
    expect_error(calculate_psi_ratio_matrix(se_tcc = se_tcc,
                                            pseudotime_matrix = pseudotime_matrix,
                                            psi_events = psi_events,
                                            window_sizes = window_sizes,
                                            window_steps = window_steps,
                                            ncpu = NULL),
                 regexp = "ncpu is not a count",
                 fixed = TRUE)

    # Testing if function returns the expected output
    set.seed(42)
    expect_silent(
        psi_mat <- calculate_psi_ratio_matrix(
            se_tcc = se_tcc,
            pseudotime_matrix = pseudotime_matrix,
            psi_events = psi_events,
            window_sizes = window_sizes,
            window_steps = window_steps,
            n_perm = 1)
    )
    expect_true(is.matrix(psi_mat))
    expect_true(is.numeric(psi_mat))
    expect_identical(dim(psi_mat), c(6L, 191L))
    expect_identical(rownames(psi_mat), psi_events)
    expect_identical(
        unique(sapply(strsplit(colnames(psi_mat), "\\."), "[", 1)),
        colnames(pseudotime_matrix)
    )
})
