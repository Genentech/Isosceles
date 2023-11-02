test_that("add_psi_counts works as expected", {

    # Preparing test data
    se_tcc <- readRDS(system.file("extdata", "se_tcc_mouse_e18.rds",
                                  package = "Isosceles"))
    pseudotime_matrix <- readRDS(system.file(
        "extdata", "pseudotime_matrix_mouse_e18.rds", package = "Isosceles"
    ))
    pseudotime <- pseudotime_matrix[, 1]
    se_window_tcc <- pseudotime_tcc(
        se_tcc = se_tcc,
        pseudotime = pseudotime,
        window_size = 30,
        window_step = 15
    )
    se_window_gene <- tcc_to_gene(
        se_tcc = se_window_tcc
    )
    se_window_transcript <- tcc_to_transcript(
        se_tcc = se_window_tcc,
        use_length_normalization = FALSE
    )
    se_window_psi <- transcript_to_psi(
        se = se_window_transcript,
    )

    # Testing if function throws the expected errors
    expect_error(add_psi_counts(se_psi = NULL),
                 regexp = "methods::is(object = se_psi, class2 =",
                 fixed = TRUE)
    se_copy <- se_window_psi
    SummarizedExperiment::assay(se_copy, "psi") <- NULL
    expect_error(add_psi_counts(se_psi = se_copy),
                 regexp = 'is.element(el = "psi",',
                 fixed = TRUE)
    se_copy <- se_window_psi
    SummarizedExperiment::rowData(se_copy)$gene_id <- NULL
    expect_error(add_psi_counts(se_psi = se_copy),
                 regexp = 'is.element(el = "gene_id", set = colnames',
                 fixed = TRUE)
    expect_error(add_psi_counts(se_psi = se_window_psi,
                                se_gene = NULL),
                 regexp = "methods::is(object = se_gene, class2 =",
                 fixed = TRUE)
    se_copy <- se_window_gene
    SummarizedExperiment::assay(se_copy, "counts") <- NULL
    expect_error(add_psi_counts(se_psi = se_window_psi,
                                se_gene = se_copy),
                 regexp = 'is.element(el = "counts",',
                 fixed = TRUE)
    se_copy <- se_window_gene
    SummarizedExperiment::rowData(se_copy)$gene_id <- NULL
    expect_error(add_psi_counts(se_psi = se_window_psi,
                                se_gene = se_copy),
                 regexp = 'is.element(el = "gene_id", set = colnames',
                 fixed = TRUE)
    expect_error(add_psi_counts(se_psi = se_window_psi,
                                se_gene = se_window_gene[, 1:10]),
                 regexp = "ncol(se_psi) not identical to ncol(se_gene)",
                 fixed = TRUE)

    # Testing if function returns the expected output
    expect_silent(
        se <- add_psi_counts(
            se_psi = se_window_psi,
            se_gene = se_window_gene
        )
    )
    assertthat::assert_that(methods::is(se, "SummarizedExperiment"))
    se_copy <- se
    SummarizedExperiment::assay(se_copy, "counts") <- NULL
    SummarizedExperiment::assay(se_copy, "other_counts") <- NULL
    expect_identical(se_copy, se_window_psi)
    expect_true(class(SummarizedExperiment::assay(se, "counts")) == "dgCMatrix")
    expect_true(class(SummarizedExperiment::assay(se, "other_counts")) == "dgCMatrix")
    expect_identical(
        round(sum(Matrix::colSums(SummarizedExperiment::assay(se, "counts")))),
        219332
    )
    expect_identical(
        round(sum(Matrix::colSums(SummarizedExperiment::assay(se, "other_counts")))),
        254412
    )
})
