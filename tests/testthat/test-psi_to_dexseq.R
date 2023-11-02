test_that("psi_to_dexseq works as expected", {

    # Preparing test data
    se_tcc <- readRDS(system.file("extdata", "se_tcc_mouse_e18.rds",
                                  package = "Isosceles"))
    sce_psi <- readRDS(system.file("extdata", "sce_psi_mouse_e18.rds",
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
    se_window_psi <- add_psi_counts(
        se_psi = se_window_psi,
        se_gene = se_window_gene
    )
    window_pseudotime <- scale(se_window_psi$pseudotime)

    # Testing if function throws the expected errors
    expect_error(psi_to_dexseq(se_psi = NULL),
                 regexp = "methods::is(object = se_psi, class2 =",
                 fixed = TRUE)
    se_copy <- se_window_psi
    SummarizedExperiment::assay(se_copy, "counts") <- NULL
    expect_error(psi_to_dexseq(se_psi = se_copy),
                 regexp = 'is.element(el = "counts",',
                 fixed = TRUE)
    se_copy <- se_window_psi
    SummarizedExperiment::assay(se_copy, "other_counts") <- NULL
    expect_error(psi_to_dexseq(se_psi = se_copy),
                 regexp = 'is.element(el = "other_counts",',
                 fixed = TRUE)
    expect_error(psi_to_dexseq(se_psi = se_window_psi,
                               condition = window_pseudotime[1:10]),
                 regexp = "length(condition) not equal to ncol(se_psi)",
                 fixed = TRUE)
    expect_error(psi_to_dexseq(se_psi = se_window_psi,
                               condition = data.frame(
                                   pseudotime = window_pseudotime[1:10]
                               )),
                 regexp = "nrow(condition) not equal to ncol(se_psi)",
                 fixed = TRUE)
    expect_error(psi_to_dexseq(se_psi = se_window_psi,
                               condition = window_pseudotime,
                               design = NULL),
                 regexp = 'methods::is(object = design, class2 = "formula")',
                 fixed = TRUE)
    expect_error(psi_to_dexseq(se_psi = se_window_psi,
                               condition = window_pseudotime,
                               psi_events = 1:10),
                 regexp = "psi_events is not a character vector",
                 fixed = TRUE)
    expect_error(psi_to_dexseq(se_psi = se_window_psi,
                               condition = window_pseudotime,
                               psi_events = "jabberwocky"),
                 regexp = "length(psi_events) not greater than 1",
                 fixed = TRUE)
    expect_error(psi_to_dexseq(se_psi = se_window_psi,
                               condition = window_pseudotime,
                               remove_redundant_psi = NULL),
                 regexp = "remove_redundant_psi is not a flag",
                 fixed = TRUE)

    # Testing if function returns the expected output
    # (default parameters)
    expect_warning(
        dxd <- psi_to_dexseq(
            se_psi = se_window_psi,
            condition = window_pseudotime
        ),
        regexp = "some variables",
        fixed = TRUE
    )
    expect_true(class(dxd) == "DEXSeqDataSet")
    expect_identical(dim(dxd), c(254L, 42L))
    expect_true(all(grepl("^ENSMUSG", rownames(dxd))))
    expect_identical(
        colnames(SummarizedExperiment::rowData(dxd)),
        c("featureID", "groupID", "exonBaseMean", "exonBaseVar")
    )
    expect_identical(
        colnames(SummarizedExperiment::colData(dxd)),
        c("sample", "se_psi_colname", "condition", "exon")
    )
    expect_identical(
        SummarizedExperiment::assayNames(dxd), "counts"
    )
    expect_identical(
        sum(SummarizedExperiment::assay(dxd, "counts")), 324459L
    )

    # Testing if function returns the expected output
    # (condition argument is a data frame, custom design formula)
    expect_warning(
        dxd <- psi_to_dexseq(
            se_psi = se_window_psi,
            condition = data.frame(pseudotime = window_pseudotime),
            design = ~ sample + exon + pseudotime:exon
        ),
        regexp = "some variables",
        fixed = TRUE
    )
    expect_true(class(dxd) == "DEXSeqDataSet")
    expect_identical(dim(dxd), c(254L, 42L))
    expect_true(all(grepl("^ENSMUSG", rownames(dxd))))
    expect_identical(
        colnames(SummarizedExperiment::rowData(dxd)),
        c("featureID", "groupID", "exonBaseMean", "exonBaseVar")
    )
    expect_identical(
        colnames(SummarizedExperiment::colData(dxd)),
        c("sample", "se_psi_colname", "pseudotime", "exon")
    )
    expect_identical(
        SummarizedExperiment::assayNames(dxd), "counts"
    )
    expect_identical(
        sum(SummarizedExperiment::assay(dxd, "counts")), 324459L
    )

    # Testing if function returns the expected output
    # (remove_redundant_psi = FALSE)
    expect_warning(
        dxd <- psi_to_dexseq(
            se_psi = se_window_psi,
            condition = window_pseudotime,
            remove_redundant_psi = FALSE
        ),
        regexp = "some variables",
        fixed = TRUE
    )
    expect_true(class(dxd) == "DEXSeqDataSet")
    expect_identical(dim(dxd), c(410L, 42L))
    expect_true(all(grepl("^ENSMUSG", rownames(dxd))))
    expect_identical(
        colnames(SummarizedExperiment::rowData(dxd)),
        c("featureID", "groupID", "exonBaseMean", "exonBaseVar")
    )
    expect_identical(
        colnames(SummarizedExperiment::colData(dxd)),
        c("sample", "se_psi_colname", "condition", "exon")
    )
    expect_identical(
        SummarizedExperiment::assayNames(dxd), "counts"
    )
    expect_identical(
        sum(SummarizedExperiment::assay(dxd, "counts")), 473744L
    )

    # Testing if function returns the expected output
    # (psi_events = rownames(sce_psi))
    expect_warning(
        dxd <- psi_to_dexseq(
            se_psi = se_window_psi,
            condition = window_pseudotime,
            psi_events = rownames(sce_psi)
        ),
        regexp = "some variables",
        fixed = TRUE
    )
    expect_true(class(dxd) == "DEXSeqDataSet")
    expect_identical(dim(dxd), c(30L, 42L))
    expect_true(all(grepl("^ENSMUSG", rownames(dxd))))
    expect_identical(
        colnames(SummarizedExperiment::rowData(dxd)),
        c("featureID", "groupID", "exonBaseMean", "exonBaseVar")
    )
    expect_identical(
        colnames(SummarizedExperiment::colData(dxd)),
        c("sample", "se_psi_colname", "condition", "exon")
    )
    expect_identical(
        SummarizedExperiment::assayNames(dxd), "counts"
    )
    expect_identical(
        sum(SummarizedExperiment::assay(dxd, "counts")), 44758L
    )
})
