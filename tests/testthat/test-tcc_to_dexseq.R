test_that("tcc_to_dexseq works as expected", {

    # Preparing test data
    se_tcc <- readRDS(system.file("extdata", "se_tcc_mouse_e18.rds",
                                  package = "Isosceles"))
    sce_psi <- readRDS(system.file("extdata", "sce_psi_mouse_e18.rds",
                                   package = "Isosceles"))
    pseudotime_matrix <- readRDS(system.file(
        "extdata", "pseudotime_matrix_mouse_e18.rds", package = "Isosceles"
    ))
    pseudotime <- pseudotime_matrix[, 1]

    # Testing if function throws the expected errors
    expect_error(tcc_to_dexseq(se_tcc = NULL),
                 regexp = "methods::is(object = se_tcc, class2 =",
                 fixed = TRUE)
    se_copy <- se_tcc
    SummarizedExperiment::assay(se_copy, "counts") <- NULL
    expect_error(tcc_to_dexseq(se_tcc = se_copy),
                 regexp = 'is.element(el = "counts",',
                 fixed = TRUE)
    expect_error(tcc_to_dexseq(se_tcc = se_tcc,
                               pseudotime = NULL),
                 regexp = "pseudotime is not a numeric or integer vector",
                 fixed = TRUE)
    expect_error(tcc_to_dexseq(se_tcc = se_tcc,
                               pseudotime = pseudotime[1:10]),
                 regexp = "length(pseudotime) not identical to ncol(se_tcc)",
                 fixed = TRUE)
    expect_error(tcc_to_dexseq(se_tcc = se_tcc,
                               pseudotime = pseudotime,
                               psi_events = 1:10),
                 regexp = "psi_events is not a character vector",
                 fixed = TRUE)
    expect_error(tcc_to_dexseq(se_tcc = se_tcc,
                               pseudotime = pseudotime,
                               psi_events = "jabberwocky"),
                 regexp = "length(psi_events) not greater than 1",
                 fixed = TRUE)
    expect_error(tcc_to_dexseq(se_tcc = se_tcc,
                               pseudotime = pseudotime,
                               trim = NULL),
                 regexp = "trim is not a number",
                 fixed = TRUE)
    expect_error(tcc_to_dexseq(se_tcc = se_tcc,
                               pseudotime = pseudotime,
                               trim = -1),
                 regexp = "trim not greater than or equal to 0",
                 fixed = TRUE)
    expect_error(tcc_to_dexseq(se_tcc = se_tcc,
                               pseudotime = pseudotime,
                               trim = 1),
                 regexp = "trim not less than 0.5",
                 fixed = TRUE)
    expect_error(tcc_to_dexseq(se_tcc = se_tcc,
                               pseudotime = pseudotime,
                               window_size = NULL),
                 regexp = "window_size is not a count",
                 fixed = TRUE)
    expect_error(tcc_to_dexseq(se_tcc = se_tcc,
                               pseudotime = pseudotime,
                               window_step = NULL),
                 regexp = "window_step is not a count",
                 fixed = TRUE)
    expect_error(tcc_to_dexseq(se_tcc = se_tcc,
                               pseudotime = pseudotime,
                               remove_redundant_psi = NULL),
                 regexp = "remove_redundant_psi is not a flag",
                 fixed = TRUE)
    expect_error(tcc_to_dexseq(se_tcc = se_tcc,
                               pseudotime = pseudotime,
                               scale_pseudotime = NULL),
                 regexp = "scale_pseudotime is not a flag",
                 fixed = TRUE)
    expect_error(tcc_to_dexseq(se_tcc = se_tcc,
                               pseudotime = pseudotime,
                               ncpu = NULL),
                 regexp = "ncpu is not a count",
                 fixed = TRUE)

    # Testing if function returns the expected output
    # (default parameters)
    expect_warning(
        dxd <- tcc_to_dexseq(
            se_tcc = se_tcc, pseudotime = pseudotime
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
        c("sample", "window_name", "pseudotime", "exon")
    )
    expect_identical(
        SummarizedExperiment::assayNames(dxd), "counts"
    )
    expect_identical(
        sum(SummarizedExperiment::assay(dxd, "counts")), 324459L
    )

    # Testing if function returns the expected output
    # (scale_pseudotime = FALSE)
    expect_message(
        expect_warning(
            dxd <- tcc_to_dexseq(
                se_tcc = se_tcc, pseudotime = pseudotime,
                scale_pseudotime = FALSE
            ),
            regexp = "some variables",
            fixed = TRUE
        ),
        regexp = "standard deviation larger than 5",
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
        c("sample", "window_name", "pseudotime", "exon")
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
        dxd <- tcc_to_dexseq(
            se_tcc = se_tcc, pseudotime = pseudotime,
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
        c("sample", "window_name", "pseudotime", "exon")
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
        dxd <- tcc_to_dexseq(
            se_tcc = se_tcc, pseudotime = pseudotime,
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
        c("sample", "window_name", "pseudotime", "exon")
    )
    expect_identical(
        SummarizedExperiment::assayNames(dxd), "counts"
    )
    expect_identical(
        sum(SummarizedExperiment::assay(dxd, "counts")), 44758L
    )

    # Testing if function returns the expected output
    # (window size = 20, window step = 10, trim = 0.05)
    expect_warning(
        dxd <- tcc_to_dexseq(
            se_tcc = se_tcc, pseudotime = pseudotime,
            window_size = 20, window_step = 10, trim = 0.05,
        ),
        regexp = "some variables",
        fixed = TRUE
    )
    expect_true(class(dxd) == "DEXSeqDataSet")
    expect_identical(dim(dxd), c(250L, 58L))
    expect_true(all(grepl("^ENSMUSG", rownames(dxd))))
    expect_identical(
        colnames(SummarizedExperiment::rowData(dxd)),
        c("featureID", "groupID", "exonBaseMean", "exonBaseVar")
    )
    expect_identical(
        colnames(SummarizedExperiment::colData(dxd)),
        c("sample", "window_name", "pseudotime", "exon")
    )
    expect_identical(
        SummarizedExperiment::assayNames(dxd), "counts"
    )
    expect_identical(
        sum(SummarizedExperiment::assay(dxd, "counts")), 294011L
    )
})
