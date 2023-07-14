test_that("prepare_pseudobulk_se works as expected", {

    # Preparing test data
    bam_file <- system.file(
        "extdata", "molecule.tags.GE.bam",
        package = "Isosceles"
    )
    bam_files <- c(Sample = bam_file)
    gtf_file <- system.file(
        "extdata", "chr4.gtf.gz",
        package = "Isosceles"
    )
    genome_fasta_file <- system.file(
        "extdata", "chr4.fa.gz",
        package = "Isosceles"
    )
    bam_parsed <- extract_read_structures(bam_files)
    transcript_data <- prepare_transcripts(
        gtf_file = gtf_file, genome_fasta_file = genome_fasta_file,
        bam_parsed = bam_parsed, min_bam_splice_read_count = 1,
        min_bam_splice_fraction = 0.01
    )
    se_tcc <- prepare_tcc_se(
        bam_files = bam_files, transcript_data = transcript_data,
        is_single_cell = TRUE, barcode_tag = "BC",
        run_mode = "de_novo_loose", min_relative_expression = 0
    )
    set.seed(42)
    cell_labels <- sample(1:2, ncol(se_tcc), replace = TRUE)

    # Testing if function throws the expected errors
    expect_error(prepare_pseudobulk_se(se_tcc = NULL),
                 regexp = "x = class(se_tcc)) is not TRUE",
                 fixed = TRUE)
    se_copy <- se_tcc
    SummarizedExperiment::assay(se_copy, "counts") <- NULL
    expect_error(prepare_pseudobulk_se(se_tcc = se_copy),
                 regexp = 'is.element(el = "counts",',
                 fixed = TRUE)
    expect_error(prepare_pseudobulk_se(se_tcc = se_tcc,
                                       cell_labels = NULL),
                 regexp = "length(cell_labels) not equal to ncol(se_tcc)",
                 fixed = TRUE)

    # Testing if function returns the expected output
    expect_silent(
        se <- prepare_pseudobulk_se(
            se_tcc = se_tcc, cell_labels = cell_labels
        )
    )
    expect_true(class(se) == "SummarizedExperiment")
    expect_identical(dim(se), c(nrow(se_tcc), length(unique(cell_labels))))
    expect_identical(colnames(se), levels(factor(cell_labels)))
    expect_identical(rownames(se), rownames(se_tcc))
    expect_identical(SummarizedExperiment::rowRanges(se),
                     SummarizedExperiment::rowRanges(se_tcc))
    expect_identical(SummarizedExperiment::rowData(se),
                     SummarizedExperiment::rowData(se_tcc))
    expect_identical(dim(SummarizedExperiment::colData(se)),
                     c(length(unique(cell_labels)), 0L))
    expect_identical(SummarizedExperiment::assayNames(se),
                     c("counts", "tpm", "relative_expression"))
    expect_true(is.matrix(SummarizedExperiment::assay(se, "counts")))
    expect_true(is.matrix(SummarizedExperiment::assay(se, "tpm")))
    expect_true(is.matrix(SummarizedExperiment::assay(se, "relative_expression")))
    expect_identical(round(Matrix::rowSums(SummarizedExperiment::assay(se, "counts"))),
                     round(Matrix::rowSums(SummarizedExperiment::assay(se_tcc, "counts"))))
    expect_true(all(round(Matrix::colSums(SummarizedExperiment::assay(se, "tpm"))) == 1e6))
    expect_true(all(round(Matrix::colSums(SummarizedExperiment::assay(se, "relative_expression"))) == 1))
    expect_identical(S4Vectors::metadata(se),
                     S4Vectors::metadata(se_tcc))
})
