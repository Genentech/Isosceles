test_that("calculate_tpm works as expected", {

    # Preparing test data (bulk RNA-Seq data)
    bam_file <- system.file(
        "extdata", "bulk_rnaseq.bam",
        package = "Isosceles"
    )
    bam_files <- c(Sample_1 = bam_file, Sample_2 = bam_file)
    gtf_file <- system.file(
        "extdata", "bulk_rnaseq.gtf",
        package = "Isosceles"
    )
    genome_fasta_file <- system.file(
        "extdata", "bulk_rnaseq.fa",
        package = "Isosceles"
    )
    transcript_data <- prepare_transcripts(
        gtf_file = gtf_file, genome_fasta_file = genome_fasta_file,
        bam_parsed = NULL
    )
    se <- bam_to_tcc(
        bam_files = bam_files, transcript_data = transcript_data
    )

    # Testing if function throws the expected errors
    expect_error(calculate_tpm(se = NULL),
                 regexp = "methods::is(object = se, class2 =",
                 fixed = TRUE)
    se_copy <- se
    SummarizedExperiment::assay(se_copy, "counts") <- NULL
    expect_error(calculate_tpm(se = se_copy),
                 regexp = 'is.element(el = "counts",',
                 fixed = TRUE)

    # Testing if function returns the expected output (bulk RNA-Seq data)
    expect_silent(
        tpm_matrix <- calculate_tpm(se = se)
    )
    expect_true(is.matrix(tpm_matrix))
    expect_identical(dim(tpm_matrix), dim(se))
    expect_identical(colnames(tpm_matrix), colnames(se))
    expect_identical(rownames(tpm_matrix), rownames(se))
    expect_true(all(round(colSums(tpm_matrix)) == 1e6))
    expect_identical(tpm_matrix[, 1], tpm_matrix[, 2])

    # Preparing test data (scRNA-Seq data)
    bam_file <- system.file(
        "extdata", "scrnaseq.bam",
        package = "Isosceles"
    )
    bam_files <- c(Sample = bam_file)
    gtf_file <- system.file(
        "extdata", "scrnaseq.gtf.gz",
        package = "Isosceles"
    )
    genome_fasta_file <- system.file(
        "extdata", "scrnaseq.fa.gz",
        package = "Isosceles"
    )
    transcript_data <- prepare_transcripts(
        gtf_file = gtf_file, genome_fasta_file = genome_fasta_file,
        bam_parsed = NULL
    )
    se <- bam_to_tcc(
        bam_files = bam_files, transcript_data = transcript_data,
        is_single_cell = TRUE, barcode_tag = "BC"
    )
    sce <- methods::as(se, "SingleCellExperiment")

    # Testing if function returns the expected output (scRNA-Seq data)
    expect_silent(
        tpm_matrix <- calculate_tpm(se = sce)
    )
    expect_true(class(tpm_matrix) == "dgCMatrix")
    expect_identical(dim(tpm_matrix), dim(sce))
    expect_identical(colnames(tpm_matrix), colnames(sce))
    expect_identical(rownames(tpm_matrix), rownames(sce))
    expect_true(all(round(Matrix::colSums(tpm_matrix)) == 1e6))
})
