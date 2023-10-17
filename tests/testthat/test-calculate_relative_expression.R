test_that("calculate_relative_expression works as expected", {

    # Preparing test data (bulk RNA-Seq data)
    bam_file <- system.file(
        "extdata", "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam",
        package = "Isosceles"
    )
    bam_files <- c(Sample_1 = bam_file, Sample_2 = bam_file)
    gtf_file <- system.file(
        "extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf",
        package = "Isosceles"
    )
    genome_fasta_file <- system.file(
        "extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa",
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
    expect_error(calculate_relative_expression(se = NULL),
                 regexp = "methods::is(object = se, class2 =",
                 fixed = TRUE)
    se_copy <- se
    SummarizedExperiment::assay(se_copy, "tpm") <- NULL
    expect_error(calculate_relative_expression(se = se_copy),
                 regexp = 'is.element(el = "tpm",',
                 fixed = TRUE)
    se_copy <- se
    SummarizedExperiment::rowData(se_copy)$gene_id <- NULL
    expect_error(calculate_relative_expression(se = se_copy),
                 regexp = 'is.element(el = "gene_id",',
                 fixed = TRUE)

    # Testing if function returns the expected output (bulk RNA-Seq data)
    expect_silent(
        rel_expr_matrix <- calculate_relative_expression(se = se)
    )
    expect_true(is.matrix(rel_expr_matrix))
    expect_identical(dim(rel_expr_matrix), dim(se))
    expect_identical(colnames(rel_expr_matrix), colnames(se))
    expect_identical(rownames(rel_expr_matrix), rownames(se))
    expect_identical(round(colSums(rel_expr_matrix)),
                     c(Sample_1 = 4, Sample_2 = 4))
    expect_identical(rel_expr_matrix[, 1], rel_expr_matrix[, 2])

    # Preparing test data (scRNA-Seq data)
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
        rel_expr_matrix <- calculate_relative_expression(se = sce)
    )
    expect_true(class(rel_expr_matrix) == "dgCMatrix")
    expect_identical(dim(rel_expr_matrix), dim(sce))
    expect_identical(colnames(rel_expr_matrix), colnames(sce))
    expect_identical(rownames(rel_expr_matrix), rownames(sce))
    expect_true(all(rel_expr_matrix >= 0))
    expect_true(all(rel_expr_matrix <= 1))
    expect_true(all(round(Matrix::colSums(rel_expr_matrix)) == 1))
})
