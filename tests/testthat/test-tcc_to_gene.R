test_that("tcc_to_gene works as expected", {

    # Preparing test data (bulk RNA-Seq data)
    bam_file <- system.file(
        "extdata", "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam",
        package = "Isosceles"
    )
    bam_files <- c(Sample = bam_file)
    gtf_file <- system.file(
        "extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf",
        package = "Isosceles"
    )
    genome_fasta_file <- system.file(
        "extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa",
        package = "Isosceles"
    )
    bam_parsed <- bam_to_read_structures(bam_files)
    transcript_data <- prepare_transcripts(
        gtf_file = gtf_file, genome_fasta_file = genome_fasta_file,
        bam_parsed = bam_parsed, min_bam_splice_read_count = 2,
        min_bam_splice_fraction = 0.01
    )
    se_tcc <- bam_to_tcc(
        bam_files = bam_files, transcript_data = transcript_data,
        run_mode = "de_novo_loose", min_relative_expression = 0
    )

    # Testing if function throws the expected errors
    expect_error(tcc_to_gene(se_tcc = NULL),
                 regexp = "methods::is(object = se_tcc, class2 =",
                 fixed = TRUE)
    se_copy <- se_tcc
    SummarizedExperiment::assay(se_copy, "counts") <- NULL
    expect_error(tcc_to_gene(se_tcc = se_copy),
                 regexp = 'is.element(el = "counts",',
                 fixed = TRUE)
    se_copy <- se_tcc
    SummarizedExperiment::assay(se_copy, "tpm") <- NULL
    expect_error(tcc_to_gene(se_tcc = se_copy),
                 regexp = 'is.element(el = "tpm",',
                 fixed = TRUE)
    se_copy <- se_tcc
    SummarizedExperiment::assay(se_copy, "relative_expression") <- NULL
    expect_error(tcc_to_gene(se_tcc = se_copy),
                 regexp = 'is.element(el = "relative_expression",',
                 fixed = TRUE)
    se_copy <- se_tcc
    SummarizedExperiment::rowData(se_copy)$gene_id <- NULL
    expect_error(tcc_to_gene(se_tcc = se_copy),
                 regexp = 'is.element(el = "gene_id",',
                 fixed = TRUE)
    se_copy <- se_tcc
    SummarizedExperiment::rowData(se_copy)$gene_name <- NULL
    expect_error(tcc_to_gene(se_tcc = se_copy),
                 regexp = 'is.element(el = "gene_name",',
                 fixed = TRUE)

    # Testing if function returns the expected output (bulk RNA-Seq data)
    expect_silent(
        se <- tcc_to_gene(se_tcc = se_tcc)
    )
    expect_true(class(se) == "SummarizedExperiment")
    expect_identical(
        dim(se),
        c(length(unique(S4Vectors::metadata(se_tcc)$transcript_df$gene_id)), ncol(se_tcc))
    )
    expect_identical(colnames(se), colnames(se_tcc))
    expect_identical(
        rownames(se),
        levels(factor(S4Vectors::metadata(se_tcc)$transcript_df$gene_id))
    )
    expect_identical(dim(SummarizedExperiment::colData(se)), c(ncol(se), 0L))
    expect_identical(colnames(SummarizedExperiment::rowData(se)),
                     c("gene_id", "gene_name"))
    expect_true(is.null(SummarizedExperiment::rowRanges(se)))
    expect_identical(SummarizedExperiment::assayNames(se),
                     c("counts", "tpm", "relative_expression"))
    expect_identical(class(SummarizedExperiment::assay(se, "counts")),
                     class(SummarizedExperiment::assay(se_tcc, "counts")))
    expect_identical(class(SummarizedExperiment::assay(se, "tpm")),
                     class(SummarizedExperiment::assay(se_tcc, "tpm")))
    expect_identical(class(SummarizedExperiment::assay(se, "relative_expression")),
                     class(SummarizedExperiment::assay(se_tcc, "relative_expression")))
    expect_identical(round(colSums(SummarizedExperiment::assay(se, "counts"))),
                     round(colSums(SummarizedExperiment::assay(se_tcc, "counts"))))
    expect_identical(round(colSums(SummarizedExperiment::assay(se, "tpm"))),
                     round(colSums(SummarizedExperiment::assay(se_tcc, "tpm"))))
    expect_identical(round(colSums(SummarizedExperiment::assay(se, "relative_expression"))),
                     round(colSums(SummarizedExperiment::assay(se_tcc, "relative_expression"))))
    expect_identical(S4Vectors::metadata(se), list())

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
    bam_parsed <- bam_to_read_structures(bam_files)
    transcript_data <- prepare_transcripts(
        gtf_file = gtf_file, genome_fasta_file = genome_fasta_file,
        bam_parsed = bam_parsed, min_bam_splice_read_count = 1,
        min_bam_splice_fraction = 0.01
    )
    se_tcc <- bam_to_tcc(
        bam_files = bam_files, transcript_data = transcript_data,
        is_single_cell = TRUE, barcode_tag = "BC",
        run_mode = "de_novo_loose", min_relative_expression = 0
    )

    # Testing if function returns the expected output (scRNA-Seq data)
    expect_silent(
        se <- tcc_to_gene(se_tcc = se_tcc)
    )
    expect_true(class(se) == "SummarizedExperiment")
    expect_identical(
        dim(se),
        c(length(unique(S4Vectors::metadata(se_tcc)$transcript_df$gene_id)), ncol(se_tcc))
    )
    expect_identical(colnames(se), colnames(se_tcc))
    expect_identical(
        rownames(se),
        levels(factor(S4Vectors::metadata(se_tcc)$transcript_df$gene_id))
    )
    expect_identical(dim(SummarizedExperiment::colData(se)), c(ncol(se), 0L))
    expect_identical(colnames(SummarizedExperiment::rowData(se)),
                     c("gene_id", "gene_name"))
    expect_true(is.null(SummarizedExperiment::rowRanges(se)))
    expect_identical(SummarizedExperiment::assayNames(se),
                     c("counts", "tpm", "relative_expression"))
    expect_identical(class(SummarizedExperiment::assay(se, "counts")),
                     class(SummarizedExperiment::assay(se_tcc, "counts")))
    expect_identical(class(SummarizedExperiment::assay(se, "tpm")),
                     class(SummarizedExperiment::assay(se_tcc, "tpm")))
    expect_identical(class(SummarizedExperiment::assay(se, "relative_expression")),
                     class(SummarizedExperiment::assay(se_tcc, "relative_expression")))
    expect_identical(round(Matrix::colSums(SummarizedExperiment::assay(se, "counts"))),
                     round(Matrix::colSums(SummarizedExperiment::assay(se_tcc, "counts"))))
    expect_identical(round(Matrix::colSums(SummarizedExperiment::assay(se, "tpm"))),
                     round(Matrix::colSums(SummarizedExperiment::assay(se_tcc, "tpm"))))
    expect_identical(round(Matrix::colSums(SummarizedExperiment::assay(se, "relative_expression"))),
                     round(Matrix::colSums(SummarizedExperiment::assay(se_tcc, "relative_expression"))))
    expect_identical(S4Vectors::metadata(se), list())
})
