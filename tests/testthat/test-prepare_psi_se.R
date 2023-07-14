test_that("prepare_psi_se works as expected", {

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
    bam_parsed <- extract_read_structures(bam_files)
    transcript_data <- prepare_transcripts(
        gtf_file = gtf_file, genome_fasta_file = genome_fasta_file,
        bam_parsed = bam_parsed, min_bam_splice_read_count = 2,
        min_bam_splice_fraction = 0.01
    )
    se_tcc <- prepare_tcc_se(
        bam_files = bam_files, transcript_data = transcript_data,
        run_mode = "de_novo_loose", min_relative_expression = 0
    )
    se_transcript <- prepare_transcript_se(
        se_tcc = se_tcc, use_length_normalization = TRUE
    )

    # Testing if function throws the expected errors
    expect_error(prepare_psi_se(se = NULL),
                 regexp = "methods::is(object = se, class2 =",
                 fixed = TRUE)
    se_copy <- se_transcript
    SummarizedExperiment::assay(se_copy, "relative_expression") <- NULL
    expect_error(prepare_psi_se(se = se_copy),
                 regexp = 'is.element(el = "relative_expression",',
                 fixed = TRUE)
    se_copy <- se_transcript
    SummarizedExperiment::rowData(se_copy)$gene_id <- NULL
    expect_error(prepare_psi_se(se = se_copy),
                 regexp = 'is.element(el = "gene_id",',
                 fixed = TRUE)
    se_copy <- se_transcript
    SummarizedExperiment::rowData(se_copy)$transcript_id <- NULL
    expect_error(prepare_psi_se(se = se_copy),
                 regexp = 'is.element(el = "transcript_id",',
                 fixed = TRUE)
    se_copy <- se_transcript
    SummarizedExperiment::rowData(se_copy)$position <- NULL
    expect_error(prepare_psi_se(se = se_copy),
                 regexp = 'is.element(el = "position",',
                 fixed = TRUE)
    se_copy <- se_transcript
    SummarizedExperiment::rowRanges(se_copy) <- NULL
    expect_error(prepare_psi_se(se = se_copy),
                 regexp = "rowRanges(se))) is not TRUE",
                 fixed = TRUE)

    # Testing if function returns the expected output (bulk RNA-Seq data)
    expect_silent(
        se <- prepare_psi_se(se = se_transcript)
    )
    expect_true(class(se) == "RangedSummarizedExperiment")
    expect_identical(dim(se), c(95L, 1L))
    expect_identical(colnames(se), colnames(se_transcript))
    expect_identical(colnames(SummarizedExperiment::rowData(se)),
                     c("gene_id", "type"))
    expect_identical(length(unique(SummarizedExperiment::rowData(se)$gene_id)),
                     3L)
    expect_identical(names(table(SummarizedExperiment::rowData(se)$type)),
                     c("A3", "A5", "CE", "RI", "TES", "TSS"))
    expect_identical(as.numeric(table(SummarizedExperiment::rowData(se)$type)),
                     c(5, 5, 20, 32, 14, 19))
    expect_true(grepl("GRanges", class(SummarizedExperiment::rowRanges(se))))
    expect_identical(SummarizedExperiment::assayNames(se),
                     c("psi"))
    expect_identical(class(SummarizedExperiment::assay(se, "psi")),
                     class(SummarizedExperiment::assay(se_transcript, "relative_expression")))
    expect_identical(round(colSums(SummarizedExperiment::assay(se, "psi"))),
                     c(Sample = 33))
    expect_true(all(SummarizedExperiment::assay(se, "psi") >= 0))
    expect_true(all(SummarizedExperiment::assay(se, "psi") <= 1))
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
    se_tcc <- se_tcc[, 1:5]
    se_transcript <- prepare_transcript_se(
        se_tcc = se_tcc, use_length_normalization = FALSE
    )
    se_transcript <- se_transcript[
        Matrix::rowSums(SummarizedExperiment::assay(se_transcript, "counts")) > 0,
    ]

    # Testing if function returns the expected output (scRNA-Seq data)
    expect_silent(
        se <- prepare_psi_se(se = se_transcript)
    )
    expect_true(class(se) == "RangedSummarizedExperiment")
    expect_identical(dim(se), c(11L, 5L))
    expect_identical(colnames(se), colnames(se_transcript))
    expect_identical(colnames(SummarizedExperiment::rowData(se)),
                     c("gene_id", "type"))
    expect_identical(length(unique(SummarizedExperiment::rowData(se)$gene_id)),
                     1L)
    expect_identical(names(table(SummarizedExperiment::rowData(se)$type)),
                     c("A5", "CE", "TES", "TSS"))
    expect_identical(as.numeric(table(SummarizedExperiment::rowData(se)$type)),
                     c(1, 6, 1, 3))
    expect_true(grepl("GRanges", class(SummarizedExperiment::rowRanges(se))))
    expect_identical(SummarizedExperiment::assayNames(se),
                     c("psi"))
    expect_identical(class(SummarizedExperiment::assay(se, "psi")),
                     class(SummarizedExperiment::assay(se_transcript, "relative_expression")))
    expect_identical(round(sum(Matrix::colSums(SummarizedExperiment::assay(se, "psi")))),
                     34)
    expect_true(all(SummarizedExperiment::assay(se, "psi") >= 0))
    expect_true(all(SummarizedExperiment::assay(se, "psi") <= 1))
    expect_identical(S4Vectors::metadata(se), list())
})
