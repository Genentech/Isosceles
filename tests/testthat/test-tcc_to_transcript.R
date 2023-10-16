test_that("tcc_to_transcript works as expected", {

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

    # Testing if function throws the expected errors
    expect_error(tcc_to_transcript(se_tcc = NULL),
                 regexp = "methods::is(object = se_tcc, class2 =",
                 fixed = TRUE)
    se_copy <- se_tcc
    SummarizedExperiment::assay(se_copy, "counts") <- NULL
    expect_error(tcc_to_transcript(se_tcc = se_copy),
                 regexp = 'is.element(el = "counts",',
                 fixed = TRUE)
    se_copy <- se_tcc
    SummarizedExperiment::rowData(se_copy)$gene_id <- NULL
    expect_error(tcc_to_transcript(se_tcc = se_copy),
                 regexp = 'is.element(el = "gene_id",',
                 fixed = TRUE)
    expect_error(tcc_to_transcript(se_tcc = se_tcc,
                                   em.maxiter = NULL),
                 regexp = "em.maxiter is not a count",
                 fixed = TRUE)
    expect_error(tcc_to_transcript(se_tcc = se_tcc,
                                   em.conv = NULL),
                 regexp = "em.conv is not a number",
                 fixed = TRUE)
    expect_error(tcc_to_transcript(se_tcc = se_tcc,
                                   ncpu = NULL),
                 regexp = "ncpu is not a count",
                 fixed = TRUE)
    expect_error(tcc_to_transcript(se_tcc = se_tcc,
                                   use_length_normalization = NULL),
                 regexp = "use_length_normalization is not a flag",
                 fixed = TRUE)

    # Testing if function returns the expected output (bulk RNA-Seq data, length normalization)
    expect_silent(
        se <- tcc_to_transcript(
            se_tcc = se_tcc, use_length_normalization = TRUE
        )
    )
    expect_true(class(se) == "RangedSummarizedExperiment")
    expect_identical(
        dim(se),
        c(nrow(S4Vectors::metadata(se_tcc)$transcript_df), ncol(se_tcc))
    )
    expect_identical(colnames(se), colnames(se_tcc))
    expect_identical(
        rownames(se),
        S4Vectors::metadata(se_tcc)$transcript_df$transcript_id
    )
    expect_identical(dim(SummarizedExperiment::colData(se)), c(ncol(se), 0L))
    expect_identical(colnames(SummarizedExperiment::rowData(se)),
                     colnames(S4Vectors::metadata(se_tcc)$transcript_df))
    expect_identical(length(unique(SummarizedExperiment::rowData(se)$gene_id)),
                     length(unique(S4Vectors::metadata(se_tcc)$transcript_df$gene_id)))
    expect_true(grepl("GRangesList", class(SummarizedExperiment::rowRanges(se))))
    expect_identical(length(SummarizedExperiment::rowRanges(se)),
                     length(S4Vectors::metadata(se_tcc)$transcript_exon_granges_list))
    expect_identical(length(unlist(SummarizedExperiment::rowRanges(se))),
                     length(unlist(S4Vectors::metadata(se_tcc)$transcript_exon_granges_list)))
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

    # Testing if function returns the expected output (bulk RNA-Seq data, no length normalization)
    expect_silent(
        se <- tcc_to_transcript(
            se_tcc = se_tcc, use_length_normalization = FALSE
        )
    )
    expect_true(class(se) == "RangedSummarizedExperiment")
    expect_identical(
        dim(se),
        c(nrow(S4Vectors::metadata(se_tcc)$transcript_df), ncol(se_tcc))
    )
    expect_identical(colnames(se), colnames(se_tcc))
    expect_identical(
        rownames(se),
        S4Vectors::metadata(se_tcc)$transcript_df$transcript_id
    )
    expect_identical(dim(SummarizedExperiment::colData(se)), c(ncol(se), 0L))
    expect_identical(colnames(SummarizedExperiment::rowData(se)),
                     colnames(S4Vectors::metadata(se_tcc)$transcript_df))
    expect_identical(length(unique(SummarizedExperiment::rowData(se)$gene_id)),
                     length(unique(S4Vectors::metadata(se_tcc)$transcript_df$gene_id)))
    expect_true(grepl("GRangesList", class(SummarizedExperiment::rowRanges(se))))
    expect_identical(length(SummarizedExperiment::rowRanges(se)),
                     length(S4Vectors::metadata(se_tcc)$transcript_exon_granges_list))
    expect_identical(length(unlist(SummarizedExperiment::rowRanges(se))),
                     length(unlist(S4Vectors::metadata(se_tcc)$transcript_exon_granges_list)))
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

    # Testing if function returns the expected output (scRNA-Seq, length normalization)
    expect_silent(
        se <- tcc_to_transcript(
            se_tcc = se_tcc, use_length_normalization = TRUE
        )
    )
    expect_true(class(se) == "RangedSummarizedExperiment")
    expect_identical(
        dim(se),
        c(nrow(S4Vectors::metadata(se_tcc)$transcript_df), ncol(se_tcc))
    )
    expect_identical(colnames(se), colnames(se_tcc))
    expect_identical(
        rownames(se),
        S4Vectors::metadata(se_tcc)$transcript_df$transcript_id
    )
    expect_identical(dim(SummarizedExperiment::colData(se)), c(ncol(se), 0L))
    expect_identical(colnames(SummarizedExperiment::rowData(se)),
                     colnames(S4Vectors::metadata(se_tcc)$transcript_df))
    expect_identical(length(unique(SummarizedExperiment::rowData(se)$gene_id)),
                     length(unique(S4Vectors::metadata(se_tcc)$transcript_df$gene_id)))
    expect_true(grepl("GRangesList", class(SummarizedExperiment::rowRanges(se))))
    expect_identical(length(SummarizedExperiment::rowRanges(se)),
                     length(S4Vectors::metadata(se_tcc)$transcript_exon_granges_list))
    expect_identical(length(unlist(SummarizedExperiment::rowRanges(se))),
                     length(unlist(S4Vectors::metadata(se_tcc)$transcript_exon_granges_list)))
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

    # Testing if function returns the expected output (scRNA-Seq, no length normalization)
    expect_silent(
        se <- tcc_to_transcript(
            se_tcc = se_tcc, use_length_normalization = FALSE
        )
    )
    expect_true(class(se) == "RangedSummarizedExperiment")
    expect_identical(
        dim(se),
        c(nrow(S4Vectors::metadata(se_tcc)$transcript_df), ncol(se_tcc))
    )
    expect_identical(colnames(se), colnames(se_tcc))
    expect_identical(
        rownames(se),
        S4Vectors::metadata(se_tcc)$transcript_df$transcript_id
    )
    expect_identical(dim(SummarizedExperiment::colData(se)), c(ncol(se), 0L))
    expect_identical(colnames(SummarizedExperiment::rowData(se)),
                     colnames(S4Vectors::metadata(se_tcc)$transcript_df))
    expect_identical(length(unique(SummarizedExperiment::rowData(se)$gene_id)),
                     length(unique(S4Vectors::metadata(se_tcc)$transcript_df$gene_id)))
    expect_true(grepl("GRangesList", class(SummarizedExperiment::rowRanges(se))))
    expect_identical(length(SummarizedExperiment::rowRanges(se)),
                     length(S4Vectors::metadata(se_tcc)$transcript_exon_granges_list))
    expect_identical(length(unlist(SummarizedExperiment::rowRanges(se))),
                     length(unlist(S4Vectors::metadata(se_tcc)$transcript_exon_granges_list)))
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
