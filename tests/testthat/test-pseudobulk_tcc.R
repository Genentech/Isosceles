test_that("pseudobulk_tcc works as expected", {

    # Preparing test data
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
    set.seed(42)
    cell_labels <- sample(1:2, ncol(se_tcc), replace = TRUE)

    # Testing if function throws the expected errors
    expect_error(pseudobulk_tcc(se_tcc = NULL),
                 regexp = "methods::is(object = se_tcc, class2 =",
                 fixed = TRUE)
    se_copy <- se_tcc
    SummarizedExperiment::assay(se_copy, "counts") <- NULL
    expect_error(pseudobulk_tcc(se_tcc = se_copy),
                 regexp = 'is.element(el = "counts",',
                 fixed = TRUE)
    expect_error(pseudobulk_tcc(se_tcc = se_tcc,
                                cell_labels = NULL),
                 regexp = "length(cell_labels) not equal to ncol(se_tcc)",
                 fixed = TRUE)

    # Testing if function returns the expected output
    expect_silent(
        se <- pseudobulk_tcc(
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
    expect_true(is.list(S4Vectors::metadata(se)))
    expect_identical(names(S4Vectors::metadata(se)),
                     c("compatibility_matrix", "transcript_df",
                       "transcript_exon_granges_list", "mean_read_length"))
    expect_true(class(S4Vectors::metadata(se)$compatibility_matrix) == "dgCMatrix")
    expect_identical(dim(S4Vectors::metadata(se)$compatibility_matrix),
                     c(8L, 4111L))
    expect_true(is.data.frame(S4Vectors::metadata(se)$transcript_df))
    expect_identical(dim(S4Vectors::metadata(se)$transcript_df),
                     c(4111L, 12L))
    expect_identical(colnames(S4Vectors::metadata(se)$transcript_df),
                     colnames(transcript_data$tx_df))
    expect_true(all(S4Vectors::metadata(se)$transcript_df$fivethree_support_level == "FL"))
    expect_true(all(S4Vectors::metadata(se)$transcript_df$splicing_support_level %in%
                        c("AP", "EC", "NC", "DN")))
    expect_true(grepl("GRangesList", class(S4Vectors::metadata(se)$transcript_exon_granges_list)))
    expect_identical(length(S4Vectors::metadata(se)$transcript_exon_granges_list),
                     4111L)
    expect_identical(length(unlist(S4Vectors::metadata(se)$transcript_exon_granges_list)),
                     29282L)
    expect_identical(round(S4Vectors::metadata(se)$mean_read_length), 1014)
})
