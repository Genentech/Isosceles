test_that("pseudotime_tcc works as expected", {

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
    pseudotime <- seq(ncol(se_tcc)) / ncol(se_tcc)

    # Testing if function throws the expected errors
    expect_error(pseudotime_tcc(se_tcc = NULL),
                 regexp = "methods::is(object = se_tcc, class2 =",
                 fixed = TRUE)
    se_copy <- se_tcc
    SummarizedExperiment::assay(se_copy, "counts") <- NULL
    expect_error(pseudotime_tcc(se_tcc = se_copy),
                 regexp = 'is.element(el = "counts",',
                 fixed = TRUE)
    expect_error(pseudotime_tcc(se_tcc = se_tcc,
                                pseudotime = NULL),
                 regexp = "pseudotime is not a numeric or integer vector",
                 fixed = TRUE)
    expect_error(pseudotime_tcc(se_tcc = se_tcc,
                                pseudotime = pseudotime[1:10]),
                 regexp = "length(pseudotime) not identical to ncol(se_tcc)",
                 fixed = TRUE)
    expect_error(pseudotime_tcc(se_tcc = se_tcc,
                                pseudotime = pseudotime,
                                trim = NULL),
                 regexp = "trim is not a number",
                 fixed = TRUE)
    expect_error(pseudotime_tcc(se_tcc = se_tcc,
                                pseudotime = pseudotime,
                                trim = -1),
                 regexp = "trim not greater than or equal to 0",
                 fixed = TRUE)
    expect_error(pseudotime_tcc(se_tcc = se_tcc,
                                pseudotime = pseudotime,
                                trim = 1),
                 regexp = "trim not less than 0.5",
                 fixed = TRUE)
    expect_error(pseudotime_tcc(se_tcc = se_tcc,
                                pseudotime = pseudotime,
                                window_size = NULL),
                 regexp = "window_size is not a count",
                 fixed = TRUE)
    expect_error(pseudotime_tcc(se_tcc = se_tcc,
                                pseudotime = pseudotime,
                                window_size = ncol(se_tcc) + 10),
                 regexp = "window_size not less than or equal to ncol(se_tcc)",
                 fixed = TRUE)
    expect_error(pseudotime_tcc(se_tcc = se_tcc,
                                pseudotime = pseudotime,
                                window_step = NULL),
                 regexp = "window_step is not a count",
                 fixed = TRUE)
    expect_error(pseudotime_tcc(se_tcc = se_tcc,
                                pseudotime = pseudotime,
                                window_size = 3,
                                window_step = 6),
                 regexp = "window_step not less than or equal to window_size",
                 fixed = TRUE)

    # Testing if function returns the expected output
    # (window size = 3, window step = 3, trim = 0)
    expect_silent(
        se <- pseudotime_tcc(
            se_tcc = se_tcc, pseudotime = pseudotime, trim = 0,
            window_size = 3, window_step = 3
        )
    )
    assertthat::assert_that(methods::is(se, "SummarizedExperiment"))
    expect_identical(dim(se), c(8L, 19L))
    expect_true(all(grepl("^Cells_", colnames(se))))
    expect_identical(dim(SummarizedExperiment::colData(se)), c(19L, 1L))
    expect_identical(colnames(SummarizedExperiment::colData(se)),
                     "pseudotime")
    expect_true(all(SummarizedExperiment::colData(se)$pseudotime > 0))
    expect_true(all(SummarizedExperiment::colData(se)$pseudotime < 1))
    expect_identical(colnames(SummarizedExperiment::rowData(se)),
                     c("ec_id", "gene_id", "gene_name"))
    expect_identical(length(unique(SummarizedExperiment::rowData(se)$gene_id)),
                     1L)
    expect_true(is.null(SummarizedExperiment::rowRanges(se)))
    expect_identical(SummarizedExperiment::assayNames(se),
                     c("counts", "tpm", "relative_expression"))
    expect_true(class(SummarizedExperiment::assay(se, "counts")) == "dgCMatrix")
    expect_true(class(SummarizedExperiment::assay(se, "tpm")) == "dgCMatrix")
    expect_true(class(SummarizedExperiment::assay(se, "relative_expression")) == "dgCMatrix")
    expect_identical(round(sum(Matrix::colSums(SummarizedExperiment::assay(se, "counts")))),
                     round(sum(Matrix::colSums(SummarizedExperiment::assay(se_tcc, "counts")))))
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

    # Testing if function returns the expected output
    # (window size = 6, window step = 3, trim = 0)
    expect_silent(
        se <- pseudotime_tcc(
            se_tcc = se_tcc, pseudotime = pseudotime, trim = 0,
            window_size = 6, window_step = 3
        )
    )
    expect_true(class(se) == "SummarizedExperiment")
    expect_identical(dim(se), c(8L, 18L))
    expect_true(all(grepl("^Cells_", colnames(se))))
    expect_identical(dim(SummarizedExperiment::colData(se)), c(18L, 1L))
    expect_identical(colnames(SummarizedExperiment::colData(se)),
                     "pseudotime")
    expect_true(all(SummarizedExperiment::colData(se)$pseudotime > 0))
    expect_true(all(SummarizedExperiment::colData(se)$pseudotime < 1))
    expect_identical(colnames(SummarizedExperiment::rowData(se)),
                     c("ec_id", "gene_id", "gene_name"))
    expect_identical(length(unique(SummarizedExperiment::rowData(se)$gene_id)),
                     1L)
    expect_true(is.null(SummarizedExperiment::rowRanges(se)))
    expect_identical(SummarizedExperiment::assayNames(se),
                     c("counts", "tpm", "relative_expression"))
    expect_true(class(SummarizedExperiment::assay(se, "counts")) == "dgCMatrix")
    expect_true(class(SummarizedExperiment::assay(se, "tpm")) == "dgCMatrix")
    expect_true(class(SummarizedExperiment::assay(se, "relative_expression")) == "dgCMatrix")
    expect_identical(
        round(sum(Matrix::colSums(SummarizedExperiment::assay(se, "counts")))),
        2 * round(sum(Matrix::colSums(SummarizedExperiment::assay(se_tcc, "counts")))) -
            round(sum(Matrix::colSums(SummarizedExperiment::assay(se_tcc[, 1:3], "counts")))) -
            round(sum(Matrix::colSums(SummarizedExperiment::assay(se_tcc[, 55:57], "counts"))))
    )
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

    # Testing if function returns the expected output
    # (window size = 3, window step = 3, trim = 0.05)
    expect_silent(
        se <- pseudotime_tcc(
            se_tcc = se_tcc, pseudotime = pseudotime, trim = 0.05,
            window_size = 3, window_step = 3
        )
    )
    expect_true(class(se) == "SummarizedExperiment")
    expect_identical(dim(se), c(8L, 17L))
    expect_true(all(grepl("^Cells_", colnames(se))))
    expect_identical(dim(SummarizedExperiment::colData(se)), c(17L, 1L))
    expect_identical(colnames(SummarizedExperiment::colData(se)),
                     "pseudotime")
    expect_true(all(SummarizedExperiment::colData(se)$pseudotime > 0))
    expect_true(all(SummarizedExperiment::colData(se)$pseudotime < 1))
    expect_identical(colnames(SummarizedExperiment::rowData(se)),
                     c("ec_id", "gene_id", "gene_name"))
    expect_identical(length(unique(SummarizedExperiment::rowData(se)$gene_id)),
                     1L)
    expect_true(is.null(SummarizedExperiment::rowRanges(se)))
    expect_identical(SummarizedExperiment::assayNames(se),
                     c("counts", "tpm", "relative_expression"))
    expect_true(class(SummarizedExperiment::assay(se, "counts")) == "dgCMatrix")
    expect_true(class(SummarizedExperiment::assay(se, "tpm")) == "dgCMatrix")
    expect_true(class(SummarizedExperiment::assay(se, "relative_expression")) == "dgCMatrix")
    expect_identical(
        round(sum(Matrix::colSums(SummarizedExperiment::assay(se, "counts")))),
        round(sum(Matrix::colSums(SummarizedExperiment::assay(se_tcc, "counts")))) -
            round(sum(Matrix::colSums(SummarizedExperiment::assay(se_tcc[, 1:3], "counts")))) -
            round(sum(Matrix::colSums(SummarizedExperiment::assay(se_tcc[, 55:57], "counts"))))
    )
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
