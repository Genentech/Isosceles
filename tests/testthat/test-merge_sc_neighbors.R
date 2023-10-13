test_that("merge_sc_neighbors works as expected", {

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
    sce_tcc <- methods::as(se_tcc, "SingleCellExperiment")
    sce_tcc <- scuttle::computeLibraryFactors(sce_tcc)
    sce_tcc <- scuttle::logNormCounts(sce_tcc)
    set.seed(42)
    sce_tcc <- scater::runPCA(sce_tcc, ncomponents = 2)
    pca_mat <- SingleCellExperiment::reducedDim(sce_tcc, "PCA")

    # Testing if function throws the expected errors
    expect_error(merge_sc_neighbors(se_tcc = NULL),
                 regexp = "methods::is(object = se_tcc, class2 =",
                 fixed = TRUE)
    se_copy <- se_tcc
    SummarizedExperiment::assay(se_copy, "counts") <- NULL
    expect_error(merge_sc_neighbors(se_tcc = se_copy),
                 regexp = 'is.element(el = "counts",',
                 fixed = TRUE)
    expect_error(merge_sc_neighbors(se_tcc = se_tcc,
                                    pca_mat = NULL),
                 regexp = "pca_mat is not a matrix",
                 fixed = TRUE)
    expect_error(merge_sc_neighbors(se_tcc = se_tcc,
                                    pca_mat = pca_mat[1:10,]),
                 regexp = "colnames(se_tcc) not identical to rownames(pca_mat)",
                 fixed = TRUE)
    expect_error(merge_sc_neighbors(se_tcc = se_tcc,
                                    pca_mat = pca_mat,
                                    k = NULL),
                 regexp = "k is not a count",
                 fixed = TRUE)
    expect_error(merge_sc_neighbors(se_tcc = se_tcc,
                                    pca_mat = pca_mat,
                                    use_annoy = NULL),
                 regexp = "use_annoy is not a flag",
                 fixed = TRUE)
    expect_error(merge_sc_neighbors(se_tcc = se_tcc,
                                    pca_mat = pca_mat,
                                    ncpu = NULL),
                 regexp = "ncpu is not a count",
                 fixed = TRUE)

    # Testing if function returns the expected output
    expect_warning(
        se <- merge_sc_neighbors(
            se_tcc = se_tcc, pca_mat = pca_mat, k = 2,
            use_annoy = FALSE, ncpu = 1
        ),
        regexp = "tied distances",
        fixed = TRUE
    )
    expect_true(class(se) == "SummarizedExperiment")
    expect_identical(dim(se), c(8L, 57L))
    expect_true(all(grepl("^Sample\\.", colnames(se))))
    expect_identical(dim(SummarizedExperiment::colData(se)), c(57L, 0L))
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
    expect_true(all(
        Matrix::colSums(SummarizedExperiment::assay(se, "counts")) >
            Matrix::colSums(SummarizedExperiment::assay(se_tcc, "counts"))
    ))
    expect_identical(round(sum(Matrix::colSums(SummarizedExperiment::assay(se, "counts")))),
                     234)
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
