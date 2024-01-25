test_that("bam_to_tcc works as expected", {

    # Preparing test data (bulk RNA-Seq data)
    bam_file <- system.file(
        "extdata", "bulk_rnaseq.bam",
        package = "Isosceles"
    )
    bam_files <- c(Sample = bam_file)
    gtf_file <- system.file(
        "extdata", "bulk_rnaseq.gtf",
        package = "Isosceles"
    )
    genome_fasta_file <- system.file(
        "extdata", "bulk_rnaseq.fa",
        package = "Isosceles"
    )
    bam_parsed <- bam_to_read_structures(bam_files)
    transcript_data <- prepare_transcripts(
        gtf_file = gtf_file, genome_fasta_file = genome_fasta_file,
        bam_parsed = bam_parsed, min_bam_splice_read_count = 2,
        min_bam_splice_fraction = 0.01
    )

    # Testing if function throws the expected errors
    expect_error(bam_to_tcc(bam_files = NULL),
                 regexp = "bam_files is not a character vector",
                 fixed = TRUE)
    expect_error(bam_to_tcc(bam_files = character(0)),
                 regexp = "length(bam_files) not greater than 0",
                 fixed = TRUE)
    expect_error(bam_to_tcc(bam_files = "jabberwocky"),
                 regexp = "!is.null(names(bam_files)) is not TRUE",
                 fixed = TRUE)
    expect_error(bam_to_tcc(bam_files = c(Sample = "jabberwocky")),
                 regexp = "Elements 1 of file.exists(bam_files) are not true",
                 fixed = TRUE)
    expect_error(bam_to_tcc(bam_files = bam_files,
                            transcript_data = NULL),
                 regexp = "transcript_data is not a list",
                 fixed = TRUE)
    expect_error(bam_to_tcc(bam_files = bam_files,
                            transcript_data = list()),
                 regexp = "names(transcript_data) not identical to",
                 fixed = TRUE)
    expect_error(bam_to_tcc(bam_files = bam_files,
                            transcript_data = list(
                                tx_df = 42,
                                tx_granges = transcript_data$tx_granges,
                                tx_exon_granges_list = transcript_data$tx_exon_granges_list,
                                tx_intron_granges_list = transcript_data$tx_intron_granges_list
                            )),
                 regexp = "transcript_data$tx_df is not a data frame",
                 fixed = TRUE)
    expect_error(bam_to_tcc(bam_files = bam_files,
                            transcript_data = list(
                                tx_df = transcript_data$tx_df,
                                tx_granges = 42,
                                tx_exon_granges_list = transcript_data$tx_exon_granges_list,
                                tx_intron_granges_list = transcript_data$tx_intron_granges_list
                            )),
                 regexp = 'class(transcript_data$tx_granges) not equal to "GRanges"',
                 fixed = TRUE)
    expect_error(bam_to_tcc(bam_files = bam_files,
                            transcript_data = list(
                                tx_df = transcript_data$tx_df,
                                tx_granges = transcript_data$tx_granges,
                                tx_exon_granges_list = 42,
                                tx_intron_granges_list = transcript_data$tx_intron_granges_list
                            )),
                 regexp = 'grepl(pattern = "GRangesList", x = class',
                 fixed = TRUE)
    expect_error(bam_to_tcc(bam_files = bam_files,
                            transcript_data = list(
                                tx_df = transcript_data$tx_df,
                                tx_granges = transcript_data$tx_granges,
                                tx_exon_granges_list = transcript_data$tx_exon_granges_list,
                                tx_intron_granges_list = 42
                            )),
                 regexp = 'grepl(pattern = "GRangesList", x = class',
                 fixed = TRUE)
    expect_error(bam_to_tcc(bam_files = bam_files,
                            transcript_data = transcript_data,
                            run_mode = NULL),
                 regexp = "run_mode is not a string",
                 fixed = TRUE)
    expect_error(bam_to_tcc(bam_files = bam_files,
                            transcript_data = transcript_data,
                            run_mode = "jabberwocky"),
                 regexp = "`%in%`(x = run_mode, table = c(",
                 fixed = TRUE)
    expect_error(bam_to_tcc(bam_files = bam_files,
                            transcript_data = transcript_data,
                            min_read_count = NULL),
                 regexp = "min_read_count is not a count",
                 fixed = TRUE)
    expect_error(bam_to_tcc(bam_files = bam_files,
                            transcript_data = transcript_data,
                            min_relative_expression = NULL),
                 regexp = "min_relative_expression is not a numeric or integer vector",
                 fixed = TRUE)
    expect_error(bam_to_tcc(bam_files = bam_files,
                            transcript_data = transcript_data,
                            min_relative_expression = -42),
                 regexp = "min_relative_expression not greater than or equal to 0",
                 fixed = TRUE)
    expect_error(bam_to_tcc(bam_files = bam_files,
                            transcript_data = transcript_data,
                            extend_spliced_transcripts = -1),
                 regexp = "extend_spliced_transcripts not equal to 0 or extend_spliced_transcripts is not a count",
                 fixed = TRUE)
    expect_error(bam_to_tcc(bam_files = bam_files,
                            transcript_data = transcript_data,
                            is_single_cell = NULL),
                 regexp = "is_single_cell is not a flag",
                 fixed = TRUE)
    expect_error(bam_to_tcc(bam_files = bam_files,
                            transcript_data = transcript_data,
                            barcode_tag = NULL),
                 regexp = "barcode_tag is not a string",
                 fixed = TRUE)
    expect_error(bam_to_tcc(bam_files = bam_files,
                            transcript_data = transcript_data,
                            chunk_size = NULL),
                 regexp = "chunk_size is not a count",
                 fixed = TRUE)
    expect_error(bam_to_tcc(bam_files = bam_files,
                            transcript_data = transcript_data,
                            ncpu = NULL),
                 regexp = "ncpu is not a count",
                 fixed = TRUE)

    # Testing if function returns the expected output (bulk RNA-Seq data, strict mode)
    expect_message(
        se <- bam_to_tcc(
            bam_files = bam_files, transcript_data = transcript_data,
            run_mode = "strict"
        ),
        regexp = "read_id",
        fixed = TRUE
    )
    expect_true(class(se) == "SummarizedExperiment")
    expect_identical(dim(se), c(18L, 1L))
    expect_identical(colnames(se), names(bam_files))
    expect_identical(dim(SummarizedExperiment::colData(se)), c(1L, 0L))
    expect_identical(colnames(SummarizedExperiment::rowData(se)),
                     c("ec_id", "gene_id", "gene_name"))
    expect_identical(length(unique(SummarizedExperiment::rowData(se)$gene_id)),
                     4L)
    expect_true(is.null(SummarizedExperiment::rowRanges(se)))
    expect_identical(SummarizedExperiment::assayNames(se),
                     c("counts", "tpm", "relative_expression"))
    expect_true(is.matrix(SummarizedExperiment::assay(se, "counts")))
    expect_true(is.matrix(SummarizedExperiment::assay(se, "tpm")))
    expect_true(is.matrix(SummarizedExperiment::assay(se, "relative_expression")))
    expect_identical(round(colSums(SummarizedExperiment::assay(se, "counts"))),
                     c(Sample = 58))
    expect_identical(round(colSums(SummarizedExperiment::assay(se, "tpm"))),
                     c(Sample = 1e6))
    expect_identical(round(colSums(SummarizedExperiment::assay(se, "relative_expression"))),
                     c(Sample = 4))
    expect_true(is.list(S4Vectors::metadata(se)))
    expect_identical(names(S4Vectors::metadata(se)),
                     c("compatibility_matrix", "transcript_df",
                       "transcript_exon_granges_list", "mean_read_length"))
    expect_true(class(S4Vectors::metadata(se)$compatibility_matrix) == "dgCMatrix")
    expect_identical(dim(S4Vectors::metadata(se)$compatibility_matrix),
                     c(18L, 105L))
    expect_true(is.data.frame(S4Vectors::metadata(se)$transcript_df))
    expect_identical(dim(S4Vectors::metadata(se)$transcript_df),
                     c(105L, 12L))
    expect_identical(colnames(S4Vectors::metadata(se)$transcript_df),
                     colnames(transcript_data$tx_df))
    expect_true(all(S4Vectors::metadata(se)$transcript_df$fivethree_support_level == "FL"))
    expect_true(all(S4Vectors::metadata(se)$transcript_df$splicing_support_level == "AP"))
    expect_true(grepl("GRangesList", class(S4Vectors::metadata(se)$transcript_exon_granges_list)))
    expect_identical(length(S4Vectors::metadata(se)$transcript_exon_granges_list),
                     105L)
    expect_identical(length(unlist(S4Vectors::metadata(se)$transcript_exon_granges_list)),
                     831L)
    expect_identical(round(S4Vectors::metadata(se)$mean_read_length), 1315)

    # Testing if function returns the expected output (bulk RNA-Seq data, de_novo_strict mode)
    expect_message(
        se <- bam_to_tcc(
            bam_files = bam_files, transcript_data = transcript_data,
            run_mode = "de_novo_strict", min_relative_expression = 0
        ),
        regexp = "read_id",
        fixed = TRUE
    )
    expect_true(class(se) == "SummarizedExperiment")
    expect_identical(dim(se), c(20L, 1L))
    expect_identical(colnames(se), names(bam_files))
    expect_identical(dim(SummarizedExperiment::colData(se)), c(1L, 0L))
    expect_identical(colnames(SummarizedExperiment::rowData(se)),
                     c("ec_id", "gene_id", "gene_name"))
    expect_identical(length(unique(SummarizedExperiment::rowData(se)$gene_id)),
                     4L)
    expect_true(is.null(SummarizedExperiment::rowRanges(se)))
    expect_identical(SummarizedExperiment::assayNames(se),
                     c("counts", "tpm", "relative_expression"))
    expect_true(is.matrix(SummarizedExperiment::assay(se, "counts")))
    expect_true(is.matrix(SummarizedExperiment::assay(se, "tpm")))
    expect_true(is.matrix(SummarizedExperiment::assay(se, "relative_expression")))
    expect_identical(round(colSums(SummarizedExperiment::assay(se, "counts"))),
                     c(Sample = 60))
    expect_identical(round(colSums(SummarizedExperiment::assay(se, "tpm"))),
                     c(Sample = 1e6))
    expect_identical(round(colSums(SummarizedExperiment::assay(se, "relative_expression"))),
                     c(Sample = 4))
    expect_true(is.list(S4Vectors::metadata(se)))
    expect_identical(names(S4Vectors::metadata(se)),
                     c("compatibility_matrix", "transcript_df",
                       "transcript_exon_granges_list", "mean_read_length"))
    expect_true(class(S4Vectors::metadata(se)$compatibility_matrix) == "dgCMatrix")
    expect_identical(dim(S4Vectors::metadata(se)$compatibility_matrix),
                     c(20L, 107L))
    expect_true(is.data.frame(S4Vectors::metadata(se)$transcript_df))
    expect_identical(dim(S4Vectors::metadata(se)$transcript_df),
                     c(107L, 12L))
    expect_identical(colnames(S4Vectors::metadata(se)$transcript_df),
                     colnames(transcript_data$tx_df))
    expect_true(all(S4Vectors::metadata(se)$transcript_df$fivethree_support_level == "FL"))
    expect_true(all(S4Vectors::metadata(se)$transcript_df$splicing_support_level %in%
                        c("AP", "EC", "NC")))
    expect_true(grepl("GRangesList", class(S4Vectors::metadata(se)$transcript_exon_granges_list)))
    expect_identical(length(S4Vectors::metadata(se)$transcript_exon_granges_list),
                     107L)
    expect_identical(length(unlist(S4Vectors::metadata(se)$transcript_exon_granges_list)),
                     855L)
    expect_identical(round(S4Vectors::metadata(se)$mean_read_length), 1315)

    # Testing if function returns the expected output (bulk RNA-Seq data, de_novo_loose mode)
    expect_message(
        se <- bam_to_tcc(
            bam_files = bam_files, transcript_data = transcript_data,
            run_mode = "de_novo_loose", min_relative_expression = 0
        ),
        regexp = "read_id",
        fixed = TRUE
    )
    expect_true(class(se) == "SummarizedExperiment")
    expect_identical(dim(se), c(26L, 1L))
    expect_identical(colnames(se), names(bam_files))
    expect_identical(dim(SummarizedExperiment::colData(se)), c(1L, 0L))
    expect_identical(colnames(SummarizedExperiment::rowData(se)),
                     c("ec_id", "gene_id", "gene_name"))
    expect_identical(length(unique(SummarizedExperiment::rowData(se)$gene_id)),
                     4L)
    expect_true(is.null(SummarizedExperiment::rowRanges(se)))
    expect_identical(SummarizedExperiment::assayNames(se),
                     c("counts", "tpm", "relative_expression"))
    expect_true(is.matrix(SummarizedExperiment::assay(se, "counts")))
    expect_true(is.matrix(SummarizedExperiment::assay(se, "tpm")))
    expect_true(is.matrix(SummarizedExperiment::assay(se, "relative_expression")))
    expect_identical(round(colSums(SummarizedExperiment::assay(se, "counts"))),
                     c(Sample = 69))
    expect_identical(round(colSums(SummarizedExperiment::assay(se, "tpm"))),
                     c(Sample = 1e6))
    expect_identical(round(colSums(SummarizedExperiment::assay(se, "relative_expression"))),
                     c(Sample = 4))
    expect_true(is.list(S4Vectors::metadata(se)))
    expect_identical(names(S4Vectors::metadata(se)),
                     c("compatibility_matrix", "transcript_df",
                       "transcript_exon_granges_list", "mean_read_length"))
    expect_true(class(S4Vectors::metadata(se)$compatibility_matrix) == "dgCMatrix")
    expect_identical(dim(S4Vectors::metadata(se)$compatibility_matrix),
                     c(26L, 112L))
    expect_true(is.data.frame(S4Vectors::metadata(se)$transcript_df))
    expect_identical(dim(S4Vectors::metadata(se)$transcript_df),
                     c(112L, 12L))
    expect_identical(colnames(S4Vectors::metadata(se)$transcript_df),
                     colnames(transcript_data$tx_df))
    expect_true(all(S4Vectors::metadata(se)$transcript_df$fivethree_support_level == "FL"))
    expect_true(all(S4Vectors::metadata(se)$transcript_df$splicing_support_level %in%
                        c("AP", "EC", "NC", "DN")))
    expect_true(grepl("GRangesList", class(S4Vectors::metadata(se)$transcript_exon_granges_list)))
    expect_identical(length(S4Vectors::metadata(se)$transcript_exon_granges_list),
                     112L)
    expect_identical(length(unlist(S4Vectors::metadata(se)$transcript_exon_granges_list)),
                     909L)
    expect_identical(round(S4Vectors::metadata(se)$mean_read_length), 1315)

    # Testing if function returns the expected output (bulk RNA-Seq data, de_novo_full mode)
    expect_warning(
        se <- bam_to_tcc(
            bam_files = bam_files, transcript_data = transcript_data,
            run_mode = "de_novo_full", min_relative_expression = 0
        ),
        regexp = "experimental feature",
        fixed = TRUE
    )
    expect_true(class(se) == "SummarizedExperiment")
    expect_identical(dim(se), c(24L, 1L))
    expect_identical(colnames(se), names(bam_files))
    expect_identical(dim(SummarizedExperiment::colData(se)), c(1L, 0L))
    expect_identical(colnames(SummarizedExperiment::rowData(se)),
                     c("ec_id", "gene_id", "gene_name"))
    expect_identical(length(unique(SummarizedExperiment::rowData(se)$gene_id)),
                     2L)
    expect_true(is.null(SummarizedExperiment::rowRanges(se)))
    expect_identical(SummarizedExperiment::assayNames(se),
                     c("counts", "tpm", "relative_expression"))
    expect_true(is.matrix(SummarizedExperiment::assay(se, "counts")))
    expect_true(is.matrix(SummarizedExperiment::assay(se, "tpm")))
    expect_true(is.matrix(SummarizedExperiment::assay(se, "relative_expression")))
    expect_identical(round(colSums(SummarizedExperiment::assay(se, "counts"))),
                     c(Sample = 59))
    expect_identical(round(colSums(SummarizedExperiment::assay(se, "tpm"))),
                     c(Sample = 1e6))
    expect_identical(round(colSums(SummarizedExperiment::assay(se, "relative_expression"))),
                     c(Sample = 2))
    expect_true(is.list(S4Vectors::metadata(se)))
    expect_identical(names(S4Vectors::metadata(se)),
                     c("compatibility_matrix", "transcript_df",
                       "transcript_exon_granges_list", "mean_read_length"))
    expect_true(class(S4Vectors::metadata(se)$compatibility_matrix) == "dgCMatrix")
    expect_identical(dim(S4Vectors::metadata(se)$compatibility_matrix),
                     c(24L, 17L))
    expect_true(is.data.frame(S4Vectors::metadata(se)$transcript_df))
    expect_identical(dim(S4Vectors::metadata(se)$transcript_df),
                     c(17L, 12L))
    expect_identical(colnames(S4Vectors::metadata(se)$transcript_df),
                     colnames(transcript_data$tx_df))
    expect_true(all(S4Vectors::metadata(se)$transcript_df$fivethree_support_level == "FL"))
    expect_true(all(S4Vectors::metadata(se)$transcript_df$splicing_support_level %in%
                        c("PC", "EC", "NC", "DN")))
    expect_true(grepl("GRangesList", class(S4Vectors::metadata(se)$transcript_exon_granges_list)))
    expect_identical(length(S4Vectors::metadata(se)$transcript_exon_granges_list),
                     17L)
    expect_identical(length(unlist(S4Vectors::metadata(se)$transcript_exon_granges_list)),
                     175L)
    expect_identical(round(S4Vectors::metadata(se)$mean_read_length), 1315)

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
    bam_parsed <- bam_to_read_structures(bam_files)
    transcript_data <- prepare_transcripts(
        gtf_file = gtf_file, genome_fasta_file = genome_fasta_file,
        bam_parsed = bam_parsed, min_bam_splice_read_count = 1,
        min_bam_splice_fraction = 0.01
    )

    # Testing if function returns the expected output (scRNA-Seq data, strict mode)
    expect_message(
        se <- bam_to_tcc(
            bam_files = bam_files, transcript_data = transcript_data,
            is_single_cell = TRUE, barcode_tag = "BC",
            run_mode = "strict"
        ),
        regexp = "read_id",
        fixed = TRUE
    )
    expect_true(class(se) == "SummarizedExperiment")
    expect_identical(dim(se), c(6L, 54L))
    expect_true(all(grepl("^Sample\\.", colnames(se))))
    expect_identical(dim(SummarizedExperiment::colData(se)), c(54L, 0L))
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
                     67)
    expect_true(all(round(Matrix::colSums(SummarizedExperiment::assay(se, "tpm"))) == 1e6))
    expect_true(all(round(Matrix::colSums(SummarizedExperiment::assay(se, "relative_expression"))) == 1))
    expect_true(is.list(S4Vectors::metadata(se)))
    expect_identical(names(S4Vectors::metadata(se)),
                     c("compatibility_matrix", "transcript_df",
                       "transcript_exon_granges_list", "mean_read_length"))
    expect_true(class(S4Vectors::metadata(se)$compatibility_matrix) == "dgCMatrix")
    expect_identical(dim(S4Vectors::metadata(se)$compatibility_matrix),
                     c(6L, 4109L))
    expect_true(is.data.frame(S4Vectors::metadata(se)$transcript_df))
    expect_identical(dim(S4Vectors::metadata(se)$transcript_df),
                     c(4109L, 12L))
    expect_identical(colnames(S4Vectors::metadata(se)$transcript_df),
                     colnames(transcript_data$tx_df))
    expect_true(all(S4Vectors::metadata(se)$transcript_df$fivethree_support_level == "FL"))
    expect_true(all(S4Vectors::metadata(se)$transcript_df$splicing_support_level == "AP"))
    expect_true(grepl("GRangesList", class(S4Vectors::metadata(se)$transcript_exon_granges_list)))
    expect_identical(length(S4Vectors::metadata(se)$transcript_exon_granges_list),
                     4109L)
    expect_identical(length(unlist(S4Vectors::metadata(se)$transcript_exon_granges_list)),
                     29269L)
    expect_identical(round(S4Vectors::metadata(se)$mean_read_length), 1014)

    # Testing if function returns the expected output (scRNA-Seq data, de_novo_strict mode)
    expect_message(
        se <- bam_to_tcc(
            bam_files = bam_files, transcript_data = transcript_data,
            is_single_cell = TRUE, barcode_tag = "BC",
            run_mode = "de_novo_strict", min_relative_expression = 0
        ),
        regexp = "read_id",
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
    expect_identical(round(sum(Matrix::colSums(SummarizedExperiment::assay(se, "counts")))),
                     72)
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
                        c("AP", "EC", "NC")))
    expect_true(grepl("GRangesList", class(S4Vectors::metadata(se)$transcript_exon_granges_list)))
    expect_identical(length(S4Vectors::metadata(se)$transcript_exon_granges_list),
                     4111L)
    expect_identical(length(unlist(S4Vectors::metadata(se)$transcript_exon_granges_list)),
                     29282L)
    expect_identical(round(S4Vectors::metadata(se)$mean_read_length), 1014)

    # Testing if function returns the expected output (scRNA-Seq data, de_novo_loose mode)
    expect_message(
        se <- bam_to_tcc(
            bam_files = bam_files, transcript_data = transcript_data,
            is_single_cell = TRUE, barcode_tag = "BC",
            run_mode = "de_novo_loose", min_relative_expression = 0
        ),
        regexp = "read_id",
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
    expect_identical(round(sum(Matrix::colSums(SummarizedExperiment::assay(se, "counts")))),
                     72)
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

    # Testing if function returns the expected output (scRNA-Seq data, de_novo_full mode)
    expect_warning(
        se <- bam_to_tcc(
            bam_files = bam_files, transcript_data = transcript_data,
            is_single_cell = TRUE, barcode_tag = "BC",
            run_mode = "de_novo_full", min_relative_expression = 0
        ),
        regexp = "experimental feature",
        fixed = TRUE
    )
    expect_true(class(se) == "SummarizedExperiment")
    expect_identical(dim(se), c(7L, 55L))
    expect_true(all(grepl("^Sample\\.", colnames(se))))
    expect_identical(dim(SummarizedExperiment::colData(se)), c(55L, 0L))
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
                     66)
    expect_true(all(round(Matrix::colSums(SummarizedExperiment::assay(se, "tpm"))) == 1e6))
    expect_true(all(round(Matrix::colSums(SummarizedExperiment::assay(se, "relative_expression"))) == 1))
    expect_true(is.list(S4Vectors::metadata(se)))
    expect_identical(names(S4Vectors::metadata(se)),
                     c("compatibility_matrix", "transcript_df",
                       "transcript_exon_granges_list", "mean_read_length"))
    expect_true(class(S4Vectors::metadata(se)$compatibility_matrix) == "dgCMatrix")
    expect_identical(dim(S4Vectors::metadata(se)$compatibility_matrix),
                     c(7L, 6L))
    expect_true(is.data.frame(S4Vectors::metadata(se)$transcript_df))
    expect_identical(dim(S4Vectors::metadata(se)$transcript_df),
                     c(6L, 12L))
    expect_identical(colnames(S4Vectors::metadata(se)$transcript_df),
                     colnames(transcript_data$tx_df))
    expect_true(all(S4Vectors::metadata(se)$transcript_df$fivethree_support_level == "FL"))
    expect_true(all(S4Vectors::metadata(se)$transcript_df$splicing_support_level %in%
                        c("PC", "EC", "NC", "DN")))
    expect_true(grepl("GRangesList", class(S4Vectors::metadata(se)$transcript_exon_granges_list)))
    expect_identical(length(S4Vectors::metadata(se)$transcript_exon_granges_list),
                     6L)
    expect_identical(length(unlist(S4Vectors::metadata(se)$transcript_exon_granges_list)),
                     35L)
    expect_identical(round(S4Vectors::metadata(se)$mean_read_length), 1014)
})
