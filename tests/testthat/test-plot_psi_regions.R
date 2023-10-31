test_that("plot_psi_regions works as expected", {

    # Preparing test data
    sce_transcript <- readRDS(system.file("extdata", "sce_transcript_cav1.rds",
                                          package = "Isosceles"))
    se_psi <- transcript_to_psi(sce_transcript)

    # Testing if function throws the expected errors
    expect_error(plot_psi_regions(se_psi = NULL),
                 regexp = "methods::is(object = se_psi, class2 =",
                 fixed = TRUE)
    se_copy <- se_psi
    SummarizedExperiment::assay(se_copy, "psi") <- NULL
    expect_error(plot_psi_regions(se_psi = se_copy),
                 regexp = 'is.element(el = "psi",',
                 fixed = TRUE)
    se_copy <- se_psi
    SummarizedExperiment::rowData(se_copy)$gene_id <- NULL
    expect_error(plot_psi_regions(se_psi = se_copy),
                 regexp = 'is.element(el = "gene_id",',
                 fixed = TRUE)
    expect_error(plot_psi_regions(se_psi = se_psi,
                                  se_transcript = NULL),
                 regexp = "methods::is(object = se_transcript, class2 =",
                 fixed = TRUE)
    se_copy <- sce_transcript
    SummarizedExperiment::assay(se_copy, "tpm") <- NULL
    expect_error(plot_psi_regions(se_psi = se_psi,
                                  se_transcript = se_copy),
                 regexp = 'is.element(el = "tpm",',
                 fixed = TRUE)
    se_copy <- sce_transcript
    SummarizedExperiment::rowData(se_copy)$gene_id <- NULL
    expect_error(plot_psi_regions(se_psi = se_psi,
                                  se_transcript = se_copy),
                 regexp = 'is.element(el = "gene_id",',
                 fixed = TRUE)
    expect_error(plot_psi_regions(se_psi = se_psi,
                                  se_transcript = sce_transcript,
                                  gene_id = NULL),
                 regexp = "gene_id is not a string",
                 fixed = TRUE)
    expect_error(plot_psi_regions(se_psi = se_psi,
                                  se_transcript = sce_transcript,
                                  gene_id = "jabberwocky"),
                 regexp = "rowData(se_psi)$gene_id) is not TRUE",
                 fixed = TRUE)
    se_copy <- se_psi
    SummarizedExperiment::rowData(se_copy)$gene_id <- "jabberwocky"
    expect_error(plot_psi_regions(se_psi = se_copy,
                                  se_transcript = sce_transcript,
                                  gene_id = "jabberwocky"),
                 regexp = "rowData(se_transcript)$gene_id) is not TRUE",
                 fixed = TRUE)
    expect_error(plot_psi_regions(se_psi = se_psi,
                                  se_transcript = sce_transcript,
                                  gene_id = "ENSG00000105974",
                                  max_transcripts = NULL),
                 regexp = "max_transcripts is not a count",
                 fixed = TRUE)
    expect_error(plot_psi_regions(se_psi = se_psi,
                                  se_transcript = sce_transcript,
                                  gene_id = "ENSG00000105974",
                                  max_intron_length = NA),
                 regexp = "max_intron_length is not a count",
                 fixed = TRUE)
    expect_error(plot_psi_regions(se_psi = se_psi,
                                  se_transcript = sce_transcript,
                                  gene_id = "ENSG00000105974",
                                  region_colors = NA),
                 regexp = "region_colors is not a character vector",
                 fixed = TRUE)
    expect_error(plot_psi_regions(se_psi = se_psi,
                                  se_transcript = sce_transcript,
                                  gene_id = "ENSG00000105974",
                                  region_colors = (TSS = "#FF83FF")),
                 regexp = "sort(names(region_colors)) not identical to",
                 fixed = TRUE)

    # Testing if function returns the expected output
    suppressWarnings(expect_warning(
        plot_data <- plot_psi_regions(
            se_psi = se_psi,
            se_transcript = sce_transcript,
            gene_id = "ENSG00000105974"
        ),
        regexp = "Ignoring unknown aesthetics",
        fixed = TRUE
    ))
    expect_true(methods::is(plot_data, "Tracks"))
    
    # Testing if function returns the expected output (intron shrinking)
    suppressWarnings(expect_warning(
      plot_data <- plot_psi_regions(
        se_psi = se_psi,
        se_transcript = sce_transcript,
        gene_id = "ENSG00000105974",
        max_intron_length = 1000
      ),
      regexp = "Ignoring unknown aesthetics",
      fixed = TRUE
    ))
    expect_true(methods::is(plot_data, "Tracks"))
})
