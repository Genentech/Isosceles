test_that("plot_psi_heatmap works as expected", {

    # Preparing test data
    sce_transcript <- readRDS(system.file("extdata", "sce_transcript_cav1.rds",
                                          package = "Isosceles"))
    se_psi <- transcript_to_psi(sce_transcript)

    # Testing if function throws the expected errors
    expect_error(plot_psi_heatmap(se_psi = NULL),
                 regexp = "methods::is(object = se_psi, class2 =",
                 fixed = TRUE)
    se_copy <- se_psi
    SummarizedExperiment::assay(se_copy, "psi") <- NULL
    expect_error(plot_psi_heatmap(se_psi = se_copy),
                 regexp = 'is.element(el = "psi",',
                 fixed = TRUE)
    se_copy <- se_psi
    SummarizedExperiment::rowData(se_copy)$gene_id <- NULL
    expect_error(plot_psi_heatmap(se_psi = se_copy),
                 regexp = 'is.element(el = "gene_id",',
                 fixed = TRUE)
    expect_error(plot_psi_heatmap(se_psi = se_psi,
                                  gene_id = NULL),
                 regexp = "gene_id is not a string",
                 fixed = TRUE)
    expect_error(plot_psi_heatmap(se_psi = se_psi,
                                  gene_id = "jabberwocky"),
                 regexp = "is.element(el = gene_id, ",
                 fixed = TRUE)
    expect_error(plot_psi_heatmap(se_psi = se_psi,
                                  gene_id = "ENSG00000105974",
                                  heatmap_colors = NULL),
                 regexp = "heatmap_colors is not a character vector",
                 fixed = TRUE)
    expect_error(plot_psi_heatmap(se_psi = se_psi,
                                  gene_id = "ENSG00000105974",
                                  region_colors = NA),
                 regexp = "region_colors is not a character vector",
                 fixed = TRUE)
    expect_error(plot_psi_heatmap(se_psi = se_psi,
                                  gene_id = "ENSG00000105974",
                                  region_colors = NA),
                 regexp = "region_colors is not a character vector",
                 fixed = TRUE)
    expect_error(plot_psi_heatmap(se_psi = se_psi,
                                  gene_id = "ENSG00000105974",
                                  region_colors = (TSS = "#FF83FF")),
                 regexp = "sort(names(region_colors)) not identical to",
                 fixed = TRUE)

    # Testing if function returns the expected output
    expect_silent(
        plot_data <- plot_psi_heatmap(
            se_psi = se_psi, gene_id = "ENSG00000105974"
        )
    )
    expect_true(methods::is(plot_data, "pheatmap"))
})
