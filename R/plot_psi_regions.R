#' Plot PSI regions
#'
#' Creates a plot showing PSI regions and transcript structures for the given
#' gene. Individual transcript structures are colored by their relative
#' expression, calculated from the overall TPM values and expressed in
#' percentages. For better visualization, introns can be shrinked using the
#' max_intron_length argument
#'
#' @param se_psi A PSI SummarizedExperiment object returned by the
#' \code{\link{prepare_psi_se}} function
#' @param se_transcript A transcript-level SummarizedExperiment object returned
#' by the \code{\link{prepare_transcript_se}} function
#' @param gene_id A string containing the identifier of the gene to plot
#' @param max_transcripts An integer scalar specifying the maximum number of
#' transcripts with the highest relative expression to plot
#' @param max_intron_length An integer scalar specifying the maximum intron
#' length after shrinking. If set to NULL, no shrinking is performed
#' @param region_colors A named character vector of colors for the PSI region
#' types
#' @return A plot object
#' @export
plot_psi_regions <- function(se_psi,
                             se_transcript,
                             gene_id,
                             max_transcripts = Inf,
                             max_intron_length = NULL,
                             region_colors = NULL) {

    # Check arguments
    assertthat::assert_that(methods::is(se_psi, "SummarizedExperiment"))
    assertthat::assert_that(is.element(
        "psi", SummarizedExperiment::assayNames(se_psi)
    ))
    assertthat::assert_that(is.element(
        "gene_id", colnames(SummarizedExperiment::rowData(se_psi))
    ))
    assertthat::assert_that(methods::is(se_transcript, "SummarizedExperiment"))
    assertthat::assert_that(is.element(
        "tpm", SummarizedExperiment::assayNames(se_transcript)
    ))
    assertthat::assert_that(is.element(
        "gene_id", colnames(SummarizedExperiment::rowData(se_transcript))
    ))
    assertthat::assert_that(assertthat::is.string(gene_id))
    assertthat::assert_that(is.element(
        gene_id, SummarizedExperiment::rowData(se_psi)$gene_id
    ))
    assertthat::assert_that(is.element(
        gene_id, SummarizedExperiment::rowData(se_transcript)$gene_id
    ))
    assertthat::assert_that(assertthat::is.count(max_transcripts))
    if (!is.null(max_intron_length)) {
        assertthat::assert_that(assertthat::is.count(max_intron_length))
    }
    if (!is.null(region_colors)) {
        assertthat::assert_that(is.character(region_colors))
        assertthat::assert_that(identical(
            sort(names(region_colors)),
            c("A3", "A5", "CE", "RI", "TES", "TSS")
        ))
    }

    # Filter the SummarizedExperiment objects
    se_psi <- se_psi[
        SummarizedExperiment::rowData(se_psi)$gene_id == gene_id,
    ]
    se_transcript <- se_transcript[
        SummarizedExperiment::rowData(se_transcript)$gene_id == gene_id,
    ]

    # Calculate the overall relative expression of transcripts
    tx_tpm_sums <- se_transcript %>%
        SummarizedExperiment::assay("tpm") %>%
        Matrix::rowSums()
    tx_rel_exprs <- tx_tpm_sums / sum(tx_tpm_sums) * 100
    tx_rel_exprs <- sort(tx_rel_exprs, decreasing = TRUE)
    tx_rel_exprs <- utils::head(tx_rel_exprs, n = max_transcripts)

    # Prepare the transcript GenomicRangesList object for visualization
    tx_granges_list <- SummarizedExperiment::rowRanges(se_transcript)
    tx_granges_list <- tx_granges_list[names(tx_rel_exprs)]
    tx_granges_list <- unlist(tx_granges_list)
    S4Vectors::mcols(tx_granges_list)$plot_type <- "exon"
    S4Vectors::mcols(tx_granges_list)$rel_expr <-
        tx_rel_exprs[names(tx_granges_list)]
    tx_granges_list <- GenomicRanges::split(
        tx_granges_list, names(tx_granges_list)
    )

    # Prepare the PSI regions GenomicRanges objects for visualization
    region_granges <- SummarizedExperiment::rowRanges(se_psi)
    bin_region_types <- c("TSS", "TES")
    other_region_types <- c("CE", "RI", "A5", "A3")
    region_granges_bin <- region_granges[
        S4Vectors::mcols(region_granges)$type %in% bin_region_types
    ]
    region_granges_other <- region_granges[
        S4Vectors::mcols(region_granges)$type %in% other_region_types
    ]
    S4Vectors::mcols(region_granges_bin)$plot_type <- "exon"
    S4Vectors::mcols(region_granges_other)$plot_type <- "exon"

    # Shrink introns if requested
    if (!is.null(max_intron_length)) {
        # Prepare the shrinkage function
        region_granges_list <- GenomicRanges::GRangesList(
            PSI = SummarizedExperiment::rowRanges(se_psi)
        )
        shrink_func <- biovizBase::shrinkageFun(
            IRanges::gaps(unlist(c(tx_granges_list, region_granges_list))),
            max.gap = max_intron_length
        )

        # Shrink the transcript GenomicRangesList object
        tx_granges_list <- S4Vectors::endoapply(tx_granges_list, function(gr) {
            gr <- shrink_func(gr)
            S4Vectors::mcols(gr)$.ori <- NULL
            return(gr)
        })

        # Shrink the PSI regions GenomicRanges objects
        region_granges_bin <- GenomicRanges::GRanges(
            shrink_func(region_granges_bin)
        )
        S4Vectors::mcols(region_granges_bin)$.ori <- NULL
        region_granges_other <- GenomicRanges::GRanges(
            shrink_func(region_granges_other)
        )
        S4Vectors::mcols(region_granges_other)$.ori <- NULL
    }

    # Create the plot
    ## TSS / TES sites
    bin_region_plot <- ggbio::autoplot(
        GenomicRanges::GRanges(region_granges_bin),
        ggplot2::aes(type = .data$plot_type,
                     col = .data$type,
                     fill = .data$type)
    )
    if (!is.null(region_colors)) {
        bin_region_plot <- bin_region_plot +
            ggplot2::scale_fill_manual(
                values = region_colors[bin_region_types]
            ) +
            ggplot2::scale_colour_manual(
                values = region_colors[bin_region_types]
            )
    }
    bin_region_plot <- bin_region_plot +
        ggplot2::theme(
            legend.key.size = ggplot2::unit(0.4, 'cm'),
            legend.title = ggplot2::element_blank()
        )
    ## Other PSI regions
    other_region_plot <- ggbio::autoplot(
        GenomicRanges::GRanges(region_granges_other),
        ggplot2::aes(type = .data$plot_type,
                     col = .data$type,
                     fill = .data$type)
    )
    if (!is.null(region_colors)) {
        other_region_plot <- other_region_plot +
            ggplot2::scale_fill_manual(
                values = region_colors[other_region_types]
            ) +
            ggplot2::scale_colour_manual(
                values = region_colors[other_region_types]
            )
    }
    other_region_plot <- other_region_plot +
        ggplot2::theme(
            legend.key.size = ggplot2::unit(0.3, 'cm'),
            legend.title = ggplot2::element_blank()
        )
    ## Transcript structures
    tx_plot <- ggbio::autoplot(
        tx_granges_list,
        ggplot2::aes(type = .data$plot_type,
                     col = .data$rel_expr,
                     fill = .data$rel_expr)
    )
    ## Combined plot
    psi_plot <- ggbio::tracks(
        `TSS/TES\nsites` = bin_region_plot,
        tx_plot,
        `PSI\nregions` = other_region_plot,
        heights = c(0.1, 0.8, 0.1),
        label.text.cex = 0.7
    )

    return(psi_plot)
}
