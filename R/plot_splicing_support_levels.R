#' Plot spliced read counts by splicing support level
#'
#' Creates a plot showing counts of spliced reads (extracted from the BAM files
#' using the \code{\link{bam_to_read_structures}} function and processed by the
#' \code{\link{prepare_transcripts}} function) by their splicing support level.
#'
#' @param transcript_data A named list containing transcript data returned by
#' the \code{\link{prepare_transcripts}} function.
#' @return A plot object.
#' @export
plot_splicing_support_levels <- function(transcript_data) {

    # Check arguments
    assertthat::assert_that(is.list(transcript_data))
    assertthat::assert_that(identical(
        names(transcript_data),
        c("tx_df","tx_granges", "tx_exon_granges_list","tx_intron_granges_list")
    ))
    assertthat::assert_that(is.data.frame(transcript_data$tx_df))
    assertthat::assert_that(
        !all(is.na(transcript_data$tx_df$read_count)),
        msg = "transcript_data contains reference transcripts only"
    )
    assertthat::assert_that(class(transcript_data$tx_granges) == "GRanges")
    assertthat::assert_that(grepl(
        "GRangesList", class(transcript_data$tx_exon_granges_list)
    ))
    assertthat::assert_that(grepl(
        "GRangesList", class(transcript_data$tx_intron_granges_list)
    ))

    # Prepare plot data
    tx_df <- transcript_data$tx_df
    tx_df <- tx_df[!is.na(tx_df$read_count),]
    plot_data <- tapply(
        tx_df$read_count, tx_df$splicing_support_level, sum, na.rm = TRUE
    )
    plot_data <- tibble::enframe(plot_data)
    colnames(plot_data) <- c("class", "count")
    plot_data$class <- forcats::fct_rev(factor(
        plot_data$class, levels = c("PC", "EC", "NC", "DN", "AF", "AS", "AX")
    ))

    # Create the plot
    bar_plot <- ggplot2::ggplot(plot_data,
                                mapping = ggplot2::aes(x = .data$count,
                                                       y = .data$class)) +
        ggplot2::geom_col() +
        ggplot2::labs(
            x = "Read count",
            y = "Splicing support level"
        ) +
        ggplot2::theme_bw()

    return(bar_plot)
}
