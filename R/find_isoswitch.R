#' Find isoform switching events
#'
#' Identifies isoform switching events by comparing every pair of cell groups
#' using the \code{\link{findMarkers}} function from the scran package and
#' searching for transcripts of the same gene showing statistically significant
#' differences in opposite directions.
#'
#' @param se A transcript-level SummarizedExperiment object returned by the
#' \code{\link{tcc_to_transcript}} function. The object must contain normalized
#' data stored in the 'logcounts' assay, which can be prepared using functions
#' from the scuttle package.
#' @param cell_labels A vector or a factor containing cell labels acting as a
#' grouping variable.
#' @param min_fdr A numeric scalar specifying the FDR threshold for filtering
#' the results.
#' @param ncpu An integer scalar specifying the number of cores to use for
#' multicore parallelization.
#' @return A data frame containing the following columns:
#' \describe{
#'   \item{transcript_id}{Isosceles transcript ID}
#'   \item{compatible_tx}{comma-separated list of annotated transcript IDs compatible with the Isosceles transcript}
#'   \item{gene_id}{gene ID}
#'   \item{gene_name}{gene symbol}
#'   \item{pvalue}{p-value from the Wilcoxon test performed by the findMarkers function}
#'   \item{fdr}{false discovery rate (FDR) value from the Wilcoxon test performed by the findMarkers function}
#'   \item{auc}{area under the curve (AUC) value from the Wilcoxon test performed by the findMarkers function}
#'   \item{group_1}{label of the cell group in which the transcript is upregulated}
#'   \item{group_2}{label of the cell group compared to which the transcript is upregulated}
#'   \item{contrast}{label of the compared cell group pair}
#' }
#' @export
find_isoswitch <- function(se,
                           cell_labels,
                           min_fdr = 0.05,
                           ncpu = 1) {

    # Check arguments
    assertthat::assert_that(methods::is(se, "SummarizedExperiment"))
    assertthat::assert_that(is.element(
        "logcounts", SummarizedExperiment::assayNames(se)
    ))
    assertthat::assert_that(
        grepl("GRangesList", class(SummarizedExperiment::rowRanges(se)))
    )
    assertthat::assert_that(length(cell_labels) == ncol(se))
    assertthat::assert_that(assertthat::is.number(min_fdr))
    assertthat::assert_that(min_fdr < 1)
    assertthat::assert_that(min_fdr > 0)
    assertthat::assert_that(assertthat::is.count(ncpu))

    # Prepare the parallelization backend
    BPPARAM <- BiocParallel::MulticoreParam(ncpu)

    # Identify marker transcripts between each pair of cell groups
    group_contrasts <- utils::combn(
        levels(as.factor(cell_labels)), 2, simplify = FALSE
    )
    marker_results_list <- BiocParallel::bplapply(group_contrasts,
                                                  function(group_contrast) {
        scran::findMarkers(se,
                           groups = cell_labels,
                           restrict = group_contrast,
                           test.type = "wilcox",
                           pval.type = "any",
                           direction = "up",
                           row.data = rowData(se))
    }, BPPARAM = BPPARAM)
    names(marker_results_list) <- sapply(group_contrasts, paste0,
                                         collapse = "__")

    # Process the marker transcript data
    marker_df_list <- lapply(names(marker_results_list),
                             function(contrast_name) {
        marker_results <- marker_results_list[[contrast_name]]
        group_1_df <- marker_results[[1]] %>%
            as.data.frame() %>%
            dplyr::transmute(
                transcript_id = .data$transcript_id,
                compatible_tx = .data$compatible_tx,
                gene_id = .data$gene_id,
                gene_name = .data$gene_name,
                pvalue = .data$p.value,
                fdr = .data$FDR,
                auc = .data$summary.AUC,
                group_1 = names(marker_results)[1],
                group_2 = names(marker_results)[2],
                contrast = contrast_name
            )
        rownames(group_1_df) <- NULL
        group_2_df <- marker_results[[2]] %>%
            as.data.frame() %>%
            dplyr::transmute(
                transcript_id = .data$transcript_id,
                compatible_tx = .data$compatible_tx,
                gene_id = .data$gene_id,
                gene_name = .data$gene_name,
                pvalue = .data$p.value,
                fdr = .data$FDR,
                auc = .data$summary.AUC,
                group_1 = names(marker_results)[2],
                group_2 = names(marker_results)[1],
                contrast = contrast_name
            )
        rownames(group_2_df) <- NULL
        return(rbind(group_1_df, group_2_df))
    })
    marker_df <- do.call(rbind, marker_df_list)

    # Detect isoform switching events
    isoswitch_df <- marker_df %>%
        dplyr::filter(.data$fdr <= min_fdr) %>%
        dplyr::group_by(.data$gene_id, .data$contrast) %>%
        dplyr::mutate(
            n_groups = length(unique(.data$group_1)),
            n_transcripts = length(unique(.data$transcript_id))
        ) %>%
        dplyr::ungroup() %>%
        dplyr::filter(.data$n_groups > 1, .data$n_transcripts > 1) %>%
        dplyr::select(-"n_groups", -"n_transcripts") %>%
        dplyr::arrange(.data$contrast, .data$gene_id) %>%
        as.data.frame()

    return(isoswitch_df)
}
