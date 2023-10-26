#' Calculate PSI count to mean permuted PSI count ratio matrix
#'
#' Calculates PSI count to mean permuted PSI count ratio matrix for pseudotime
#' window data. This function is designed for preparing data to be visualized
#' as a heatmap, and might take a long time to run - see the vignettes for an
#' example.
#'
#' @param se_tcc A TCC SummarizedExperiment object returned by the
#' \code{\link{bam_to_tcc}} function.
#' @param pseudotime_matrix A numeric matrix containing the pseudotime values
#' for each cell (rows) in different trajectories (columns). Cells not
#' belonging to given trajectory should be denoted using NA values.
#' @param psi_events A character vector specifying the PSI events to calculate
#' the ratios for.
#' @param window_sizes A named integer vector specifying the window size for
#' each trajectory.
#' @param window_steps A named integer vector specifying the window step for
#' each trajectory.
#' @param trim A numeric scalar specifying the fraction (0 to 0.5) of cells
#' to be trimmed from each end of the pseudotime spectrum for each trajectory.
#' @param ncpu An integer scalar specifying the number of cores to use for
#' multicore parallelization.
#' @param n_perm An integer scalar specifying the number of PSI count
#' permutations to calculate.
#' @return A numeric matrix containing the PSI count to mean permuted PSI count
#' ratio values.
#' @export
calculate_psi_ratio_matrix <- function(se_tcc,
                                       pseudotime_matrix,
                                       psi_events,
                                       window_sizes,
                                       window_steps,
                                       trim = 0,
                                       n_perm = 100,
                                       ncpu = 1) {

    # Check arguments
    assertthat::assert_that(methods::is(se_tcc, "SummarizedExperiment"))
    assertthat::assert_that(is.element(
        "counts", SummarizedExperiment::assayNames(se_tcc)
    ))
    assertthat::assert_that(is.element(
        "gene_id", colnames(SummarizedExperiment::rowData(se_tcc))
    ))
    assertthat::assert_that(is.matrix(pseudotime_matrix))
    assertthat::assert_that(identical(nrow(pseudotime_matrix), ncol(se_tcc)))
    assertthat::assert_that(is.character(psi_events))
    assertthat::assert_that(length(psi_events) > 1)
    assertthat::assert_that(is.numeric(window_sizes))
    assertthat::assert_that(identical(
        sort(names(window_sizes)), sort(colnames(pseudotime_matrix))
    ))
    assertthat::assert_that(is.numeric(window_steps))
    assertthat::assert_that(identical(
        sort(names(window_steps)), sort(colnames(pseudotime_matrix))
    ))
    assertthat::assert_that(assertthat::is.number(trim))
    assertthat::assert_that(trim >= 0)
    assertthat::assert_that(trim < 0.5)
    assertthat::assert_that(assertthat::is.count(n_perm))
    assertthat::assert_that(assertthat::is.count(ncpu))

    # Prepare gene data
    se_gene <- tcc_to_gene(se_tcc = se_tcc)
    psi_gene_ids <- sapply(strsplit(psi_events, ":"), "[", 1)
    gene_ids <- unique(psi_gene_ids)

    # Filter the TCC SE object
    tcc_row_selector <- SummarizedExperiment::rowData(se_tcc)$gene_id %in% gene_ids
    tcc_tx_selector <- S4Vectors::metadata(se_tcc)$transcript_df$gene_id %in% gene_ids
    filtered_se_tcc <- se_tcc[tcc_row_selector,]
    S4Vectors::metadata(filtered_se_tcc)$compatibility_matrix <-
        S4Vectors::metadata(filtered_se_tcc)$compatibility_matrix[tcc_row_selector, tcc_tx_selector]
    S4Vectors::metadata(filtered_se_tcc)$transcript_df <-
        S4Vectors::metadata(filtered_se_tcc)$transcript_df[tcc_tx_selector,]
    S4Vectors::metadata(filtered_se_tcc)$transcript_exon_granges_list <-
        S4Vectors::metadata(filtered_se_tcc)$transcript_exon_granges_list[tcc_tx_selector,]

    # Prepare pseudotime windows SE objects
    se_window_tcc_list <- lapply(colnames(pseudotime_matrix), function(traj_name) {
        se_window_tcc <- pseudotime_tcc(
            se_tcc = filtered_se_tcc,
            pseudotime = pseudotime_matrix[, traj_name],
            trim = trim,
            window_size = window_sizes[traj_name],
            window_step = window_steps[traj_name]
        )
    })
    names(se_window_tcc_list) <- colnames(pseudotime_matrix)
    se_window_gene_list <- lapply(se_window_tcc_list, function(se_window_tcc) {
        se_window_gene <- tcc_to_gene(
            se_tcc = se_window_tcc
        )
    })
    se_window_transcript_list <- lapply(se_window_tcc_list, function(se_window_tcc) {
        se_window_transcript <- tcc_to_transcript(
            se_tcc = se_window_tcc,
            use_length_normalization = FALSE, ncpu = ncpu
        )
    })
    se_window_psi_list <- lapply(se_window_transcript_list, function(se_window_transcript) {
        se_window_psi <- transcript_to_psi(
            se = se_window_transcript,
            ncpu = ncpu
        )
    })

    # Prepare PSI count data
    psi_counts_list <- lapply(colnames(pseudotime_matrix), function(traj_name) {
        se_window_psi <- se_window_psi_list[[traj_name]]
        se_window_gene <- se_window_gene_list[[traj_name]]
        psi_values <- as.matrix(SummarizedExperiment::assay(se_window_psi, "psi"))
        psi_values <- sapply(colnames(psi_values), function(window_id) {
            psi_vector <- psi_values[, window_id]
            psi_vector <- psi_vector[psi_events]
            psi_vector[is.na(psi_vector)] <- 0
            names(psi_vector) <- psi_events
            return(psi_vector)
        })
        gene_counts <- as.matrix(
            SummarizedExperiment::assay(se_window_gene, "counts")[psi_gene_ids,]
        )
        psi_counts <- psi_values * gene_counts
    })
    names(psi_counts_list) <- colnames(pseudotime_matrix)

    # Prepare permuted PSI count data
    perm_psi_counts_list <- lapply(gene_ids, function(gene_id) {
        lapply(seq(n_perm), function(i) {
            perm_psi_counts_list_gene <- lapply(colnames(pseudotime_matrix), function(traj_name) {
                pseudotime <- pseudotime_matrix[, traj_name]
                names(pseudotime) <- colnames(se_tcc)
                pseudotime <- pseudotime[!is.na(pseudotime)]
                se_tcc_traj <- se_tcc[, names(pseudotime)]
                se_gene_traj <- se_gene[, names(pseudotime)]
                se_window_gene <- se_window_gene_list[[traj_name]]
                tcc_row_selector <- SummarizedExperiment::rowData(se_tcc_traj)$gene_id == gene_id
                tcc_tx_selector <- S4Vectors::metadata(se_tcc_traj)$transcript_df$gene_id == gene_id
                perm_se_tcc <- se_tcc_traj[tcc_row_selector,]
                S4Vectors::metadata(perm_se_tcc)$compatibility_matrix <-
                    S4Vectors::metadata(perm_se_tcc)$compatibility_matrix[tcc_row_selector, tcc_tx_selector]
                S4Vectors::metadata(perm_se_tcc)$transcript_df <-
                    S4Vectors::metadata(perm_se_tcc)$transcript_df[tcc_tx_selector,]
                S4Vectors::metadata(perm_se_tcc)$transcript_exon_granges_list <-
                    S4Vectors::metadata(perm_se_tcc)$transcript_exon_granges_list[tcc_tx_selector,]
                non_zero_indices <- which(
                    SummarizedExperiment::assay(se_gene_traj[gene_id,], "counts")[1,] != 0
                )
                if (length(non_zero_indices) > 1) {
                    cell_indices <- seq(ncol(perm_se_tcc))
                    cell_indices[non_zero_indices] <- sample(cell_indices[non_zero_indices])
                    perm_se_tcc <- perm_se_tcc[, cell_indices]
                }
                perm_se_window_tcc <- pseudotime_tcc(
                    se_tcc = perm_se_tcc,
                    pseudotime = pseudotime,
                    trim = trim,
                    window_size = window_sizes[traj_name],
                    window_step = window_steps[traj_name]
                )
                perm_se_window_transcript <- tcc_to_transcript(
                    se_tcc = perm_se_window_tcc,
                    use_length_normalization = FALSE, ncpu = ncpu
                )
                perm_se_window_psi <- transcript_to_psi(
                    perm_se_window_transcript,
                    ncpu = ncpu
                )
                perm_psi_values <- SummarizedExperiment::assay(perm_se_window_psi, "psi")
                orig_gene_counts <- SummarizedExperiment::assay(se_window_gene, "counts")[
                    SummarizedExperiment::rowData(perm_se_window_psi)$gene_id,
                ]
                return(perm_psi_values * orig_gene_counts)
            })
            names(perm_psi_counts_list_gene) <- colnames(pseudotime_matrix)
            return(perm_psi_counts_list_gene)
        })
    })
    names(perm_psi_counts_list) <- gene_ids

    # Prepare PSI count matrices
    psi_counts <- do.call(cbind, lapply(names(psi_counts_list), function(traj_name) {
        psi_counts <- psi_counts_list[[traj_name]]
        colnames(psi_counts) <- paste0(traj_name, ".", colnames(psi_counts))
        return(as.matrix(psi_counts))
    }))
    avg_perm_psi_counts_list <- lapply(colnames(pseudotime_matrix), function(traj_name) {
        avg_perm_matrix <- t(sapply(psi_events, function(psi_event) {
            gene_id <- strsplit(psi_event, ":")[[1]][1]
            perm_psi_event_counts <- sapply(seq(n_perm), function(i) {
                perm_matrix <- perm_psi_counts_list[[gene_id]][[i]][[traj_name]]
                if (psi_event %in% rownames(perm_matrix)) {
                    return(perm_matrix[psi_event,])
                } else {
                    return(stats::setNames(rep(0, ncol(perm_matrix)), colnames(perm_matrix)))
                }
            })
            avg_psi_event_counts <- apply(perm_psi_event_counts, 1, mean)
            names(avg_psi_event_counts) <- paste0(
                traj_name, ".", rownames(perm_psi_event_counts)
            )
            return(avg_psi_event_counts)
        }))
    })
    names(avg_perm_psi_counts_list) <- colnames(pseudotime_matrix)
    avg_perm_psi_counts <- do.call(cbind, avg_perm_psi_counts_list)

    # Prepare the PSI count ratio matrix
    psi_ratio_matrix <- psi_counts / avg_perm_psi_counts

    return(psi_ratio_matrix)
}
