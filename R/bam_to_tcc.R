#' Prepare a TCC SummarizedExperiment object
#'
#' Prepares a TCC (Transcript Compatibility Counts) SummarizedExperiment object
#' for the given BAM files and transcript set.
#'
#' @param bam_files A named character vector containing BAM file paths.
#' @param transcript_data A named list containing transcript data returned by
#' the \code{\link{prepare_transcripts}} function.
#' @param run_mode A string specifying the mode for choosing the transcript set
#' ('strict', 'de_novo_strict', 'de_novo_loose' or 'de_novo_full').
#' @param min_read_count An integer scalar specifying the read count threshold
#' for transcripts extracted from the BAM files.
#' @param min_relative_expression A numeric scalar specifying the relative
#' expression threshold for transcripts extracted from the BAM files.
#' @param extend_spliced_transcripts An integer scalar specifying the number of
#' base pairs by which transcript starts and ends are extended for spliced read
#' compatibility search.
#' @param is_single_cell A logical scalar specifying if the BAM files contain
#' single cell data.
#' @param barcode_tag A string specifying the name of the BAM file tag
#' containing cell barcodes.
#' @param chunk_size An integer scalar specifying the chunk size for reading
#' the BAM files.
#' @param ncpu An integer scalar specifying the number of cores to use for
#' multicore parallelization.
#' @return A SummarizedExperiment object containing TCC annotation and
#' quantification data.
#' @export
bam_to_tcc <- function(bam_files,
                       transcript_data,
                       run_mode = "strict",
                       min_read_count = 1,
                       min_relative_expression = 0.1,
                       extend_spliced_transcripts = 100,
                       is_single_cell = FALSE,
                       barcode_tag = "BC",
                       chunk_size = 1000000,
                       ncpu = 1) {

    # Check arguments
    assertthat::assert_that(is.character(bam_files))
    assertthat::assert_that(length(bam_files) > 0)
    assertthat::assert_that(!is.null(names(bam_files)))
    assertthat::assert_that(all(file.exists(bam_files)))
    assertthat::assert_that(is.list(transcript_data))
    assertthat::assert_that(identical(
        names(transcript_data),
        c("tx_df","tx_granges", "tx_exon_granges_list","tx_intron_granges_list")
    ))
    assertthat::assert_that(is.data.frame(transcript_data$tx_df))
    assertthat::assert_that(class(transcript_data$tx_granges) == "GRanges")
    assertthat::assert_that(grepl(
        "GRangesList", class(transcript_data$tx_exon_granges_list)
    ))
    assertthat::assert_that(grepl(
        "GRangesList", class(transcript_data$tx_intron_granges_list)
    ))
    assertthat::assert_that(assertthat::is.string(run_mode))
    assertthat::assert_that(run_mode %in% c("strict", "de_novo_strict",
                                            "de_novo_loose", "de_novo_full"))
    assertthat::assert_that(assertthat::is.count(min_read_count))
    assertthat::assert_that(is.numeric(min_relative_expression))
    assertthat::assert_that(min_relative_expression >= 0)
    assertthat::assert_that(
        extend_spliced_transcripts == 0 ||
            assertthat::is.count(extend_spliced_transcripts)
    )
    assertthat::assert_that(assertthat::is.flag(is_single_cell))
    assertthat::assert_that(assertthat::is.string(barcode_tag))
    assertthat::assert_that(assertthat::is.count(chunk_size))
    assertthat::assert_that(assertthat::is.count(ncpu))

    # Prepare the parallelization backend
    BPPARAM <- BiocParallel::MulticoreParam(ncpu)

    # Select full-length transcripts
    tx_fl_selector <- transcript_data$tx_df$fivethree_support_level == "FL"
    if (run_mode == "strict") {
        tx_fa_selector <- transcript_data$tx_df$splicing_support_level == "AP"
        tx_selector <- tx_fl_selector & tx_fa_selector
    }
    if (run_mode == "de_novo_strict") {
        tx_fa_selector <- transcript_data$tx_df$splicing_support_level == "AP"
        tx_denovo_selector <-
            transcript_data$tx_df$splicing_support_level %in% c("EC", "NC") &
            !is.na(transcript_data$tx_df$gene_id) &
            transcript_data$tx_df$read_count >= min_read_count &
            transcript_data$tx_df$relative_expression >= min_relative_expression
        tx_selector <- tx_fl_selector & (tx_fa_selector | tx_denovo_selector)
    }
    if (run_mode == "de_novo_loose") {
        tx_fa_selector <- transcript_data$tx_df$splicing_support_level == "AP"
        tx_denovo_selector <-
            transcript_data$tx_df$splicing_support_level %in% c("EC", "NC", "DN") &
            !is.na(transcript_data$tx_df$gene_id) &
            transcript_data$tx_df$read_count >= min_read_count &
            transcript_data$tx_df$relative_expression >= min_relative_expression
        tx_selector <- tx_fl_selector & (tx_fa_selector | tx_denovo_selector)
    }
    if (run_mode == "de_novo_full") {
        warning("The de_novo_full run mode is an experimental feature that ",
                "only uses transcripts discovered by Isosceles - use it at ",
                "your own risk!")
        tx_denovo_selector <-
            transcript_data$tx_df$splicing_support_level %in% c("PC", "EC", "NC", "DN") &
            !is.na(transcript_data$tx_df$gene_id) &
            transcript_data$tx_df$read_count >= min_read_count &
            transcript_data$tx_df$relative_expression >= min_relative_expression
        tx_selector <- tx_fl_selector & tx_denovo_selector
    }

    # Prepare selected transcripts
    tx_df <- transcript_data$tx_df[tx_selector,]
    tx_granges <- transcript_data$tx_granges[tx_selector]
    tx_exon_granges_list <- transcript_data$tx_exon_granges_list[tx_selector]
    tx_intron_granges_list <- transcript_data$tx_intron_granges_list[tx_selector]

    # Prepare transcript exon data
    tx_exon_granges <- unlist(tx_exon_granges_list)
    S4Vectors::mcols(tx_exon_granges)$tx_idx <- match(
        names(tx_exon_granges), names(tx_exon_granges_list)
    )
    names(tx_exon_granges) <- NULL

    # Prepare transcript intron index data
    tx_intron_idx <- tx_intron_granges_list %>%
        unlist() %>%
        BiocGenerics::unstrand() %>%
        as.character()
    nr_unstranded_intron_positions <- unique(tx_intron_idx)
    tx_intron_idx <- split(
        match(tx_intron_idx, nr_unstranded_intron_positions),
        names(tx_intron_idx)
    )
    tx_intron_idx <- tx_intron_idx[names(tx_intron_granges_list)]
    names(tx_intron_idx) <- NULL

    # Prepare transcript start & end positions
    tx_starts <- BiocGenerics::start(tx_granges)
    tx_ends <- BiocGenerics::end(tx_granges)

    # Prepare intron start & end positions
    nr_intron_granges <- methods::as(nr_unstranded_intron_positions, "GRanges")
    nr_intron_starts <- BiocGenerics::start(nr_intron_granges)
    nr_intron_ends <- BiocGenerics::end(nr_intron_granges)

    # Process the BAM files
    bam_list <- BiocParallel::bplapply(seq_along(bam_files), function(i) {

        ## Parse the BAM file
        bam_file <- bam_files[i]
        bam_file_con <- Rsamtools::BamFile(bam_file, yieldSize = chunk_size)
        bam_tags <- character(0)
        if (is_single_cell) {
            bam_tags <- c(bam_tags, barcode_tag)
        }
        bam_param <- Rsamtools::ScanBamParam(
            what = "qname",
            flag = Rsamtools::scanBamFlag(isSupplementaryAlignment = FALSE),
            tag = bam_tags
        )
        open(bam_file_con)
        count_df <- data.frame(
            ec_id = character(0),
            count = numeric(0)
        )
        if (is_single_cell) {
            count_df <- data.frame(
                ec_id = character(0),
                cell_barcode = character(0),
                count = numeric(0)
            )
        }
        read_length <- integer(0)
        repeat {

            ### Read the chunk of the BAM file
            chunk <- GenomicAlignments::readGAlignments(bam_file_con,
                                                        use.names = TRUE,
                                                        param = bam_param)
            if (length(chunk) == 0L)
                break

            ### Process spliced reads
            chunk_spliced <- chunk[GenomicAlignments::njunc(chunk) != 0,]
            chunk_spliced_position_df <- data.frame(
                read_id = names(chunk_spliced),
                read_start = BiocGenerics::start(chunk_spliced),
                read_end = BiocGenerics::end(chunk_spliced)
            )
            if (is_single_cell) {
                chunk_spliced_position_df$cell_barcode <-
                    S4Vectors::mcols(chunk_spliced)[, barcode_tag]
            }
            chunk_introns <- GenomicAlignments::junctions(chunk_spliced) %>%
                unlist() %>%
                BiocGenerics::unstrand()
            chunk_spliced_df <- data.frame(
                read_id = names(chunk_introns),
                position = as.character(chunk_introns)
            ) %>%
                dplyr::group_by(.data$read_id) %>%
                dplyr::summarise(
                    intron_idx = list(fastmatch::fmatch(
                        .data$position, nr_unstranded_intron_positions
                    ))
                )
            all_introns_known <- sapply(chunk_spliced_df$intron_idx, function(x) {
                !any(is.na(x))
            })
            chunk_spliced_df <- chunk_spliced_df[all_introns_known,]
            chunk_spliced_df$introns <- sapply(chunk_spliced_df$intron_idx,
                                               paste0, collapse = ",")
            chunk_spliced_df <- chunk_spliced_df %>%
                dplyr::left_join(chunk_spliced_position_df,
                                 relationship = "many-to-many")
            chunk_spliced_df$read_id <- as.numeric(as.factor(chunk_spliced_df$read_id))
            chunk_intron_idx_unique <- unique(chunk_spliced_df$intron_idx)
            names(chunk_intron_idx_unique) <- sapply(chunk_intron_idx_unique, paste0,
                                                     collapse = ",")
            chunk_spliced_df$intron_idx <- NULL
            if (length(chunk_intron_idx_unique) > 0) {
                spliced_compatibility_df <- check_splicing_compatibility(
                    chunk_intron_idx_unique, tx_intron_idx
                )
                spliced_compatibility_df <- spliced_compatibility_df %>%
                    dplyr::rename(
                        introns = "subject_idx",
                        transcript_idx = "target_idx"
                    )
                spliced_compatibility_df$transcript_start <-
                    tx_starts[spliced_compatibility_df$transcript_idx]
                spliced_compatibility_df$transcript_end <-
                    tx_ends[spliced_compatibility_df$transcript_idx]
                spliced_compatibility_df$prev_intron_end <- sapply(
                    seq(nrow(spliced_compatibility_df)), function(i) {
                        introns <- spliced_compatibility_df$introns[i]
                        transcript_idx <- spliced_compatibility_df$transcript_idx[i]
                        read_intron_idx <- chunk_intron_idx_unique[[introns]]
                        transcript_intron_idx <- tx_intron_idx[[transcript_idx]]
                        position <- match(read_intron_idx[1], transcript_intron_idx)
                        prev_intron_end <- -Inf
                        if (position != 1) {
                            prev_intron_end <- nr_intron_ends[transcript_intron_idx[position - 1]]
                        }
                        return(prev_intron_end)
                    }
                )
                spliced_compatibility_df$next_intron_start <- sapply(
                    seq(nrow(spliced_compatibility_df)), function(i) {
                        introns <- spliced_compatibility_df$introns[i]
                        transcript_idx <- spliced_compatibility_df$transcript_idx[i]
                        read_intron_idx <- chunk_intron_idx_unique[[introns]]
                        transcript_intron_idx <- tx_intron_idx[[transcript_idx]]
                        position <- match(read_intron_idx[length(read_intron_idx)],
                                          transcript_intron_idx)
                        next_intron_start <- Inf
                        if (position != length(transcript_intron_idx)) {
                            next_intron_start <- nr_intron_starts[transcript_intron_idx[position + 1]]
                        }
                        return(next_intron_start)
                    }
                )
                spliced_compatibility_df <- chunk_spliced_df %>%
                    dplyr::left_join(spliced_compatibility_df,
                                     relationship = "many-to-many")
                spliced_compatibility_df <- spliced_compatibility_df %>%
                    dplyr::filter(.data$read_start >= .data$transcript_start - extend_spliced_transcripts,
                                  .data$read_end <= .data$transcript_end + extend_spliced_transcripts,
                                  .data$read_start > .data$prev_intron_end,
                                  .data$read_end < .data$next_intron_start)
                spliced_count_df <- spliced_compatibility_df %>%
                    dplyr::group_by(dplyr::across(tidyselect::any_of(c("read_id", "cell_barcode")))) %>%
                    dplyr::summarise(
                        ec_id = paste0(sort(.data$transcript_idx), collapse = ",")
                    ) %>%
                    dplyr::ungroup() %>%
                    dplyr::group_by(dplyr::across(tidyselect::any_of(c("ec_id", "cell_barcode")))) %>%
                    dplyr::summarise(count = dplyr::n()) %>%
                    dplyr::ungroup() %>%
                    dplyr::mutate(count = as.numeric(.data$count))
            } else {
                spliced_count_df <- data.frame(
                    ec_id = character(0),
                    count = numeric(0)
                )
                if (is_single_cell) {
                    spliced_count_df <- data.frame(
                        ec_id = character(0),
                        cell_barcode = character(0),
                        count = numeric(0)
                    )
                }
            }

            ### Process unspliced reads
            chunk_unspliced <- chunk[GenomicAlignments::njunc(chunk) == 0,]
            chunk_exon_granges <- methods::as(chunk_unspliced, "GRanges")
            unspliced_hits <- GenomicRanges::findOverlaps(
                chunk_exon_granges, tx_exon_granges, ignore.strand = TRUE, type = "within"
            )
            unspliced_compatibility_df <- data.frame(
                subject_idx = S4Vectors::queryHits(unspliced_hits),
                target_idx = S4Vectors::mcols(tx_exon_granges)$tx_idx[
                    S4Vectors::subjectHits(unspliced_hits)
                ]
            )
            if (is_single_cell) {
                unspliced_compatibility_df$cell_barcode <-
                    S4Vectors::mcols(chunk_unspliced)[, barcode_tag][
                        S4Vectors::queryHits(unspliced_hits)
                    ]
            }
            unspliced_count_df <- unspliced_compatibility_df %>%
                dplyr::group_by(dplyr::across(tidyselect::any_of(c("subject_idx", "cell_barcode")))) %>%
                dplyr::summarise(
                    ec_id = paste0(sort(.data$target_idx), collapse = ",")
                ) %>%
                dplyr::ungroup() %>%
                dplyr::group_by(dplyr::across(tidyselect::any_of(c("ec_id", "cell_barcode")))) %>%
                dplyr::summarise(count = dplyr::n()) %>%
                dplyr::ungroup() %>%
                dplyr::mutate(count = as.numeric(.data$count))

            ### Merge the data
            count_df <- rbind(count_df, spliced_count_df, unspliced_count_df)
            read_length <- c(read_length, GenomicAlignments::qwidth(chunk))
        }
        close(bam_file_con)

        ## Remove redundant EC rows
        count_df <- count_df %>%
            dplyr::group_by(dplyr::across(tidyselect::any_of(c("ec_id", "cell_barcode")))) %>%
            dplyr::summarise(count = sum(.data$count)) %>%
            dplyr::ungroup()
        count_df$sample_idx <- i

        ## Remove ECs matching more than one gene
        ec_gene_count <- sapply(strsplit(count_df$ec_id, ","), function(x) {
            length(unique(tx_df$gene_id[as.numeric(x)]))
        })
        count_df <- count_df[ec_gene_count == 1,]

        return(list(
            ec_count_df = count_df,
            read_length = read_length
        ))
    }, BPPARAM = BPPARAM)
    ec_count_df <- lapply(bam_list, "[[", "ec_count_df")
    ec_count_df <- do.call(rbind, ec_count_df)

    # Calculate mean read length
    mean_read_length <- lapply(bam_list, "[[", "read_length") %>%
        unlist() %>%
        mean()

    # Prepare sample IDs
    sample_ids <- names(bam_files)
    if (is_single_cell) {
        ec_count_df$sample_id <- paste0(
            names(bam_files)[ec_count_df$sample_idx], ".", ec_count_df$cell_barcode
        )
        sample_ids <- unique(ec_count_df$sample_id)
        ec_count_df$sample_idx <- fastmatch::fmatch(ec_count_df$sample_id, sample_ids)
        sample_ids <- as.character(sample_ids)
        ec_count_df$cell_barcode <- NULL
        ec_count_df$sample_id <- NULL
    }

    # Prepare non-redundant set of ECs
    ec_ids <- unique(ec_count_df$ec_id)
    ec_gene_ids <- sapply(strsplit(ec_ids, ","), function(x) {
        unique(tx_df$gene_id[as.numeric(x)])
    })
    ec_gene_names <- sapply(strsplit(ec_ids, ","), function(x) {
        unique(tx_df$gene_name[as.numeric(x)])
    })
    ec_order_idx <- order(ec_gene_ids)
    ec_ids <- ec_ids[ec_order_idx]
    ec_gene_ids <- ec_gene_ids[ec_order_idx]
    ec_gene_names <- ec_gene_names[ec_order_idx]

    # Prepare the EC count matrix
    ec_counts <- Matrix::sparseMatrix(
        i = match(ec_count_df$ec_id, ec_ids),
        j = ec_count_df$sample_idx,
        x = ec_count_df$count,
        dims = c(length(ec_ids), length(sample_ids))
    )
    if (!is_single_cell) {
        ec_counts <- as.matrix(ec_counts)
    }
    rownames(ec_counts) <- ec_ids
    colnames(ec_counts) <- sample_ids

    # Prepare the EC metadata
    ec_metadata_df <- data.frame(
        ec_id = ec_ids,
        gene_id = ec_gene_ids,
        gene_name = ec_gene_names
    )

    # Prepare the EC compatibility matrix
    ec_compatibility_df <- data.frame(
        ec_idx = seq_along(ec_ids),
        ec_id = ec_ids
    )
    ec_compatibility_df$transcript_idx <- strsplit(ec_compatibility_df$ec_id, ",")
    ec_compatibility_df <- ec_compatibility_df %>%
        tidyr::unchop("transcript_idx") %>%
        dplyr::mutate(transcript_idx = as.integer(.data$transcript_idx))
    ec_compatibility <- Matrix::sparseMatrix(
        i = ec_compatibility_df$ec_idx,
        j = ec_compatibility_df$transcript_idx,
        x = 1,
        dims = c(length(ec_ids), nrow(tx_df))
    )
    rownames(ec_compatibility) <- ec_ids
    colnames(ec_compatibility) <- tx_df$transcript_id

    # Prepare the SummarizedExperiment object
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = ec_counts),
        rowData = ec_metadata_df
    )
    SummarizedExperiment::assay(se, "tpm") <- calculate_tpm(se)
    SummarizedExperiment::assay(se, "relative_expression") <-
        calculate_relative_expression(se)
    S4Vectors::metadata(se)$compatibility_matrix <- ec_compatibility
    S4Vectors::metadata(se)$transcript_df <- tx_df
    S4Vectors::metadata(se)$transcript_exon_granges_list <- tx_exon_granges_list
    S4Vectors::metadata(se)$mean_read_length <- mean_read_length

    return(se)
}
