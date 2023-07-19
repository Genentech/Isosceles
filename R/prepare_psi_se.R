#' Prepare a PSI SummarizedExperiment object
#'
#' Prepares a PSI (Percent Spliced In) SummarizedExperiment object
#' for the given transcript-level SummarizedExperiment object
#'
#' @param se A transcript-level SummarizedExperiment object returned by the
#' \code{\link{prepare_transcript_se}} function
#' @param ncpu An integer scalar specifying the number of cores to use for
#' multicore parallelization
#' @return A SummarizedExperiment object containing PSI annotation and
#' quantification data
#' @export
prepare_psi_se <- function(se,
                           ncpu = 1) {

    # Check arguments
    assertthat::assert_that(methods::is(se, "SummarizedExperiment"))
    assertthat::assert_that(is.element(
        "relative_expression", SummarizedExperiment::assayNames(se)
    ))
    assertthat::assert_that(is.element(
        "gene_id", colnames(SummarizedExperiment::rowData(se))
    ))
    assertthat::assert_that(is.element(
        "transcript_id", colnames(SummarizedExperiment::rowData(se))
    ))
    assertthat::assert_that(is.element(
        "position", colnames(SummarizedExperiment::rowData(se))
    ))
    assertthat::assert_that(
        grepl("GRangesList", class(SummarizedExperiment::rowRanges(se)))
    )
    assertthat::assert_that(assertthat::is.count(ncpu))

    # Fix the "No visible binding for global variable '.'" issue
    . = NULL

    # Prepare the parallelization backend
    BPPARAM <- BiocParallel::MulticoreParam(ncpu)

    # Remove unspliced transcripts from the SummarizedExperiment
    se <- se[!is.na(SummarizedExperiment::rowData(se)$intron_positions),]

    # Prepare input data
    is_bulk_rnaseq <- is.matrix(
        SummarizedExperiment::assay(se, "relative_expression")
    )
    tx_data_list <- SummarizedExperiment::rowData(se) %>%
        as.data.frame() %>%
        dplyr::mutate(transcript_idx = seq_along(.data$transcript_id)) %>%
        split(.$gene_id)

    # Prepare non-overlapping transcript regions
    regions_data <- BiocParallel::bplapply(tx_data_list, function(tx_data_df) {

        # Prepare transcript data
        tx_granges <- methods::as(tx_data_df$position, "GRanges")
        exon_granges_list <-
            SummarizedExperiment::rowRanges(se)[tx_data_df$transcript_id]
        intron_granges_list <- IRanges::psetdiff(tx_granges, exon_granges_list)
        seq_name <- unique(as.character(GenomeInfoDb::seqnames(tx_granges)))
        strand <- unique(as.character(BiocGenerics::strand(tx_granges)))
        gene_id <- unique(tx_data_df$gene_id)

        # Prepare start positions of exons/introns/etc.
        intron_starts <- BiocGenerics::start(intron_granges_list)
        intron_starts <- unname(unlist(intron_starts))
        exon_starts <- BiocGenerics::end(intron_granges_list) + 1
        exon_starts <- unname(unlist(exon_starts))
        tss <- tx_data_df$transcript_id %>%
            strsplit(":") %>%
            sapply("[", 2) %>%
            (function(x) {gsub("^s", "", x)}) %>%
            as.numeric()
        tes <- tx_data_df$transcript_id %>%
            strsplit(":") %>%
            sapply("[", 3) %>%
            (function(x) {gsub("^e", "", x)}) %>%
            as.numeric()
        if (strand == "-") {
            tss_copy <- tss
            tss <- tes
            tes <- tss_copy
        }

        # Prepare non-redundant region data
        regions_df <- data.frame(
            position =  c(intron_starts, exon_starts),
            type = c(rep("IN", length(intron_starts)),
                     rep("EX", length(exon_starts)))
        )
        regions_df <- regions_df %>%
            dplyr::distinct() %>%
            dplyr::arrange(.data$position)
        region_starts <- utils::head(regions_df$position,
                                     length(regions_df$position) - 1)
        region_ends <- utils::tail(regions_df$position,
                                   length(regions_df$position) - 1) - 1
        region_granges <- methods::as(
            glue::glue("{seq_name}:{region_starts}-{region_ends}:{strand}"),
            "GRanges"
        )
        S4Vectors::mcols(region_granges)$type <-
            sapply(seq(nrow(regions_df) - 1), function(i) {
                type_start <- regions_df$type[i]
                type_end <- regions_df$type[i + 1]
                if (type_start == "EX" && type_end == "IN" && strand == "+")
                    return("CE")
                else if (type_start == "IN" && type_end == "EX" && strand == "-")
                    return("CE")
                else if (type_start == "IN" && type_end == "EX" && strand == "+")
                    return("RI")
                else if (type_start == "EX" && type_end == "IN" && strand == "-")
                    return("RI")
                else if (type_start == "IN" && type_end == "IN" && strand == "+")
                    return("A5")
                else if (type_start == "EX" && type_end == "EX" && strand == "-")
                    return("A5")
                else if (type_start == "EX" && type_end == "EX" && strand == "+")
                    return("A3")
                else if (type_start == "IN" && type_end == "IN" && strand == "-")
                    return("A3")
            })

        # Map regions to transcripts they overlap
        exon_transcript_idx <- rep(
            tx_data_df$transcript_idx,
            S4Vectors::elementNROWS(exon_granges_list)
        )
        region_overlaps <- IRanges::findOverlaps(
            region_granges, unlist(exon_granges_list),
            type = "within"
        ) %>% as.data.frame()
        region_transcripts <- split(
            exon_transcript_idx[region_overlaps$subjectHits],
            as.character(region_granges)[region_overlaps$queryHits]
        )
        region_selector <-
            as.character(region_granges) %in% names(region_transcripts)
        region_granges <- region_granges[region_selector]

        # Prepare TSS/TES data
        tss_granges <- methods::as(
            glue::glue("{seq_name}:{tss}-{tss}:{strand}"),
            "GRanges"
        )
        S4Vectors::mcols(tss_granges)$type <- "TSS"
        tes_granges <- methods::as(
            glue::glue("{seq_name}:{tes}-{tes}:{strand}"),
            "GRanges"
        )
        S4Vectors::mcols(tes_granges)$type <- "TES"

        # Map TSS/TES to transcripts
        tss_transcripts <- split(
            tx_data_df$transcript_idx,
            as.character(tss_granges)
        )
        tes_transcripts <- split(
            tx_data_df$transcript_idx,
            as.character(tes_granges)
        )

        # Add TSS/TES data to regions
        region_granges <- c(region_granges,
                            BiocGenerics::unique(tss_granges),
                            BiocGenerics::unique(tes_granges))
        region_granges <- GenomicRanges::sort(
            region_granges, decreasing = (strand == "-")
        )
        region_transcripts <- c(region_transcripts, tss_transcripts, tes_transcripts)
        region_transcripts <- region_transcripts[as.character(region_granges)]

        # Prepare summarized region data
        regions_tbl <- tibble::tibble(
            gene_id = gene_id,
            range = as.character(region_granges),
            type = S4Vectors::mcols(region_granges)$type,
            transcripts = unname(region_transcripts)
        )
        return(regions_tbl)
    }, BPPARAM = BPPARAM)
    regions_data <- do.call(rbind, regions_data)

    # Calculate PSI values
    psi_matrix <- scuttle::sumCountsAcrossFeatures(
        SummarizedExperiment::assay(se, "relative_expression"),
        regions_data$transcripts
    )
    if (!is_bulk_rnaseq) {
        psi_matrix <- methods::as(psi_matrix, "dgCMatrix")
    }
    rownames(psi_matrix) <- paste0(
        regions_data$gene_id, ":", regions_data$range, ":", regions_data$type
    )

    # Create the SummarizedExperiment object
    se_psi <- SummarizedExperiment::SummarizedExperiment(
        assays = list(psi = psi_matrix),
        rowRanges = methods::as(regions_data$range, "GRanges")
    )
    row_data_columns <- c("gene_id", "type")
    SummarizedExperiment::rowData(se_psi) <- regions_data[, row_data_columns]

    # Remove regions with no expression
    region_selector <-
        Matrix::rowSums(SummarizedExperiment::assay(se_psi, "psi")) != 0
    se_psi <- se_psi[region_selector,]

    return(se_psi)
}
