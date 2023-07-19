#' Read structure extraction from BAM files
#'
#' Extract non-redundant read structures from one or multiple BAM files
#'
#' @param bam_files A character vector containing BAM file paths
#' @param chunk_size An integer scalar specifying the chunk size for reading
#' the BAM files
#' @param ncpu An integer scalar specifying the number of cores to use for
#' multicore parallelization
#' @return A data frame containing non-redundant read structure data obtained
#' from the BAM files
#' @export
extract_read_structures <- function(bam_files,
                                    chunk_size = 1000000,
                                    ncpu = 1) {

    # Check arguments
    assertthat::assert_that(is.character(bam_files))
    assertthat::assert_that(length(bam_files) > 0)
    assertthat::assert_that(all(file.exists(bam_files)))
    assertthat::assert_that(assertthat::is.count(chunk_size))
    assertthat::assert_that(assertthat::is.count(ncpu))

    # Prepare the parallelization backend
    BPPARAM <- BiocParallel::MulticoreParam(ncpu)

    # Helper functions
    merge_list <- function(x) {
        x %>%
            unlist %>%
            list()
    }

    # Get read structures from the BAM files
    read_summary <- BiocParallel::bplapply(bam_files, function(bam_file) {

        ## Parse the BAM file
        bam_read_summary <- NULL
        bam_file_con <- Rsamtools::BamFile(bam_file, yieldSize = chunk_size)
        bam_param <- Rsamtools::ScanBamParam(
            what = "qname",
            flag = Rsamtools::scanBamFlag(isSupplementaryAlignment = FALSE)
        )
        open(bam_file_con)
        repeat {

            ### Read the chunk of the BAM file
            chunk <- GenomicAlignments::readGAlignments(bam_file_con,
                                                        use.names = TRUE,
                                                        param = bam_param)
            if (length(chunk) == 0L)
                break

            ### Get chunk introns
            chunk_spliced <- chunk[GenomicAlignments::njunc(chunk) != 0,]
            chunk_introns <- GenomicAlignments::junctions(chunk_spliced) %>%
                unlist() %>%
                BiocGenerics::unstrand()

            ### Prepare read-level chunk intron summary
            chunk_intron_summary <- data.frame(
                read_id = names(chunk_introns),
                position = as.character(chunk_introns)
            ) %>%
                dplyr::group_by(.data$read_id) %>%
                dplyr::summarise(
                    intron_positions = paste0(.data$position, collapse = ",")
                )

            ### Prepare chunk read structure data
            chunk_read_summary <- data.frame(
                read_id = names(chunk_spliced),
                chromosome = as.character(GenomeInfoDb::seqnames(chunk_spliced)),
                start = BiocGenerics::start(chunk_spliced),
                end = BiocGenerics::end(chunk_spliced)
            )  %>%
                dplyr::left_join(chunk_intron_summary) %>%
                dplyr::select(-"read_id")
            bam_read_summary <- rbind(bam_read_summary,
                                      as.data.frame(chunk_read_summary))
        }
        close(bam_file_con)

        ## Collapse BAM file read structure data
        bam_read_summary <- bam_read_summary %>%
            dplyr::group_by(.data$intron_positions) %>%
            dplyr::summarise(
                read_count = dplyr::n(),
                chromosome = unique(.data$chromosome),
                start_positions = list(.data$start),
                end_positions = list(.data$end)
            ) %>%
            as.data.frame()

        return(bam_read_summary)
    }, BPPARAM = BPPARAM)
    read_summary <- do.call(rbind, read_summary)

    # Merge read structures from different BAM files
    read_summary <- read_summary %>%
        dplyr::group_by(.data$intron_positions) %>%
        dplyr::summarise(
            read_count = sum(.data$read_count),
            chromosome = unique(.data$chromosome),
            start_positions = merge_list(.data$start_positions),
            end_positions = merge_list(.data$end_positions)
        ) %>%
        as.data.frame()

    return(read_summary)
}
