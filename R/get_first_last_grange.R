#' First and last feature extraction from a GRangesList
#'
#' Obtain the first and last features from a GRangesList object
#'
#' @param feature_grangeslist A GRangesList object
#' @return A data frame containing following columns:
#' \describe{
#'   \item{feature_id}{character vector of GRangesList names}
#'   \item{first_grange}{character vector of first feature positions}
#'   \item{last_grange}{character vector of last feature positions}
#' }
#' @keywords internal
get_first_last_grange <- function(feature_grangeslist) {

    # Check arguments
    assertthat::assert_that(class(feature_grangeslist) %in%
                                c("SimpleGRangesList", "CompressedGRangesList"))

    # Convert GRangesList to a feature data frame
    feature_granges <- BiocGenerics::sort(unlist(feature_grangeslist))
    S4Vectors::mcols(feature_granges) <- NULL
    S4Vectors::mcols(feature_granges)$feature_id <- names(feature_granges)
    S4Vectors::mcols(feature_granges)$idx <- seq_along(feature_granges)
    feature_df <- cbind(data.frame(grange = as.character(feature_granges),
                                   strand = BiocGenerics::strand(feature_granges)),
                        S4Vectors::mcols(feature_granges))

    # Get first and last features (+ strand)
    feature_df_plus <- dplyr::filter(feature_df, .data$strand == "+")
    feature_first_plus <- feature_df_plus %>%
        dplyr::group_by(.data$feature_id) %>%
        dplyr::slice_min(.data$idx) %>%
        dplyr::select(.data$feature_id, .data$grange)
    feature_last_plus <- feature_df_plus %>%
        dplyr::group_by(.data$feature_id) %>%
        dplyr::slice_max(.data$idx) %>%
        dplyr::select(.data$feature_id, .data$grange)

    # Get first and last features (- strand)
    feature_df_minus <- dplyr::filter(feature_df, .data$strand == "-")
    feature_first_minus <- feature_df_minus %>%
        dplyr::group_by(.data$feature_id) %>%
        dplyr::slice_max(.data$idx) %>%
        dplyr::select(.data$feature_id, .data$grange)
    feature_last_minus <- feature_df_minus %>%
        dplyr::group_by(.data$feature_id) %>%
        dplyr::slice_min(.data$idx) %>%
        dplyr::select(.data$feature_id, .data$grange)

    # Merge the first and last features for both strands
    feature_first_df <- rbind(feature_first_plus, feature_first_minus) %>%
        dplyr::rename(first_grange = .data$grange)
    feature_last_df <- rbind(feature_last_plus, feature_last_minus) %>%
        dplyr::rename(last_grange = .data$grange)
    feature_first_last_df <- data.frame(
        feature_id = names(feature_grangeslist)
    ) %>%
        dplyr::left_join(feature_first_df) %>%
        dplyr::left_join(feature_last_df) %>%
        dplyr::filter(!is.na(.data$first_grange), !is.na(.data$last_grange))
    return(feature_first_last_df)
}
