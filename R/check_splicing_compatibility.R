#' Splicing compatibility check
#'
#' Check splicing compatibility between subject and target transcripts
#'
#' @param subject A list of vectors storing intron indices for subject
#' transcripts
#' @param target A list of vectors storing intron indices for target
#' transcripts
#' @return A data frame storing indices of compatible subject and target
#' transcripts
#' @keywords internal
check_splicing_compatibility <- function(subject, target) {

    # Check arguments
    assertthat::assert_that(is.list(subject))
    assertthat::assert_that(length(subject) > 0)
    assertthat::assert_that(is.list(target))
    assertthat::assert_that(length(target) > 0)

    # Find transcripts sharing at least one intron
    subject_df <- subject %>%
        tibble::enframe() %>%
        dplyr::rename(subject_idx = "name", intron_idx = "value") %>%
        tidyr::unchop("intron_idx")
    target_df <- target %>%
        tibble::enframe() %>%
        dplyr::rename(target_idx = "name", intron_idx = "value") %>%
        tidyr::unchop("intron_idx")
    compatibility_df <- subject_df %>%
        dplyr::inner_join(target_df, relationship = "many-to-many") %>%
        dplyr::select(-"intron_idx") %>%
        dplyr::distinct()

    # Check splicing compatibility
    compatibility_df$is_compatible <- mapply(
        is_ordered_subset,
        subject[compatibility_df$subject_idx],
        target[compatibility_df$target_idx]
    )

    # Prepare a data frame of compatible transcript indices
    compatibility_df <- compatibility_df %>%
        dplyr::filter(.data$is_compatible) %>%
        dplyr::select(-"is_compatible")

    return(compatibility_df)
}
