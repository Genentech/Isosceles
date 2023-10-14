#' Vector subset check
#'
#' Checks if a vector is an ordered subset of another vector.
#'
#' @param x A vector.
#' @param y A vector.
#' @return A logical scalar indicating if x is an ordered subset of y.
#' @keywords internal
is_ordered_subset <- function(x, y) {
    if (identical(x, y)) return(TRUE)
    if (length(x) >= length(y)) return(FALSE)
    if (!all(x %in% y)) return(FALSE)
    y_start_idx <- match(x[1], y)
    y_subset <- y_start_idx:(y_start_idx + length(x) - 1)
    return(identical(x, y[y_subset]))
}
