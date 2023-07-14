test_that("EM works as expected", {

    # Testing if function throws the expected errors
    counts <- diag(c(100, 250))
    compatibility_matrix <- list(c(1, 1, 0),
                                 c(0, 0, 1))
    compatibility_matrix <- do.call(rbind, compatibility_matrix)
    tx_effective_lengths <- c(1, 1, 1)

    expect_error(EM(counts = 1,
                    compatibility_matrix = compatibility_matrix,
                    use_length_normalization = TRUE,
                    tx_effective_lengths = tx_effective_lengths,
                    maxiter = 250,
                    conv = 0.01),
                 regexp = "Not a matrix",
                 fixed = TRUE)
    expect_error(EM(counts = counts,
                    compatibility_matrix = 1,
                    use_length_normalization = TRUE,
                    tx_effective_lengths = tx_effective_lengths,
                    maxiter = 250,
                    conv = 0.01),
                 regexp = "Not a matrix",
                 fixed = TRUE)
    expect_error(EM(counts = counts,
                    compatibility_matrix = compatibility_matrix,
                    use_length_normalization = NULL,
                    tx_effective_lengths = tx_effective_lengths,
                    maxiter = 250,
                    conv = 0.01),
                 regexp = "Expecting a single value",
                 fixed = TRUE)
    expect_error(EM(counts = counts,
                    compatibility_matrix = compatibility_matrix,
                    use_length_normalization = TRUE,
                    tx_effective_lengths = 1,
                    maxiter = 250,
                    conv = 0.01),
                 regexp = "element-wise division",
                 fixed = TRUE)
    expect_error(EM(counts = counts,
                    compatibility_matrix = compatibility_matrix,
                    use_length_normalization = TRUE,
                    tx_effective_lengths = tx_effective_lengths,
                    maxiter = NULL,
                    conv = 0.01),
                 regexp = "Expecting a single value",
                 fixed = TRUE)
    expect_error(EM(counts = counts,
                    compatibility_matrix = compatibility_matrix,
                    use_length_normalization = TRUE,
                    tx_effective_lengths = tx_effective_lengths,
                    maxiter = 250,
                    conv = NULL),
                 regexp = "Expecting a single value",
                 fixed = TRUE)

    # Testing if function returns the expected output (case 1)
    counts <- diag(c(100, 250))
    compatibility_matrix <- list(c(1, 1, 0),
                                 c(0, 0, 1))
    compatibility_matrix <- do.call(rbind, compatibility_matrix)
    tx_effective_lengths <- c(1, 1, 1)

    expect_silent(
        tx_counts <- EM(counts = counts,
                        compatibility_matrix = compatibility_matrix,
                        use_length_normalization = FALSE,
                        tx_effective_lengths = tx_effective_lengths,
                        maxiter = 250,
                        conv = 0.01)
    )
    expect_true(is.matrix(tx_counts))
    expect_identical(dim(tx_counts), c(1L, 3L))
    expect_identical(as.numeric(tx_counts), c(50, 50, 250))

    expect_silent(
        tx_counts <- EM(counts = counts,
                        compatibility_matrix = compatibility_matrix,
                        use_length_normalization = TRUE,
                        tx_effective_lengths = tx_effective_lengths,
                        maxiter = 250,
                        conv = 0.01)
    )
    expect_true(is.matrix(tx_counts))
    expect_identical(dim(tx_counts), c(1L, 3L))
    expect_identical(as.numeric(tx_counts), c(50, 50, 250))

    # Testing if function returns the expected output (case 2)
    counts <- diag(c(100, 100, 100))
    compatibility_matrix <- list(c(1, 0),
                                 c(0, 1),
                                 c(1, 1))
    compatibility_matrix <- do.call(rbind, compatibility_matrix)
    tx_effective_lengths <- c(1.2, 1)

    expect_silent(
        tx_counts <- EM(counts = counts,
                        compatibility_matrix = compatibility_matrix,
                        use_length_normalization = FALSE,
                        tx_effective_lengths = tx_effective_lengths,
                        maxiter = 250,
                        conv = 0.01)
    )
    expect_true(is.matrix(tx_counts))
    expect_identical(dim(tx_counts), c(1L, 2L))
    expect_identical(as.numeric(tx_counts), c(150, 150))

    expect_silent(
        tx_counts <- EM(counts = counts,
                        compatibility_matrix = compatibility_matrix,
                        use_length_normalization = TRUE,
                        tx_effective_lengths = tx_effective_lengths,
                        maxiter = 250,
                        conv = 0.01)
    )
    expect_true(is.matrix(tx_counts))
    expect_identical(dim(tx_counts), c(1L, 2L))
    expect_identical(as.numeric(round(tx_counts)), c(144, 156))
})
