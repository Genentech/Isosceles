test_that("check_splicing_compatibility works as expected", {

    # Preparing test data
    subject <- list(
        1L, 1:2, 1:3, 1:10, 15:20
    )
    target <- list(
        1:2, 1:5, 1:10, 11:22
    )

    # Testing if function throws the expected errors
    expect_error(check_splicing_compatibility(subject = NULL),
                 regexp = "subject is not a list",
                 fixed = TRUE)
    expect_error(check_splicing_compatibility(subject = list()),
                 regexp = "length(subject) not greater than 0",
                 fixed = TRUE)
    expect_error(check_splicing_compatibility(subject = subject,
                                              target = NULL),
                 regexp = "target is not a list",
                 fixed = TRUE)
    expect_error(check_splicing_compatibility(subject = subject,
                                              target = list()),
                 regexp = "length(target) not greater than 0",
                 fixed = TRUE)

    # Testing if function returns the expected output
    expect_message(
        compatibility_df <- check_splicing_compatibility(
            subject = subject, target = target
        ),
        regexp = "intron_idx", fixed = TRUE
    )
    expect_true(is.data.frame(compatibility_df))
    expect_identical(colnames(compatibility_df),
                     c("subject_idx", "target_idx"))
    expect_identical(dim(compatibility_df), c(10L, 2L))
    expect_identical(compatibility_df$subject_idx,
                     c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 4L, 5L))
    expect_identical(compatibility_df$target_idx,
                     c(1L, 2L, 3L, 1L, 2L, 3L, 2L, 3L, 3L, 4L))
})
