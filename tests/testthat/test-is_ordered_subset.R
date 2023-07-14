test_that("is_ordered_subset works as expected", {

    # Testing if function returns the expected output
    expect_true(is_ordered_subset(1L, 1:3))
    expect_true(is_ordered_subset(1:3, 1:3))
    expect_true(is_ordered_subset(1:3, 1:10))
    expect_true(is_ordered_subset(7:9, 1:10))
    expect_true(is_ordered_subset(5:10, 1:10))
    expect_false(is_ordered_subset(c(7L, 9L), 1:10))
    expect_false(is_ordered_subset(1:20, 1:10))
    expect_false(is_ordered_subset(13L, 1:10))
})
