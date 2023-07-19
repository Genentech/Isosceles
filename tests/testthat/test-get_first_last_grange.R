test_that("get_first_last_grange works as expected", {

    # Preparing test data
    granges_1 <- GenomicRanges::GRanges(seqnames = c("chr1", "chr1"),
                                        ranges = IRanges::IRanges(
                                            start = c(100, 500),
                                            end = c(200, 750)),
                                        strand = c("+", "+"))
    granges_2 <- GenomicRanges::GRanges(seqnames = c("chr2", "chr2"),
                                        ranges = IRanges::IRanges(
                                            start = c(1000, 3000),
                                            end = c(1500, 6000)),
                                        strand = c("-", "-"))
    grangeslist <- GenomicRanges::GRangesList("feature_1" = granges_1,
                                              "feature_2" = granges_2)
    first_last_df <- data.frame(
        feature_id = c("feature_1", "feature_2"),
        first_grange = c("chr1:100-200:+", "chr2:3000-6000:-"),
        last_grange = c("chr1:500-750:+", "chr2:1000-1500:-")
    )

    # Testing if function throws the expected errors
    expect_error(get_first_last_grange(NULL),
                 regexp = "methods::is(object = feature_grangeslist, class2 =",
                 fixed = TRUE)

    # Testing if function returns the expected output
    expect_identical(get_first_last_grange(grangeslist), first_last_df)
})
