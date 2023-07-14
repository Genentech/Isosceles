test_that("prepare_reference_spliced_transcripts works as expected", {

    # Preparing test data
    gtf_file <- system.file(
        "extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf",
        package = "Isosceles"
    )
    anno_data <- prepare_reference_annotations(gtf_file, is_technical = FALSE)

    # Testing if function throws the expected errors
    expect_error(prepare_reference_spliced_transcripts(anno_data = NULL),
                 regexp = "anno_data is not a list",
                 fixed = TRUE)
    expect_error(prepare_reference_spliced_transcripts(anno_data = list()),
                 regexp = "anno_data does not have all of these name(s): 'transcript_df'",
                 fixed = TRUE)
    expect_error(prepare_reference_spliced_transcripts(anno_data = list(transcript_df = 42)),
                 regexp = "anno_data$transcript_df is not a data frame",
                 fixed = TRUE)
    expect_error(prepare_reference_spliced_transcripts(anno_data = anno_data,
                                                       bin_size = NULL),
                 regexp = "bin_size is not a count",
                 fixed = TRUE)

    # Testing if function returns the expected output
    expect_silent(
        tx_list <- prepare_reference_spliced_transcripts(
            anno_data = anno_data
        )
    )
    expect_true(is.list(tx_list))
    expect_identical(length(tx_list), 4L)
    expect_identical(names(tx_list),
                     c("tx_df", "tx_granges","tx_exon_granges_list",
                       "tx_intron_granges_list"))
    expect_identical(dim(tx_list$tx_df), c(97L, 12L))
    expect_identical(colnames(tx_list$tx_df),
                     c("hash_id", "position", "intron_positions",
                       "gene_id", "gene_name", "compatible_gene_ids",
                       "compatible_gene_names", "compatible_tx",
                       "splicing_support_level", "fivethree_support_level",
                       "read_count", "relative_expression"))
    expect_identical(sum(is.na(tx_list$tx_df$hash_id)), 0L)
    expect_identical(sum(is.na(tx_list$tx_df$position)), 0L)
    expect_identical(sum(is.na(tx_list$tx_df$intron_positions)), 0L)
    expect_identical(sum(is.na(tx_list$tx_df$gene_id)), 0L)
    expect_identical(sum(is.na(tx_list$tx_df$gene_name)), 0L)
    expect_identical(sum(is.na(tx_list$tx_df$compatible_gene_ids)), 0L)
    expect_identical(sum(is.na(tx_list$tx_df$compatible_gene_names)), 0L)
    expect_identical(sum(is.na(tx_list$tx_df$compatible_tx)), 0L)
    expect_identical(names(table(tx_list$tx_df$splicing_support_level)),
                     "AP")
    expect_identical(as.numeric(table(tx_list$tx_df$splicing_support_level)),
                     97)
    expect_identical(names(table(tx_list$tx_df$fivethree_support_level)),
                     "FL")
    expect_identical(as.numeric(table(tx_list$tx_df$fivethree_support_level)),
                     97)
    expect_true(all(is.na(tx_list$tx_df$read_count)))
    expect_true(all(is.na(tx_list$tx_df$relative_expression)))
    expect_true(class(tx_list$tx_granges) == "GRanges")
    expect_identical(length(tx_list$tx_granges), 97L)
    expect_identical(
        names(table(BiocGenerics::strand(tx_list$tx_granges))),
        c("+", "-")
    )
    expect_identical(
        as.numeric(table(BiocGenerics::strand(tx_list$tx_granges))),
        c(43, 54)
    )
    expect_true(grepl("GRangesList", class(tx_list$tx_exon_granges_list)))
    expect_identical(length(tx_list$tx_exon_granges_list), 97L)
    expect_identical(length(unlist(tx_list$tx_exon_granges_list)), 823L)
    expect_true(grepl("GRangesList", class(tx_list$tx_intron_granges_list)))
    expect_identical(length(tx_list$tx_intron_granges_list), 97L)
    expect_identical(length(unlist(tx_list$tx_intron_granges_list)), 726L)
})
