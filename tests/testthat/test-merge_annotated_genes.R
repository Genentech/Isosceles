test_that("merge_annotated_genes works as expected", {

    # Preparing test data
    gene_df <- data.frame(
        gene_id = c("GeneA", "GeneB", "GeneC"),
        gene_name = c("A", "B", "C")
    )
    intron_df <- data.frame(
        gene_id = c("GeneA", "GeneB", "GeneB", "GeneC", "GeneC"),
        position = c("chr1:100-200:+", "chr2:600-750:-", "chr2:800-900:-",
                     "chr2:600-750:-", "chr2:800-1050:-")
    )

    # Testing if function throws the expected errors
    expect_error(merge_annotated_genes(gene_df = NULL),
                 regexp = "gene_df is not a data frame",
                 fixed = TRUE)
    expect_error(merge_annotated_genes(gene_df = data.frame(),
                                       intron_df = NULL),
                 regexp = "intron_df is not a data frame",
                 fixed = TRUE)
    expect_error(merge_annotated_genes(gene_df = data.frame(),
                                       intron_df = data.frame()),
                 regexp = "gene_df does not have all of these name(s): 'gene_id'",
                 fixed = TRUE)
    expect_error(merge_annotated_genes(gene_df = data.frame(gene_id = "A"),
                                       intron_df = data.frame()),
                 regexp = "gene_df does not have all of these name(s): 'gene_name'",
                 fixed = TRUE)
    expect_error(merge_annotated_genes(gene_df = data.frame(gene_id = "A",
                                                            gene_name = "A"),
                                       intron_df = data.frame()),
                 regexp = "intron_df does not have all of these name(s): 'gene_id'",
                 fixed = TRUE)
    expect_error(merge_annotated_genes(gene_df = data.frame(gene_id = "A",
                                                            gene_name = "A"),
                                       intron_df = data.frame(gene_id = "A")),
                 regexp = "intron_df does not have all of these name(s): 'position'",
                 fixed = TRUE)

    # Testing if function returns the expected output
    expect_silent(merged_df <- merge_annotated_genes(gene_df = gene_df,
                                                     intron_df = intron_df))
    expect_true(is.data.frame(merged_df))
    expect_identical(dim(merged_df), c(3L, 4L))
    expect_identical(colnames(merged_df),
                     c("gene_id", "gene_name", "merged_gene_id",
                       "merged_gene_name"))
    expect_identical(merged_df$gene_id, gene_df$gene_id)
    expect_identical(merged_df$gene_name, gene_df$gene_name)
    expect_identical(merged_df$merged_gene_id,
                     c("GeneA", "GeneB,GeneC", "GeneB,GeneC"))
    expect_identical(merged_df$merged_gene_name, c("A", "B,C", "B,C"))
})
