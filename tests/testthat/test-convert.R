library(GenomicRanges)
library(GenomicFeatures)

############
## Mock Data
############

## Multi Exon
gr_multi_exon <- GRanges(
  seqnames=rep("chr1", 5),
  strand="+",
  ranges=IRanges(start=c(1000, 1000, 1000, 5000, 9000),
                 end=c(9299, 9299, 1299, 5299, 9299)),
  type=c("gene", "transcript", "exon", "exon", "exon"),
  ID=c("gene_1", "tx_1", "exon_1", "exon_2", "exon_3"),
  Parent=c(NA, "gene_1", "tx_1", "tx_1", "tx_1"),
  gene_id="gene_1",
  tx_id=c(NA, "tx_1", "tx_1", "tx_1", "tx_1"),
  exon_id=c(NA, NA, "exon_1", "exon_2", "exon_3"))

txdb_multi_exon <- makeTxDbFromGRanges(gr_multi_exon)

## Negative Strand
txdb_multi_exon_neg <- makeTxDbFromGRanges(invertStrand(gr_multi_exon))

########
## Tests
########

test_that("conversion yields expected fields, positive strand", {
  grl <- txdbToGRangesList(txdb_multi_exon)
  expect_s4_class(grl, "GRangesList")
  expect_named(grl, c("genes", "transcripts", "exons"))
  expect_setequal(grl[['exons']]$type, c("exon"))
  expect_setequal(grl[['transcripts']]$type, c("transcript"))
  expect_setequal(grl[['genes']]$type, c("gene"))
})

test_that("conversion yields expected fields, negative strand", {
  grl <- txdbToGRangesList(txdb_multi_exon_neg)
  expect_s4_class(grl, "GRangesList")
  expect_named(grl, c("genes", "transcripts", "exons"))
  expect_setequal(grl[['exons']]$type, c("exon"))
  expect_setequal(grl[['transcripts']]$type, c("transcript"))
  expect_setequal(grl[['genes']]$type, c("gene"))
})

test_that("conversion maintains hierarchy, positive strand", {
  grl <- txdbToGRangesList(txdb_multi_exon)
  expect_setequal(grl[['transcripts']]$gene_id, grl[['genes']]$gene_id)
  expect_setequal(grl[['exons']]$gene_id, grl[['genes']]$gene_id)
  expect_setequal(grl[['exons']]$tx_name, grl[['transcripts']]$tx_name)
})

test_that("conversion maintains hierarchy, negative strand", {
  grl <- txdbToGRangesList(txdb_multi_exon_neg)
  expect_setequal(grl[['transcripts']]$gene_id, grl[['genes']]$gene_id)
  expect_setequal(grl[['exons']]$gene_id, grl[['genes']]$gene_id)
  expect_setequal(grl[['exons']]$tx_name, grl[['transcripts']]$tx_name)
})

