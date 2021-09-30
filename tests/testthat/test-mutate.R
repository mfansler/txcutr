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
grl_multi_exon <- exonsBy(txdb_multi_exon, use.names=TRUE)

gr_multi_tx <- GRanges(
  seqnames=rep("chr1", 5),
  strand="+",
  ranges=IRanges(end=5000,
                 width=c(1000, 1000, 1000, 800, 800)),
  type=c("gene", "transcript", "exon", "transcript", "exon"),
  ID=c("gene_1", "tx_1", "exon_1", "tx_2", "exon_2"),
  Parent=c(NA, "gene_1", "tx_1", "gene_1", "tx_2"),
  gene_id="gene_1",
  tx_id=c(NA, "tx_1", "tx_1", "tx_2", "tx_2"),
  exon_id=c(NA, NA, "exon_1", NA, "exon_2"))

txdb_multi_tx <- makeTxDbFromGRanges(gr_multi_tx)
grl_multi_tx <- exonsBy(txdb_multi_tx, use.names=TRUE)

########
## Tests
########

test_that("adds correct transcript names, multi-exon", {
  grl <- .mutateEach(grl_multi_exon, new_col=names(grl_multi_exon))
  for (tx in names(grl)) {
    expect_setequal(grl[[tx]]$new_col, tx)
  }
})

test_that("adds correct transcript names, multi-tx", {
  grl <- .mutateEach(grl_multi_tx, new_col=names(grl_multi_tx))
  for (tx in names(grl)) {
    expect_setequal(grl[[tx]]$new_col, tx)
  }
})

test_that("SimpleGRangesList and CompressedGRangesList match", {
  sgrl_in <- as(grl_multi_tx, "SimpleGRangesList")
  cgrl_in <- as(grl_multi_tx, "CompressedGRangesList")
  sgrl_out <- .mutateEach(sgrl_in, new_col=names(sgrl_in))
  cgrl_out <- .mutateEach(cgrl_in, new_col=names(cgrl_in))
  expect_equal(sgrl_out, cgrl_out)
})

test_that("incorrect assignment lengths are caught", {
  bad_col_long <- c(0, seq_along(names(grl_multi_tx)))
  bad_col_short <- seq_along(names(grl_multi_tx))[-1]
  empty_col <- numeric(0)
  expect_error(.mutateEach(grl_multi_tx, new_col=bad_col_long))
  expect_error(.mutateEach(grl_multi_tx, new_col=bad_col_short))
  expect_error(.mutateEach(grl_multi_tx, new_col=empty_col))
})
