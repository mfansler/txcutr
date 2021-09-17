library(GenomicRanges)
library(GenomicFeatures)

gr_single_contig_w1000 <- GRanges(
  seqnames=rep("chr1", 3),
  strand="+",
  ranges=IRanges(end=5000,
                 width=1000),
  type=c("gene", "transcript", "exon"),
  ID=c("gene_1", "tx_1", "exon_1"),
  Parent=c(NA, "gene_1", "tx_1"),
  gene_id="gene_1",
  tx_id=c(NA, "tx_1", "tx_1"),
  exon_id=c(NA, NA, "exon_1"))

txdb_single_contig_w1000 <- makeTxDbFromGRanges(gr_single_contig_w1000)

expect_equal_applied <- function (x, y, .fn) {
  for (fn in .fn) {
    expect_equal(fn(x), fn(y))
  }
}

test_that("simple truncation works", {
  txdb_res_w500 <- truncateTxome(txdb_single_contig_w1000, maxTxLength=500)

  ## correct lengths
  expect_equal(width(genes(txdb_res_w500)), c(500))
  expect_equal(width(transcripts(txdb_res_w500)), c(500))
  expect_equal(width(exons(txdb_res_w500)), c(500))

  ## correct ends
  expect_equal(end(genes(txdb_res_w500)), end(genes(txdb_single_contig_w1000)))
  expect_equal(end(transcripts(txdb_res_w500)), end(transcripts(txdb_single_contig_w1000)))
  expect_equal(end(exons(txdb_res_w500)), end(exons(txdb_single_contig_w1000)))
})

test_that("idempotent", {
  txdb_res_w500_1 <- truncateTxome(txdb_single_contig_w1000, maxTxLength=500)
  txdb_res_w500_2 <- truncateTxome(txdb_res_w500_1, maxTxLength=500)

  expect_equal_applied(txdb_res_w500_1, txdb_res_w500_2, .fn=list(
    function (x) { width(genes(x)) },
    function (x) { width(transcripts(x)) },
    function (x) { width(exons(x)) },
    function (x) { start(genes(x)) },
    function (x) { start(transcripts(x)) },
    function (x) { start(exons(x)) },
    function (x) { end(genes(x)) },
    function (x) { end(transcripts(x)) },
    function (x) { end(exons(x)) }))
})
