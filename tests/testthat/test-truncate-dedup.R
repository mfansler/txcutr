library(GenomicRanges)
library(GenomicFeatures)

############
## Mock Data
############

## Single Exon Gene
gr_contig <- GRanges(
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

txdb_contig <- makeTxDbFromGRanges(gr_contig)
txdb_contig_inv <- makeTxDbFromGRanges(invertStrand(gr_contig))

## Negative Strand
gr_contig_neg <- GRanges(
  seqnames=rep("chr1", 5),
  strand="-",
  ranges=IRanges(start=5000,
                 width=c(1000, 1000, 1000, 800, 800)),
  type=c("gene", "transcript", "exon", "transcript", "exon"),
  ID=c("gene_1", "tx_1", "exon_1", "tx_2", "exon_2"),
  Parent=c(NA, "gene_1", "tx_1", "gene_1", "tx_2"),
  gene_id="gene_1",
  tx_id=c(NA, "tx_1", "tx_1", "tx_2", "tx_2"),
  exon_id=c(NA, NA, "exon_1", NA, "exon_2"))

txdb_contig_neg <- makeTxDbFromGRanges(gr_contig_neg)
txdb_contig_neg_inv <- makeTxDbFromGRanges(invertStrand(gr_contig_neg))

########
## Tests
########

test_that("identical transcripts are merged, positive strand", {
  LENGTHS_TO_TEST <- c(100, 500)

  for (n in LENGTHS_TO_TEST) {
    txdb_res <- truncateTxome(txdb_contig, maxTxLength=n)

    n_txs <- length(transcripts(txdb_res))
    n_exons <- length(exons(txdb_res))
    expect_equal(n_txs, 1)
    expect_equal(n_exons, 1)
  }
})

test_that("identical transcripts are merged, negative strand", {
  LENGTHS_TO_TEST <- c(100, 500)

  for (n in LENGTHS_TO_TEST) {
    txdb_res <- truncateTxome(txdb_contig_neg, maxTxLength=n)

    n_txs <- length(transcripts(txdb_res))
    n_exons <- length(exons(txdb_res))
    expect_equal(n_txs, 1)
    expect_equal(n_exons, 1)
  }
})

test_that("non-identical transcripts are retained, positive strand", {
  LENGTHS_TO_TEST <- c(900)

  for (n in LENGTHS_TO_TEST) {
    txdb_res <- truncateTxome(txdb_contig, maxTxLength=n)

    n_txs <- length(transcripts(txdb_res))
    n_exons <- length(exons(txdb_res))
    expect_equal(n_txs, 2)
    expect_equal(n_exons, 2)
  }
})

test_that("non-identical transcripts are retained, negative strand", {
  LENGTHS_TO_TEST <- c(900)

  for (n in LENGTHS_TO_TEST) {
    txdb_res <- truncateTxome(txdb_contig_neg, maxTxLength=n)

    n_txs <- length(transcripts(txdb_res))
    n_exons <- length(exons(txdb_res))
    expect_equal(n_txs, 2)
    expect_equal(n_exons, 2)
  }
})

test_that("APA transcripts are retained, positive strand", {
  LENGTHS_TO_TEST <- c(100, 500)

  for (n in LENGTHS_TO_TEST) {
    txdb_res <- truncateTxome(txdb_contig_neg_inv, maxTxLength=n)

    n_txs <- length(transcripts(txdb_res))
    n_exons <- length(exons(txdb_res))
    expect_equal(n_txs, 2)
    expect_equal(n_exons, 2)
  }
})

test_that("APA transcripts are retained, negative strand", {
  LENGTHS_TO_TEST <- c(100, 500)

  for (n in LENGTHS_TO_TEST) {
    txdb_res <- truncateTxome(txdb_contig_inv, maxTxLength=n)

    n_txs <- length(transcripts(txdb_res))
    n_exons <- length(exons(txdb_res))
    expect_equal(n_txs, 2)
    expect_equal(n_exons, 2)
  }
})

