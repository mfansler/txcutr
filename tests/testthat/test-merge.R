library(GenomicRanges)
library(GenomicFeatures)

############
## Mock Data
############

## Single Exon Gene
gr_contig <- GRanges(
  seqnames=rep("chr1", 5),
  strand="+",
  ranges=IRanges(start=5000,
                 width=c(1000, 1000, 1000, 900, 900)),
  type=c("gene", "transcript", "exon", "transcript", "exon"),
  ID=c("gene_1", "tx_1", "exon_1", "tx_2", "exon_2"),
  Parent=c(NA, "gene_1", "tx_1", "gene_1", "tx_2"),
  gene_id="gene_1",
  tx_id=c(NA, "tx_1", "tx_1", "tx_2", "tx_2"),
  exon_id=c(NA, NA, "exon_1", NA, "exon_2"))

txdb_contig <- makeTxDbFromGRanges(gr_contig)

## Negative Strand
gr_contig_neg <- GRanges(
  seqnames=rep("chr1", 5),
  strand="-",
  ranges=IRanges(end=5000,
                 width=c(1000, 1000, 1000, 900, 900)),
  type=c("gene", "transcript", "exon", "transcript", "exon"),
  ID=c("gene_1", "tx_1", "exon_1", "tx_2", "exon_2"),
  Parent=c(NA, "gene_1", "tx_1", "gene_1", "tx_2"),
  gene_id="gene_1",
  tx_id=c(NA, "tx_1", "tx_1", "tx_2", "tx_2"),
  exon_id=c(NA, NA, "exon_1", NA, "exon_2"))

txdb_contig_neg <- makeTxDbFromGRanges(gr_contig_neg)

## Transitive Positive
gr_transitive <- GRanges(
  seqnames="chr1",
  strand="+",
  ranges=IRanges(start=c(5000, 5000, 5000, 5100, 5100, 5200, 5200),
                 end=c(5400, 5200, 5200, 5300, 5300, 5400, 5400)),
  type=c("gene",
         "transcript", "exon",
         "transcript", "exon",
         "transcript", "exon"),
  ID=c("gene_1",
       "tx_1", "exon_1",
       "tx_2", "exon_2",
       "tx_3", "exon_3"),
  Parent=c(NA,
           "gene_1", "tx_1",
           "gene_1", "tx_2",
           "gene_1", "tx_3"),
  gene_id="gene_1",
  tx_id=c(NA, "tx_1", "tx_1", "tx_2", "tx_2", "tx_3", "tx_3"),
  exon_id=c(NA, NA, "exon_1", NA, "exon_2", NA, "exon_3"))

txdb_transitive <- makeTxDbFromGRanges(gr_transitive)
txdb_transitive_neg <- makeTxDbFromGRanges(invertStrand(gr_transitive))

########
## Tests
########

test_that("nearby transcripts are merged, positive strand", {
  df <- generateMergeTable(txdb_contig, minDistance=200)
  n_txdb_txs <- length(transcripts(txdb_contig))
  n_txs_in <- length(unique(df$tx_in))
  n_txs_out <- length(unique(df$tx_out))
  expect_equal(n_txs_in, n_txdb_txs)
  expect_equal(n_txs_out, 1)
})

test_that("nearby transcripts are merged, negative strand", {
  df <- generateMergeTable(txdb_contig_neg, minDistance=200)
  n_txdb_txs <- length(transcripts(txdb_contig))
  n_txs_in <- length(unique(df$tx_in))
  n_txs_out <- length(unique(df$tx_out))
  expect_equal(n_txs_in, n_txdb_txs)
  expect_equal(n_txs_out, 1)
})

test_that("far transcripts are unmerged, positive strand", {
  df <- generateMergeTable(txdb_contig_neg, minDistance=50)
  n_txdb_txs <- length(transcripts(txdb_contig))
  n_txs_in <- length(unique(df$tx_in))
  n_txs_out <- length(unique(df$tx_out))
  expect_equal(n_txs_in, n_txdb_txs)
  expect_equal(n_txs_out, 2)
})

test_that("far transcripts are unmerged, negative strand", {
  df <- generateMergeTable(txdb_contig_neg, minDistance=50)
  n_txdb_txs <- length(transcripts(txdb_contig))
  n_txs_in <- length(unique(df$tx_in))
  n_txs_out <- length(unique(df$tx_out))
  expect_equal(n_txs_in, n_txdb_txs)
  expect_equal(n_txs_out, 2)
})

test_that("merging is transitive, positive strand", {
  df <- generateMergeTable(txdb_transitive, minDistance=150)
  n_txdb_txs <- length(transcripts(txdb_transitive))
  n_txs_in <- length(unique(df$tx_in))
  n_txs_out <- length(unique(df$tx_out))
  tx_distal <- get_distal_tx_name(transcripts(txdb_transitive))
  expect_equal(n_txs_in, n_txdb_txs)
  expect_equal(n_txs_out, 1)
  expect_equal(df$tx_out, rep(tx_distal, nrow(df)))
})

test_that("merging is transitive, negaative strand", {
  df <- generateMergeTable(txdb_transitive_neg, minDistance=150)
  n_txdb_txs <- length(transcripts(txdb_transitive_neg))
  n_txs_in <- length(unique(df$tx_in))
  n_txs_out <- length(unique(df$tx_out))
  tx_distal <- get_distal_tx_name(transcripts(txdb_transitive_neg))
  expect_equal(n_txs_in, n_txdb_txs)
  expect_equal(n_txs_out, 1)
  expect_equal(df$tx_out, rep(tx_distal, nrow(df)))
})
