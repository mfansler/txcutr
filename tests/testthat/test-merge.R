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

## Overlapping Genes
gr_multigene <- GRanges(
  seqnames=rep("chr1", 6),
  strand="+",
  ranges=IRanges(end=5000,
                 width=c(1000, 1000, 1000, 1200, 1200, 1200)),
  type=c("gene", "transcript", "exon",
         "gene", "transcript", "exon"),
  ID=c("gene_1", "tx_1", "exon_1",
       "gene_2", "tx_2", "exon_2"),
  Parent=c(NA, "gene_1", "tx_1",
           NA, "gene_2", "tx_2"),
  gene_id=c("gene_1", "gene_1", "gene_1",
            "gene_2", "gene_2", "gene_2"),
  tx_id=c(NA, "tx_1", "tx_1",
          NA, "tx_2", "tx_2"),
  exon_id=c(NA, NA, "exon_1",
            NA, NA, "exon_2"))

txdb_multigene <- makeTxDbFromGRanges(gr_multigene)

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
