library(GenomicRanges)
library(txdbmaker)

############
## Mock Data
############

default_meta <- data.frame(name="Truncation End", value="3prime")

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

txdb_contig <- .suppressTxDbGenomeWarning(makeTxDbFromGRanges(gr_contig, metadata = default_meta))

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

txdb_contig_neg <- .suppressTxDbGenomeWarning(makeTxDbFromGRanges(gr_contig_neg, metadata = default_meta))

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

txdb_transitive <- .suppressTxDbGenomeWarning(makeTxDbFromGRanges(gr_transitive, metadata = default_meta))
txdb_transitive_neg <- .suppressTxDbGenomeWarning(makeTxDbFromGRanges(invertStrand(gr_transitive), metadata = default_meta))

## Truncation + Merge on complex transcript structure
gr_complex <- GRanges(
  seqnames="chr1",
  strand="+",
  ranges=IRanges(start=c(1000, 
                         1010, 1010, 4300, 7650, 
                         1000, 1000, 4300, 7700, 
                         1100, 1100, 7300, 
                         1000,
                         1100, 1100, 7300),
                 end=c(8000, 
                       8000, 1300, 4700, 8000, 
                       7990, 1300, 4700, 7990, 
                       7900, 1700, 7900, 
                       8000,
                       7900, 1700, 7900)),
  type=c("gene",
         "transcript", "exon", "exon", "exon",
         "transcript", "exon", "exon", "exon",
         "transcript", "exon","exon",
         "gene",
         "transcript", "exon","exon"),
  ID=c("gene_1",
       "tx_1", "exon_1-1", "exon_1-2", "exon_1-3",
       "tx_2", "exon_2-1", "exon_2-2", "exon_2-3",
       "tx_3", "exon_3-1", "exon_3-2",
       "gene_2",
       "tx_4", "exon_4-1", "exon_4-2"),
  Parent=c(NA,
           "gene_1", "tx_1", "tx_1", "tx_1",
           "gene_1", "tx_2", "tx_2", "tx_2",
           "gene_1", "tx_3", "tx_3",
           NA,
           "gene_2", "tx_4", "tx_4"),
  gene_id=c(rep("gene_1", 12), rep("gene_2", 4)),
  tx_id=c(NA, "tx_1", "tx_1", "tx_1", "tx_1", "tx_2", "tx_2", "tx_2", "tx_2", "tx_3", "tx_3", "tx_3", NA, "tx_4", "tx_4", "tx_4"),
  exon_id=c(NA, NA, "exon_1-1", "exon_1-2", "exon_1-3", NA, "exon_2-1", "exon_2-2", "exon_2-3", NA, "exon_3-1", "exon_3-2", NA, NA, "exon_4-1", "exon_4-2"))

txdb_complex <- .suppressTxDbGenomeWarning(makeTxDbFromGRanges(gr_complex))
txdb_complex_neg <- .suppressTxDbGenomeWarning(makeTxDbFromGRanges(invertStrand(gr_complex)))

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

test_that("merging is transitive, negative strand", {
  df <- generateMergeTable(txdb_transitive_neg, minDistance=150)
  n_txdb_txs <- length(transcripts(txdb_transitive_neg))
  n_txs_in <- length(unique(df$tx_in))
  n_txs_out <- length(unique(df$tx_out))
  tx_distal <- get_distal_tx_name(transcripts(txdb_transitive_neg))
  expect_equal(n_txs_in, n_txdb_txs)
  expect_equal(n_txs_out, 1)
  expect_equal(df$tx_out, rep(tx_distal, nrow(df)))
})

test_that("mergin works after truncation, 3prime positive strand", {
  txdb_trunc <- truncate3primeTxome(txdb_complex, maxTxLength = 400, quiet = T)
  df <- generateMergeTable(txdb_trunc, minDistance=150)
  n_txdb_txs <- length(transcripts(txdb_trunc))
  n_txs_in <- length(unique(df$tx_in))
  n_txs_out <- length(unique(df$tx_out))
  expect_equal(n_txs_in, n_txdb_txs)
  expect_equal(n_txs_out, 2)
  expect_equal(df$tx_in, c("tx_1", "tx_2", "tx_3", "tx_4"))
  expect_equal(df$tx_out, c("tx_1", "tx_1", "tx_1", "tx_4"))
})

test_that("mergin works after truncation, 3prime negative strand", {
  txdb_trunc <- truncate3primeTxome(txdb_complex_neg, maxTxLength = 400, quiet = T)
  df <- generateMergeTable(txdb_trunc, minDistance=150)
  n_txdb_txs <- length(transcripts(txdb_trunc))
  n_txs_in <- length(unique(df$tx_in))
  n_txs_out <- length(unique(df$tx_out))
  expect_equal(n_txs_in, n_txdb_txs)
  expect_equal(n_txs_out, 2)
  expect_equal(df$tx_in, c("tx_1", "tx_2", "tx_3", "tx_4"))
  expect_equal(df$tx_out, c("tx_2", "tx_2", "tx_2", "tx_4"))
})

test_that("mergin works after truncation, 5prime positive strand", {
  txdb_trunc <- truncate5primeTxome(txdb_complex, maxTxLength = 400, quiet = T)
  df <- generateMergeTable(txdb_trunc, minDistance=150)
  n_txdb_txs <- length(transcripts(txdb_trunc))
  n_txs_in <- length(unique(df$tx_in))
  n_txs_out <- length(unique(df$tx_out))
  expect_equal(n_txs_in, n_txdb_txs)
  expect_equal(n_txs_out, 2)
  expect_equal(df$tx_in, c("tx_1", "tx_2", "tx_3", "tx_4"))
  expect_equal(df$tx_out, c("tx_2", "tx_2", "tx_2", "tx_4"))
})

test_that("mergin works after truncation, 5prime negative strand", {
  txdb_trunc <- truncate5primeTxome(txdb_complex_neg, maxTxLength = 400, quiet = T)
  df <- generateMergeTable(txdb_trunc, minDistance=150)
  n_txdb_txs <- length(transcripts(txdb_trunc))
  n_txs_in <- length(unique(df$tx_in))
  n_txs_out <- length(unique(df$tx_out))
  expect_equal(n_txs_in, n_txdb_txs)
  expect_equal(n_txs_out, 2)
  expect_equal(df$tx_in, c("tx_1", "tx_2", "tx_3", "tx_4"))
  expect_equal(df$tx_out, c("tx_1", "tx_1", "tx_1", "tx_4"))
})

test_that("mergin works after truncation, 3prime positive strand (smaller merge length)", {
  txdb_trunc <- truncate3primeTxome(txdb_complex, maxTxLength = 400, quiet = T)
  df <- generateMergeTable(txdb_trunc, minDistance=90)
  n_txdb_txs <- length(transcripts(txdb_trunc))
  n_txs_in <- length(unique(df$tx_in))
  n_txs_out <- length(unique(df$tx_out))
  expect_equal(n_txs_in, n_txdb_txs)
  expect_equal(n_txs_out, 3)
  expect_equal(df$tx_in, c("tx_1", "tx_2", "tx_3", "tx_4"))
  expect_equal(df$tx_out, c("tx_1", "tx_1", "tx_3", "tx_4"))
})
