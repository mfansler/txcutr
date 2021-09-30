library(GenomicRanges)
library(GenomicFeatures)

############
## Mock Data
############

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

########
## Tests
########

test_that("outputs uncompressed GTF on .gtf", {
  filename <- file.path(tempdir(), "transitive.gtf")
  on.exit(unlink(filename))
  expect_invisible(exportGTF(txdb_transitive, file=filename))
  con <- file(filename)
  on.exit(close(con), add=TRUE, after=FALSE)
  expect_identical(summary(con)$class, "file")
})

test_that("outputs compressed GTF on .gz", {
  filename <- file.path(tempdir(), "transitive.gtf.gz")
  on.exit(unlink(filename))
  expect_invisible(exportGTF(txdb_transitive, file=filename))
  con <- file(filename)
  on.exit(close(con), add=TRUE, after=FALSE)
  expect_identical(summary(con)$class, "gzfile")
})

test_that("outputs uncompressed TSV on .tsv", {
  filename <- file.path(tempdir(), "transitive.merge.tsv")
  on.exit(unlink(filename))
  expect_invisible(exportMergeTable(txdb_transitive, file=filename))
  con <- file(filename)
  on.exit(close(con), add=TRUE, after=FALSE)
  expect_identical(summary(con)$class, "file")
})

test_that("outputs compressed TSV on .tsv.gz", {
  filename <- file.path(tempdir(), "transitive.merge.tsv.gz")
  on.exit(unlink(filename))
  expect_invisible(exportMergeTable(txdb_transitive, file=filename))
  con <- file(filename)
  on.exit(close(con), add=TRUE, after=FALSE)
  expect_identical(summary(con)$class, "gzfile")
})
