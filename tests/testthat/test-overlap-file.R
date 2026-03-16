library(GenomicRanges)
library(txdbmaker)

############
## Mock Data
############

## Four transcripts with different boundaries in multiple genes (two transcripts
## overlaps after truncation)
gr_collapse_test <- GRanges(
  seqnames = rep("chr1", 11),
  strand = "+",
  ranges = IRanges(
    start = c(
      1000, # gene_1
      1000, # gene_2
      2000, 2000, # tx_1
      1000, 1000, 1500, # tx_2,
      1500, 2500, # tx_3
      1500, 2500
    ), # tx_4
    end = c(
      6000, # gene_1
      5000, # gene_2
      5000, 5000, # tx_1
      5000, 1400, 5000, # tx_2
      6000, 6000, # tx_3
      6000, 6000
    ) # tx_4
  ),
  type = c(
    "gene", "gene",
    "transcript", "exon",
    "transcript", "exon", "exon",
    "transcript", "exon",
    "transcript", "exon"
  ),
  ID = c(
    "gene_1", "gene_2",
    "tx_1", "exon_1",
    "tx_2", "exon_2", "exon_3",
    "tx_3", "exon_4",
    "tx_4", "exon_5"
  ),
  Parent = c(
    NA, NA,
    "gene_1", "tx_1",
    "gene_1", "tx_2", "tx_2",
    "gene_1", "tx_3",
    "gene_2", "tx_4"
  ),
  gene_id = c(
    "gene_1", "gene_2",
    "gene_1", "gene_1", "gene_1", "gene_1", "gene_1", "gene_1", "gene_1",
    "gene_2", "gene_2"
  ),
  tx_id = c(
    NA, NA,
    "tx_1", "tx_1",
    "tx_2", "tx_2", "tx_2",
    "tx_3", "tx_3",
    "tx_4", "tx_4"
  ),
  exon_id = c(
    NA, NA,
    NA, "exon_1",
    NA, "exon_2", "exon_3",
    NA, "exon_4",
    NA, "exon_5"
  )
)

txdb_collapse_test <- .suppressTxDbGenomeWarning(makeTxDbFromGRanges(gr_collapse_test))

## Three transcripts with identical boundaries after truncation but different
## internal exon structure
gr_exon_structure_test <- GRanges(
  seqnames = rep("chr1", 10),
  strand = "+",
  ranges = IRanges(
    start = c(
      1000, # gene_1
      2000, 2000, # tx_1 (single exon)
      2000, 2000, 4000, # tx_2 (two exons)
      2000, 2000, 3000, 4000
    ), # tx_3 (three exons)
    end = c(
      5000, # gene_1
      5000, 5000, # tx_1
      5000, 2500, 5000, # tx_2
      5000, 2500, 3500, 5000
    ) # tx_3
  ),
  type = c(
    "gene",
    "transcript", "exon",
    "transcript", "exon", "exon",
    "transcript", "exon", "exon", "exon"
  ),
  ID = c(
    "gene_1",
    "tx_1", "exon_1",
    "tx_2", "exon_2", "exon_3",
    "tx_3", "exon_4", "exon_5", "exon_6"
  ),
  Parent = c(
    NA,
    "gene_1", "tx_1",
    "gene_1", "tx_2", "tx_2",
    "gene_1", "tx_3", "tx_3", "tx_3"
  ),
  gene_id = "gene_1",
  tx_id = c(
    NA,
    "tx_1", "tx_1",
    "tx_2", "tx_2", "tx_2",
    "tx_3", "tx_3", "tx_3", "tx_3"
  ),
  exon_id = c(
    NA,
    NA, "exon_1",
    NA, "exon_2", "exon_3",
    NA, "exon_4", "exon_5", "exon_6"
  )
)

txdb_exon_structure_test <- .suppressTxDbGenomeWarning(makeTxDbFromGRanges(gr_exon_structure_test))

########
## Tests
########

test_that("identical transcripts are collapsed after truncation, positive strand 3'", {
  LENGTHS_TO_TEST <- c(100, 500)

  for (n in LENGTHS_TO_TEST) {
    temp_file <- withr::local_tempfile(fileext = ".csv")
    txdb_res <- truncate3primeTxome(txdb_collapse_test, maxTxLength = n, overlapFile = temp_file, quiet = T)

    ## correct overlaps removal
    test_overlap <- data.frame(query_transcript = "tx_2", subject_transcript = "tx_1")
    read_overlap <- read.table(temp_file, sep = ",", header = T)[, c("query_transcript", "subject_transcript")]

    expect_equal(test_overlap, read_overlap)
    expect_equal(nrow(read_overlap), 1)
  }
})

test_that("identical transcripts are collapsed after truncation, positive strand 5'", {
  LENGTHS_TO_TEST <- c(100, 500)

  for (n in LENGTHS_TO_TEST) {
    temp_file <- withr::local_tempfile(fileext = ".csv")
    txdb_res <- truncate5primeTxome(txdb_collapse_test, maxTxLength = n, overlapFile = temp_file, quiet = T)

    ## correct overlaps removal
    test_overlap <- data.frame(query_transcript = logical(0), subject_transcript = logical(0))
    read_overlap <- read.table(temp_file, sep = ",", header = T)[, c("query_transcript", "subject_transcript")]

    expect_equal(test_overlap, read_overlap)
    expect_equal(nrow(read_overlap), 0)
  }
})

test_that("transcripts with same coordinates but different exon structure", {
  LENGTHS_TO_TEST <- c(100, 500)

  for (n in LENGTHS_TO_TEST) {
    temp_file <- withr::local_tempfile(fileext = ".csv")
    txdb_res <- truncate3primeTxome(txdb_exon_structure_test, maxTxLength = n, overlapFile = temp_file, quiet = TRUE)

    # Before truncation: 3 transcripts with different exon structures
    expect_equal(length(transcripts(txdb_exon_structure_test)), 3)

    # Read overlap file
    read_overlap <- read.table(temp_file, sep = ",", header = TRUE)[, c("query_transcript", "subject_transcript")]

    # Should detect 3 overlaps (tx_2 vs tx_1, tx_3 vs tx_1 and tx_3 vs tx_2)
    expect_equal(nrow(read_overlap), 3)

    # All three transcript IDs should appear in the overlap file
    all_tx_in_overlaps <- unique(c(read_overlap$query_transcript, read_overlap$subject_transcript))
    expect_setequal(all_tx_in_overlaps, c("tx_1", "tx_2", "tx_3"))

    # After collapsing, should have only 1 transcript remaining
    expect_equal(length(transcripts(txdb_res)), 1)
  }
})
