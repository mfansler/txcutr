#' @rdname truncateTxome
#' @export
setGeneric("truncateTxome", function(txdb, maxTxLength = 500, txEnd = "3prime", overlapFile = NULL, BPPARAM = bpparam(), ...) {
  standardGeneric("truncateTxome")
}
)

#' @rdname truncateTxome
#' @export
setGeneric("truncate3primeTxome", function(txdb, maxTxLength = 500, overlapFile = NULL, BPPARAM = bpparam(), quiet = FALSE, ...) {
  standardGeneric("truncate3primeTxome")
} 
)

#' @rdname truncateTxome
#' @export
setMethod("truncate3primeTxome", "TxDb", function(txdb, maxTxLength = 500, overlapFile = NULL, 
                                                  BPPARAM = bpparam(), quiet = FALSE, ...) {
  if (quiet) {
    suppressMessages(
      truncateTxome(txdb, maxTxLength = maxTxLength, txEnd = "3prime", 
                    overlapFile = overlapFile, BPPARAM = BPPARAM, ...)
    )
  } else {
    truncateTxome(txdb, maxTxLength = maxTxLength, txEnd = "3prime", 
                  overlapFile = overlapFile, BPPARAM = BPPARAM, ...)
  }
})


#' @rdname truncateTxome
#' @export
setGeneric("truncate5primeTxome", function(txdb, maxTxLength = 300, overlapFile = NULL, 
                                           BPPARAM = bpparam(), quiet = FALSE, ...) {
  standardGeneric("truncate5primeTxome")
}
)

#' @rdname truncateTxome
#' @export
setMethod("truncate5primeTxome", "TxDb", function(txdb, maxTxLength = 300, overlapFile = NULL, 
                                                  BPPARAM = bpparam(), quiet = FALSE, ...) {
  if (quiet) {
    suppressMessages(
      truncateTxome(txdb, maxTxLength = maxTxLength, txEnd = "5prime", 
                    overlapFile = overlapFile, BPPARAM = BPPARAM, ...)
    )
  } else {
    truncateTxome(txdb, maxTxLength = maxTxLength, txEnd = "5prime", 
                  overlapFile = overlapFile, BPPARAM = BPPARAM, ...)
  }
})

setMethod("truncateTxome", "TxDb", function(txdb){})
                                            

#' Truncate Transcriptome
#'
#' Truncate transcripts to a specific maximum length from either the 3' or 5'
#' end, keeping only the terminal portion of each transcript.
#'
#' @param txdb a \code{TxDb} object representing the transcriptome annotation
#' @param maxTxLength the maximum length of transcripts. Defaults to 500 bp
#' @param txEnd transcript truncation end, either `3prime` (default) or
#'   `5prime`.
#' @param overlapFile optional path to export a TSV file containing transcript
#'   overlaps (query_transcript, subject_transcript) post-truncation. If NULL
#'   (default), no file is exported.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether
#'   and how the method should be parallelized.
#' @param quiet suppress progress messages. Only available for
#'   \code{truncate3primeTxome} and \code{truncate5primeTxome}. Defaults to
#'   FALSE
#' @return a \code{TxDb} object
#'
#' @details \code{truncate3primeTxome} and \code{truncate5primeTxome} are
#'   wrappers that call \code{truncateTxome} with \code{txEnd} preset to
#'   "3prime" or "5prime" respectively. They also provide a \code{quiet}
#'   parameter to suppress messages.
#'
#' The function performs the following steps:
#' \enumerate{
#'   \item Truncates each transcript to the specified maximum length from the chosen end
#'   \item Identifies duplicate transcripts (transcripts
#'         with identical coordinates belonging to the same gene) and removes them to avoid redundancy
#'   \item Rebuilds the TxDb object with updated gene, transcript, and exon ranges
#' }
#'
#' @examples
#' library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
#'
#' ## load annotation
#' txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
#'
#' ## restrict to 'chrI' transcripts
#' seqlevels(txdb) <- c("chrI")
#'
#' ## last 500 nts per tx
#' txdb_w500 <- truncateTxome(txdb)
#' txdb_w500
#'
#' ## last 100 nts per tx
#' txdb_w100 <- truncateTxome(txdb, maxTxLength = 100)
#' txdb_w100
#'
#' ## first 500 nts per tx (5' truncation)
#' txdb_5p_w500 <- truncateTxome(txdb, txEnd = "5prime")
#' txdb_5p_w500
#' 
#' ## using convenience wrapper. Same as truncateTxome(..., txEnd = "3prime")
#' txdb_3p <- truncate3primeTxome(txdb, maxTxLength = 500)
#'
#' ## Suppress messages with quiet parameter
#' txdb_quiet <- truncate3primeTxome(txdb, quiet = TRUE)
#'
#' ## Export overlap information
#' txdb_w500 <- truncateTxome(txdb, overlapFile = "overlaps.tsv")
#' 
#' @importFrom GenomicRanges GRangesList GRanges mcols
#' @importFrom GenomicFeatures exonsBy
#' @importFrom txdbmaker makeTxDbFromGRanges
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom AnnotationDbi select taxonomyId
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom methods setMethod
#' @export
#' @rdname truncateTxome
setMethod("truncateTxome", "TxDb", function(txdb,
                                            maxTxLength = 500,
                                            txEnd = "3prime",
                                            overlapFile = NULL,
                                            BPPARAM = bpparam()) {
  ############################################################################
  # Ensure correct values of `txEnd`
  valid_3prime <- c("3", "3'", "3p", "3prime", "3_prime")
  valid_5prime <- c("5", "5'", "5p", "5prime", "5_prime")
  
  if (txEnd %in% valid_3prime) txEnd <- "3prime"
  if (txEnd %in% valid_5prime) txEnd <- "5prime"
  
  if (!txEnd %in% c("3prime", "5prime")) stop("txEnd parameter not valid - only '3prime' or '5prime' parameters are accepted.")
  
  ############################################################################
  # Split exons by transcripts and create a mapping dictionary from
  # transcript_id to gene_id
  grlExons <- exonsBy(txdb, use.names = TRUE)
  dfTxGene <- suppressMessages(select(txdb, keys = names(grlExons), keytype = "TXNAME", columns = "GENEID"))
  mapTxToGene <- setNames(dfTxGene$GENEID, dfTxGene$TXNAME)
  
  ############################################################################
  # Transcript truncation
  message("Truncating transcripts...")
  clipped <- bplapply(grlExons, .clipTranscript,
                      maxTxLength = maxTxLength, txEnd = txEnd,
                      BPPARAM = BPPARAM
  )
  clipped <- GRangesList(clipped)
  message("Done.")
  
  ############################################################################
  # Remove overlapping transcripts
  message("Checking for duplicate transcripts...")
  overlaps <- findOverlaps(clipped,
                           minoverlap = maxTxLength,
                           ignore.strand = FALSE,
                           drop.self = TRUE, drop.redundant = TRUE
  )
  
  ## ensure genes match
  if (length(overlaps) > 0) {
    idx_genes_match <- mapply(function(idx1, idx2) {
      mapTxToGene[names(clipped[idx1])] == mapTxToGene[names(clipped[idx2])]
    }, idx = queryHits(overlaps), idx2 = subjectHits(overlaps))
    overlaps <- overlaps[idx_genes_match]
  }
  
  ## export overlap data.frame
  if (!is.null(overlapFile) && overlapFile != "") {
    ### create overlap_df with the names of the transcripts
    overlap_df <- data.frame(
      query_transcript = names(clipped)[queryHits(overlaps)],
      subject_transcript = names(clipped)[subjectHits(overlaps)],
      stringsAsFactors = FALSE
    )
    
    ### create parent directory
    output_dir <- dirname(overlapFile)
    if (output_dir != "." && !dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    ### write in disk
    write.csv(overlap_df, overlapFile, row.names = F, quote = F)
    message(sprintf("Post-truncation transcript overlaps exported to: %s", overlapFile))
  }
  
  ## get duplicate indices
  duplicates <- unique(queryHits(overlaps))
  if (length(duplicates) > 0) {
    clipped <- clipped[-duplicates]
  }
  message(sprintf("Removed %d duplicates.", length(duplicates)))
  
  ############################################################################
  # Create the final exon ranges
  message("Creating exon ranges...")
  
  ## flatten with tx_id in metadata
  grExons <- unlist(.mutateEach(clipped, transcript_id = names(clipped)))
  names(grExons) <- NULL
  mcols(grExons)["type"] <- "exon"
  
  ## add gene id
  mcols(grExons)["gene_id"] <- mapTxToGene[mcols(grExons)$transcript_id]
  
  ## reindex exon info
  grExons <- sort(grExons)
  mcols(grExons)["exon_id"] <- seq_along(grExons)
  mcols(grExons)["exon_name"] <- NULL
  ## TODO: include `exon_rank`
  
  message("Done.")
  
  ############################################################################
  # Create the final transcript ranges
  message("Creating tx ranges...")
  ## generate transcripts GRanges with clipped bounds
  grTxs <- unlist(GRangesList(bplapply(clipped, .fillReduce,
                                       BPPARAM = BPPARAM)))
  mcols(grTxs)["transcript_id"] <- names(grTxs)
  mcols(grTxs)["type"] <- "transcript"
  
  ## add gene id
  mcols(grTxs)["gene_id"] <- mapTxToGene[grTxs$transcript_id]
  
  message("Done.")
  
  ############################################################################
  # Create the final gene ranges
  message("Creating gene ranges...")
  grGenes <- unlist(GRangesList(bplapply(split(grTxs, grTxs$gene_id),
                                         .fillReduce,
                                         BPPARAM = BPPARAM
  )))
  mcols(grGenes)["gene_id"] <- names(grGenes)
  mcols(grGenes)["type"] <- "gene"
  message("Done.")
  
  ############################################################################
  # Generate the final TxDb object
  dfMetadata <- data.frame(
    name=c("Truncated by", "Maximum Transcript Length", "Truncation End"),
    value=c("txendcutr", maxTxLength, txEnd)
  )
  
  .suppressTxDbGenomeWarning(
    makeTxDbFromGRanges(c(grGenes, grTxs, grExons),
                        taxonomyId = taxonomyId(txdb),
                        metadata = dfMetadata
    )
  )
})

#' Clip Transcript to Given Length
#'
#' Internal function for operating on individual \code{GRanges}, where ranges
#' represent exons in a transcript. This is designed to be used in an
#' \code{*apply} function over a \code{GRangesList} object.
#'
#' @param gr a \code{GRanges} object
#' @param maxTxLength a positive integer
#' @param txEnd transcript truncation end
#' @param BPPARAM A \linkS4class{BiocParallelParam} object for parallelization
#'
#' @return the clipped \code{GRanges} object
#'
#' @importFrom GenomicRanges GRanges width strand start end intersect invertStrand
#' @importFrom IRanges IRanges
#'
.clipTranscript <- function(gr, maxTxLength, txEnd) {
  if (sum(width(gr)) <= maxTxLength) { ## already short enough
    gr
  } else { ## need to adjust
    ## adjustment is directed
    txStrand <- strand(gr)
    virtual_txStrand <- if (txEnd == "3prime") txStrand else invertStrand(txStrand)

    if (all(virtual_txStrand == "+")) {
      ## order txs
      idx <- order(-end(gr))

      ## compute cumulative lengths
      cumLength <- cumsum(width(gr[idx]))

      ## index of exon that exceeds maximum length
      idxLast <- min(which(cumLength > maxTxLength))

      ## compute cutoff (genomic position)
      startNew <- start(gr[idx[idxLast]]) + (cumLength[idxLast] - maxTxLength)

      ## new transcript interval
      grMask <- GRanges(seqnames(gr[1]),
        IRanges(startNew, max(end(gr))),
        strand = "+"
      )

      if (txEnd == "5prime") grMask <- invertStrand(grMask)

      ## clip exons with interval
      intersect(gr, grMask)
    } else if (all(virtual_txStrand == "-")) {
      ## order txs
      idx <- order(start(gr))

      ## compute cumulative lengths
      cumLength <- cumsum(width(gr[idx]))

      ## index of exon that exceeds maximum length
      idxLast <- min(which(cumLength > maxTxLength))

      ## compute cutoff (genomic position)
      endNew <- end(gr[idx[idxLast]]) - (cumLength[idxLast] - maxTxLength)

      ## new transcript interval
      grMask <- GRanges(seqnames(gr[1]),
        IRanges(min(start(gr)), endNew),
        strand = "-"
      )

      if (txEnd == "5prime") grMask <- invertStrand(grMask)

      ## clip exons with interval
      intersect(gr, grMask)
    } else {
      warning("Skipping Transcript: Encountered inconsistent strand annotation!", gr)
      gr
    }
  }
}


#' Convert GRanges to Single Range
#'
#' @param gr a \code{GRanges} with ranges to be merged.
#' @param validate logical determining whether entries should be checked for compatible
#' seqnames and strands.
#'
#' @return \code{GRanges} with single interval
#'
#' @details The validation assumes seqnames and strand are \code{Rle} objects.
#'
#' @importFrom GenomicRanges seqnames start end strand reduce start<- end<-
#' @importFrom S4Vectors nrun
.fillReduce <- function(gr, validate = TRUE) {
  if (validate) {
    stopifnot(
      nrun(seqnames(gr)) == 1,
      nrun(strand(gr)) == 1
    )
  }

  ## TODO: Check if faster to construct new GRanges
  ## Current implementation makes retention of seqinfo simple.
  start(gr) <- min(start(gr))
  end(gr) <- max(end(gr))
  reduce(gr)
}
