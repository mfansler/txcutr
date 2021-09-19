#' @rdname truncateTxome
#' @param txdb an object representing a transcriptome
#' @param maxTxLength the maximum length of resulting transcripts
#' @param ... additional arguments
#'
#' @return a \code{TxDb} object
#' @export
#'
#' @importFrom methods setGeneric
setGeneric("truncateTxome", signature=c("txdb", "maxTxLength"),
           function(txdb, maxTxLength=500, ...) standardGeneric("truncateTxome")
)

#' Truncate Transcriptome
#'
#' @rdname truncateTxome
#'
#' @param txdb a \code{TxDb} object
#' @param maxTxLength the maximum length of transcripts
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether
#'   and how the method should be parallelized.
#' @return a \code{TxDb} object
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
#' txdb_w100 <- truncateTxome(txdb, maxTxLength=100)
#' txdb_w100
#'
#' @importFrom GenomicRanges GRangesList mcols
#' @importFrom GenomicFeatures exonsBy
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom AnnotationDbi select
#' @importFrom methods setMethod
#' @export
setMethod("truncateTxome", "TxDb", function(txdb,
                                            maxTxLength=500,
                                            BPPARAM=bpparam()) {
    grlExons <- exonsBy(txdb, use.names=TRUE)
    dfTxGene <- select(txdb, keys=names(grlExons),
                       keytype="TXNAME", columns="GENEID")
    mapTxToGene <- setNames(dfTxGene$GENEID, dfTxGene$TXNAME)

    message("Truncating transcripts...")
    clipped <- bplapply(grlExons, .clipTranscript, maxTxLength=maxTxLength,
                        BPPARAM=BPPARAM)
    clipped <- GRangesList(clipped)
    message("Done.")

    message("Checking for duplicate transcripts...")
    overlaps <- findOverlaps(clipped, minoverlap=maxTxLength,
                             ignore.strand=FALSE,
                             drop.self=TRUE, drop.redundant=TRUE)
    ## ensure genes match
    if (length(overlaps) > 0) {
      idx_genes_match <- mapply(function (idx1, idx2) {
        mapTxToGene[names(clipped[idx1])] == mapTxToGene[names(clipped[idx2])]
      }, idx=queryHits(overlaps), idx2=subjectHits(overlaps))
      overlaps <- overlaps[idx_genes_match]
    }

    ## get duplicate indices
    duplicates <- unique(queryHits(overlaps))
    if (length(duplicates) > 0) {
        clipped <- clipped[-duplicates]
    }
    message(sprintf("Removed %d duplicates.", length(duplicates)))

    message("Creating exon ranges...")
    ## flatten with tx_id in metadata
    grExons <- unlist(.mutateEach(clipped, transcript_id=names(clipped)))
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

    message("Creating tx ranges...")
    ## generate transcripts GRanges with clipped bounds
    grTxs <- unlist(GRangesList(bplapply(clipped, .fillReduce,
                                         BPPARAM=BPPARAM)))
    mcols(grTxs)['transcript_id'] <- names(grTxs)
    mcols(grTxs)["type"] <- "transcript"

    ## add gene id
    mcols(grTxs)["gene_id"] <- mapTxToGene[grTxs$transcript_id]

    message("Done.")

    message("Creating gene ranges...")
    grGenes <- unlist(GRangesList(bplapply(split(grTxs, grTxs$gene_id),
                                           .fillReduce, BPPARAM=BPPARAM)))
    mcols(grGenes)['gene_id'] <- names(grGenes)
    mcols(grGenes)["type"] <- "gene"
    message("Done.")

    dfMetadata <- data.frame(
      name=c("Truncated by", "Maximum Transcript Length"),
      value=c("txcutr", maxTxLength)
    )

    makeTxDbFromGRanges(c(grGenes, grTxs, grExons),
                        taxonomyId=taxonomyId(txdb),
                        metadata=dfMetadata)
}
)

#' Clip Transcript to Given Length
#'
#' Internal function for operating on individual \code{GRanges}, where ranges
#' represent exons in a transcript. This is designed to be used in an
#' \code{*apply} function over a \code{GRangesList} object.
#'
#' @param gr a \code{GRanges} object
#' @param maxTxLength a positive integer
#'
#' @return the clipped \code{GRanges} object
#'
#' @importFrom GenomicRanges GRanges width strand start end intersect
#' @importFrom IRanges IRanges
#'
.clipTranscript <- function (gr, maxTxLength) {
    if (sum(width(gr)) <= maxTxLength) { ## already short enough
      gr
    } else { ## need to adjust
        ## adjustment is directed
        txStrand <- strand(gr)
        if (all(txStrand == "+")) {
          ## order by 3' ends
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
                            strand="+")

          ## clip exons with interval
          intersect(gr, grMask)
      } else if (all(txStrand == "-")) {
          ## order by 3' ends
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
                            strand="-")

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
.fillReduce <- function (gr, validate=TRUE) {
  if (validate) {
    stopifnot(nrun(seqnames(gr)) ==  1,
              nrun(strand(gr)) ==  1)
  }

  ## TODO: Check if faster to construct new GRanges
  ## Current implementation makes retention of seqinfo simple.
  start(gr) <- min(start(gr))
  end(gr) <- max(end(gr))
  reduce(gr)
}

