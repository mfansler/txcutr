#' Export GTF
#'
#' Exports a TxDb annotation to a GTF file
#'
#' @param txdb transcriptome to be output
#' @param file a string or \code{\link{connection}} to output GTF file.
#'   Automatically recognizes strings ending with "**.gz**" for zipped output.
#' @param source a string to go in the \code{source} column
#' @return The \code{txdb} argument is invisibly returned.
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
#'
#' ## export uncompressed
#' exportGTF(txdb_w500, "sacCer3.sgdGene.w500.gtf")
#'
#' ## export compressed
#' exportGTF(txdb_w500, "sacCer3.sgdGene.w500.gtf.gz")
#'
#' @importFrom rtracklayer export
#' @importFrom BiocGenerics which paste
#' @importFrom S4Vectors mcols<-
#' @export
exportGTF <- function (txdb, file, source="txcutr") {
    grl <- txdbToGRangesList(txdb)

    ## update gene fields
    mcols(grl$genes)[["source"]] <- source

    ## update transcript fields
    mcols(grl$transcripts)[["source"]] <- source
    mcols(grl$transcripts)[["transcript_id"]] <- grl$transcripts$tx_name
    grl$transcripts$tx_name <- NULL

    ## update exon fields
    mcols(grl$exons)[["source"]] <- source
    mcols(grl$exons)[["transcript_id"]] <- grl$exons$tx_name
    mcols(grl$exons)[["exon_id"]] <- paste(grl$exons$tx_name,
                                           grl$exons$exon_rank, sep=":")
    grl$exons$tx_name <- NULL
    grl$exons$exon_rank <- NULL

    idxFwd <- which(strand(grl$exons) == "+")
    idxRev <- which(strand(grl$exons) == "-")

    ## export
    if (!inherits(file, "connection")) {
      if (grepl(".gz$", file)) {
        file <- gzfile(file, "wab")
        on.exit(close(file))
      }
    }
    export(grl$genes, file, format="GTF")
    export(grl$transcripts, file, format="GTF", append=TRUE)
    if (length(idxFwd) > 0) {
      export(grl$exons[idxFwd], file, format="GTF", append=TRUE)
    }
    if (length(idxRev) > 0) {
      export(sort(grl$exons[idxRev], decreasing=TRUE), file, format="GTF",
             append=TRUE)
    }

    invisible(txdb)
}


#' Export Transcriptome as FASTA
#'
#' @param txdb a \code{TxDb} object representing a transcriptome annotation
#' @param genome a \code{BSgenome} object from which to extract sequences
#' @param file a string for output FASTA file. File names ending in "**.gz**"
#'     will automatically use gzip compression.
#' @param ... additional arguments to pass through to
#'     \code{\link[Biostrings]{writeXStringSet}}
#' @return The \code{txdb} argument is invisibly returned.
#'
#' @examples
#' library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
#' library(BSgenome.Scerevisiae.UCSC.sacCer3)
#'
#' ## load annotation and genome
#' txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
#' sacCer3 <- BSgenome.Scerevisiae.UCSC.sacCer3
#'
#' ## restrict to 'chrI' transcripts (makes for briefer example runtime)
#' seqlevels(txdb) <- c("chrI")
#'
#' ## last 500 nts per tx
#' txdb_w500 <- truncateTxome(txdb)
#'
#' ## export uncompressed
#' exportFASTA(txdb_w500, sacCer3, "sacCer3.sgdGene.w500.fa")
#'
#' ## export compressed
#' exportFASTA(txdb_w500, sacCer3, "sacCer3.sgdGene.w500.fa.gz")
#'
#' @importFrom GenomicFeatures extractTranscriptSeqs
#' @importFrom Biostrings writeXStringSet
#' @export
exportFASTA <- function (txdb, genome, file, ...) {
    seqs <- extractTranscriptSeqs(genome, txdb, use.names=TRUE)
    compress <- grepl(".gz$", file)
    writeXStringSet(seqs, file, format="fasta", compress=compress, ...)
    invisible(txdb)
}

#' Export Merge Table for Transcriptome
#'
#' @param txdb a \code{TxDb} object representing a transcriptome annotation
#' @param file a string or \code{\link{connection}} to output TSV file.
#'   Automatically recognizes strings ending with "**.gz**" for zipped output.
#' @param minDistance the minimum separation to regard overlapping transcripts
#'     as unique.
#' @return The \code{txdb} argument is invisibly returned.
#'
#' @examples
#' library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
#'
#' ## load annotation
#' txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
#'
#' ## restrict to 'chrI' transcripts (makes for briefer example runtime)
#' seqlevels(txdb) <- c("chrI")
#'
#' ## last 500 nts per tx
#' txdb_w500 <- truncateTxome(txdb)
#'
#' ## export plain format
#' exportMergeTable(txdb_w500, "sacCer3.sgdGene.w500.merge.tsv")
#'
#' ## export compressed format
#' exportMergeTable(txdb_w500, "sacCer3.sgdGene.w500.merge.tsv.gz")
#'
#' @importFrom utils write.table
#' @export
exportMergeTable <- function (txdb, file, minDistance=200L) {
  df <- generateMergeTable(txdb, minDistance=minDistance)
  if (!inherits(file, "connection")) {
    if (grepl(".gz$", file)) {
      file <- gzfile(file, "wb")
      on.exit(close(file))
    }
  }
  write.table(df, file, sep="\t", row.names=FALSE, quote=FALSE)
  invisible(txdb)
}
