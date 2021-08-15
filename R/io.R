#' Export GTF
#'
#' Exports a TxDb annotation to a GTF file
#'
#' @param txdb transcriptome to be output
#' @param file location to write GTF
#' @param source a string to go in the \code{source} column
#'
#' @export
#'
#' @importFrom rtracklayer export
#' @importFrom BiocGenerics which paste
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
  mcols(grl$exons)[["exon_id"]] <- paste(grl$exons$tx_name, grl$exons$exon_rank, sep=":")
  grl$exons$tx_name <- NULL
  grl$exons$exon_rank <- NULL

  idxFwd <- which(strand(grl$exons) == "+")
  idxRev <- which(strand(grl$exons) == "-")

  ## export
  export(grl$genes, file, format="GTF")
  export(grl$transcripts, file, format="GTF", append=TRUE)
  if (length(idxFwd) > 0) {
    export(grl$exons[idxFwd], file, format="GTF", append=TRUE)
  }
  if (length(idxRev) > 0) {
    export(sort(grl$exons[idxRev], decreasing=TRUE), file, format="GTF", append=TRUE)
  }
}


#' Export Transcriptome as FASTA
#'
#' @param txdb a \code{TxDb} object representing a transcriptome annotation
#' @param genome a \code{BSgenome} object from which to extract sequences
#' @param file output FASTA file
#'
#' @export
#'
#' @importFrom GenomicFeatures extractTranscriptSeqs
#' @importFrom Biostrings writeXStringSet
exportFASTA <- function (txdb, genome, file) {
  seqs <- extractTranscriptSeqs(genome, txdb, use.names=TRUE)
  writeXStringSet(seqs, file, format="fasta")
}
