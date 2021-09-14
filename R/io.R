#' Export GTF
#'
#' Exports a TxDb annotation to a GTF file
#'
#' @param txdb transcriptome to be output
#' @param file location to write GTF
#' @param source a string to go in the \code{source} column
#' @return The \code{txdb} argument is invisibly returned.
#'
#' @examples
#' library(GenomicFeatures)
#'
#' ## create TxDb using SGD genes for sacCer3 genome from UCSC
#' txdb <- makeTxDbFromUCSC("sacCer3", "sgdGene")
#'
#' ## restrict to 'chrI' transcripts
#' seqlevels(txdb) <- c("chrI")
#'
#' ## last 500 nts per tx
#' txdb_w500 <- truncateTxome(txdb)
#'
#' exportGTF(txdb_w500, "sacCer3.sgdGene.w500.gtf")
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
#' @param file output FASTA file
#' @return The \code{txdb} argument is invisibly returned.
#'
#' @examples
#' library(GenomicFeatures)
#' library(BSGenome)
#'
#' ## create TxDb using SGD genes for sacCer3 genome from UCSC
#' txdb <- makeTxDbFromUCSC("sacCer3", "sgdGene")
#'
#' ## retrieve sacCer3 sequence as BSGenome object
#' sacCer3 <- getBSgenome("sacCer3")
#'
#' ## restrict to 'chrI' transcripts (makes for briefer example runtime)
#' seqlevels(txdb) <- c("chrI")
#'
#' ## last 500 nts per tx
#' txdb_w500 <- truncateTxome(txdb)
#'
#' exportFASTA(txdb_w500, sacCer3, "sacCer3.sgdGene.w500.fa")
#'
#' @importFrom GenomicFeatures extractTranscriptSeqs
#' @importFrom Biostrings writeXStringSet
#' @export
exportFASTA <- function (txdb, genome, file) {
    seqs <- extractTranscriptSeqs(genome, txdb, use.names=TRUE)
    writeXStringSet(seqs, file, format="fasta")
    invisible(txdb)
}
