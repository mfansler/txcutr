#' Convert TxDb object to GRangesList
#'
#' @rdname txdbToGRangesList
#'
#' @param txdb a \code{TxDb} object
#' @param geneCols names of columns to include in the \code{genes} ranges
#' @param transcriptCols names of columns to include in the \code{transcripts} ranges
#' @param exonCols names of columns to include in the \code{exons} ranges
#'
#' @return a \code{GRangesList} object with entries \code{c(genes, transcripts, exons)}
#'
#' @examples
#' library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
#'
#' ## load annotation
#' txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
#'
#' grl <- txdbToGRangesList(txdb)
#' grl
#'
#' @importFrom GenomicRanges mcols GRangesList
#' @importFrom GenomicFeatures exons transcripts genes
#' @importFrom S4Vectors expand mcols<-
#' @export
txdbToGRangesList <- function (txdb,
                               geneCols=c("gene_id"),
                               transcriptCols=c("gene_id", "tx_name"),
                               exonCols=c("gene_id", "tx_name",
                                          "exon_id", "exon_rank")) {

    grExons <- exons(txdb, columns=exonCols)
    mcols(grExons)["type"] <- "exon"
    grExons <- expand(grExons)

    grTranscripts <- transcripts(txdb, columns=transcriptCols)
    mcols(grTranscripts)["type"] <- "transcript"
    grTranscripts <- expand(grTranscripts)

    grGenes <- genes(txdb, columns=geneCols)
    mcols(grGenes)["type"] <- "gene"

    GRangesList(genes=grGenes, transcripts=grTranscripts, exons=grExons,
                compress=FALSE)
}
