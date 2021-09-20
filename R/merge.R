#' @rdname generateMergeTable
#' @param txdb an object representing a transcriptome
#' @param minDistance the minimum separation to regard overlapping transcripts
#'     as unique
#'
#' @return a \code{data.frame} with three columns
#'     - \code{tx_in} the input transcript
#'     - \code{tx_out} the transcript merged into
#'     - \code{gene_out} the gene merged into
#'
#' @importFrom methods setGeneric
#' @export
setGeneric("generateMergeTable", signature=c("txdb", "minDistance"),
           function(txdb, minDistance=200) standardGeneric("generateMergeTable")
)

#' Generate Merge Table
#'
#' @rdname generateMergeTable
#' @param txdb an object representing a transcriptome
#' @param minDistance the minimum separation to regard overlapping transcripts
#'     as unique
#'
#' @return a \code{data.frame} with three columns
#'     - \code{tx_in} the input transcript
#'     - \code{tx_out} the transcript merged into
#'     - \code{gene_out} the gene merged into
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
#' @importFrom GenomicRanges mcols resize findOverlaps
#' @importFrom GenomicFeatures transcripts
#' @importFrom methods setMethod
#' @importFrom stats setNames
#' @export
setMethod("generateMergeTable", "TxDb", function(txdb, minDistance=200L) {
    grTxs <- transcripts(txdb,
                         columns=c("gene_id", "tx_id", "tx_name"))

    ## use only ends
    grTxs <- resize(grTxs, width=minDistance, fix="end", ignore.strand=FALSE)

    ## generate raw overlaps
    overlaps <- as.data.frame(findOverlaps(grTxs, ignore.strand=FALSE,
                                           drop.self=TRUE))
    overlaps["tx_in"] <- unlist(grTxs$tx_name[overlaps$queryHits])
    overlaps["tx_out"] <- unlist(grTxs$tx_name[overlaps$subjectHits])
    overlaps["gene_in"] <- unlist(grTxs$gene_id[overlaps$queryHits])
    overlaps["gene_out"] <- unlist(grTxs$gene_id[overlaps$subjectHits])

    ## filter unmatched genes
    overlaps <- overlaps[overlaps$gene_in == overlaps$gene_out,]

    ## keep downstream
    overlaps["strand"] <- as.character(strand(grTxs))[overlaps$queryHits]
    overlaps["end_in"] <- ifelse(overlaps$strand == "+",
                                 end(grTxs[overlaps$queryHits]),
                                 -start(grTxs[overlaps$queryHits]))
    overlaps["end_out"] <- ifelse(overlaps$strand == "+",
                                  end(grTxs[overlaps$subjectHits]),
                                  -start(grTxs[overlaps$subjectHits]))
    overlaps <- overlaps[overlaps$end_in <= overlaps$end_out,]

    ## if equal, prioritize early transcript
    ## e.g., a transcript ENSMUST000000001 will replace ENSMUST000001000
    idxEqual <- which(overlaps$end_in == overlaps$end_out)
    if (length(idxEqual) > 0) {
        isLaterTx <- overlaps[idxEqual, "tx_in"] < overlaps[idxEqual, "tx_out"]
        overlaps <- overlaps[-idxEqual[isLaterTx], ]
    }

    ## pick unique out
    groupedOverlaps <- split(overlaps, overlaps$queryHits)
    singleOverlaps <- lapply(groupedOverlaps,
                             function (df) {
                                 df[with(df, order(-end_out, tx_out)),][1,]
                             })
    dfMerge <- do.call(rbind, singleOverlaps)[, c("tx_in", "tx_out")]
    dfMerge <- .propagateMap(dfMerge)

    ## include self-maps
    unmergedTxs <- grTxs$tx_name[!(grTxs$tx_name %in% dfMerge$tx_in)]
    dfMerge <- rbind(dfMerge, data.frame(tx_in=unmergedTxs, tx_out=unmergedTxs))

    ## append gene
    txToGeneMap <- setNames(object=unlist(grTxs$gene_id), nm=grTxs$tx_name)
    dfMerge['gene_out'] <- txToGeneMap[as.character(dfMerge$tx_out)]

    ## reorder and drop rownames
    dfMerge <- dfMerge[with(dfMerge, order(tx_in)),]
    rownames(dfMerge) <- NULL

    dfMerge
}
)


#' Propagate Transcript Merge Map
#'
#' @param df a \code{data.frame} with columns \code{tx_in} and \code{tx_out}
#' @param MAXITERS a numeric controlling the maximum number of iterations
#'
#' @return a converged \code{data.frame}, such that, \code{tx_out} is not
#'     present in any \code{tx_in}
#'
#' @importFrom stats setNames
.propagateMap <- function (df, MAXITERS=1000) {
    ## NOTE: This method assumes no loops. To provide execution safety, we
    ## enforce a maximum number of iterations.
    iteration <- 0

    ## while there are any outputs in the inputs,
    while(any(df$tx_out %in% df$tx_in) && iteration < MAXITERS) {
        iteration <- iteration + 1

        ## generate merge map for simple translation
        mergeMap <- setNames(object=df$tx_out, nm=df$tx_in)

        ## these will be propagated
        idxMappable <- df$tx_out %in% df$tx_in

        ## outputs to translate
        oldOut <- as.character(df$tx_out[idxMappable])

        ## translations
        newOut <- mergeMap[oldOut]

        ## replace old with new
        df[idxMappable, "tx_out"] <- newOut
    }
    ## did we fail to converge?
    stopifnot(iteration <= MAXITERS && !any(df$tx_out %in% df$tx_in))

    df
}
