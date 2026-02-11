#' @importFrom methods setGeneric
setGeneric(".mutateEach", signature=c("grl"),
           function(grl, ...) standardGeneric(".mutateEach"))

#' Efficient Metadata Columns Mutation
#'
#' @param grl a CompressedGRangesList
#' @param ... named list of vectors to insert as metadata columns on each element
#' \code{GRanges}. Each vector length must match the length of the \code{GRangesList}.
#'
#' @return a \code{CompressedGRangesList} with all element \code{GRanges} updated
#' with supplied metadata columns
#'
#' @importFrom methods setMethod slot slot<-
#' @importFrom utils capture.output
#' @importFrom S4Vectors elementNROWS
setMethod(".mutateEach", "CompressedGRangesList",
          function (grl, ...) {
              ## ensure data is a valid length
              inputLengths <- vapply(list(...), length, integer(1), USE.NAMES=TRUE)
              if (!all(inputLengths == length(grl))) {
                stop("Mismatched lengths detected:\n",
                     "\tExpected length ", length(grl),
                     ", but found length(s) ",
                     capture.output(dput(inputLengths)))
              }

              ## expand list to full length
              expandedListData <- lapply(list(...), rep, times=elementNROWS(grl))
              mcols(slot(grl, "unlistData"))[names(list(...))] <- expandedListData
              grl
          }
)

#' @importFrom methods setMethod
setMethod(".mutateEach", "SimpleGRangesList",
          function (grl, ...) .mutateEach(GRangesList(grl, compress=TRUE), ...)
)

#' Suppress `txdbmaker` warning message
#'
#' Due to missing genome information in the input data,
#' \code{txdbmaker::makeTxDbFromGRanges} always returns a warning message
#' stating 'genome version information is not available for this TxDb object'.
#' This warning does not affect functionality and can be safely suppressed. This
#' function executes a TxDb-generating expression while suppressing only this
#' specific warning.
#'
#' @param fn_call an R expression that creates a TxDb object, expected a call to
#'   \code{makeTxDbFromGRanges} or similar function.
#'
#' @return a \code{txdb} object representing a transcriptome, as returned by the
#'   expression passed to \code{fn_call}.
.suppressTxDbGenomeWarning <- function(fn_call){
  withCallingHandlers(
    {
      tmp_txdb <- fn_call
    },
    warning = function(w){
      if(grepl("genome version information is not available", conditionMessage(w))){
        invokeRestart("muffleWarning")
      }
    }
  )
  
  return(tmp_txdb)
}