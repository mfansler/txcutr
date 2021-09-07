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
#' @importFrom methods setMethod
setMethod(".mutateEach", "CompressedGRangesList",
          function (grl, ...) {
            ## ensure data is a valid length
            inputLengths <- vapply(list(...), length, integer(1), USE.NAMES=TRUE)
            if (!all(inputLengths == length(grl))) {
              stop("Mismatched lengths detected:\n",
                   "\tExpected length ", length(grl),
                   ", but found length(s) ", capture.output(dput(inputLengths)))
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

