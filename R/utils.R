#----------------------------------------------------------------------------------------------------
#' provide formatted console output as in the standard C library
#'
#' @param ... format string and corresponding variables
#'
#' @return no return value, side effect is a formatted string printed to the R console
#'
#'
#' @aliases printf
#' @rdname printf
#'
#' @export
#'
printf <- function(...)
{
   print(noquote(sprintf(...)))

} # printf
#----------------------------------------------------------------------------------------------------
#' truncate precision in data.frames
#'
#' column by column rounding of data.frame numerical values, useful for readibility
#'
#' @param tbl  a data.frame
#' @param digits the number of places of decimal precision to preserve
#' @param exponentalColumnNames so that they can be preserved in exponential notation
#'
#' @return the transformed data.frame, with rows and columns intact
#'
#' @export
roundNumericColumns <- function(tbl, digits, exponentialColumnNames=NA)
{
  tbl.exponentials <- data.frame()
  tbl.main <- tbl

  if(!(all(is.na(exponentialColumnNames)))){
     exponential.cols <- grep(exponentialColumnNames, colnames(tbl))
     stopifnot(length(exponential.cols) == length(exponentialColumnNames))
     tbl.exponentials <- tbl[, exponential.cols, drop=FALSE]
     tbl.main <- tbl[, -exponential.cols, drop=FALSE]
     }

  numeric_columns <- sapply(tbl.main, mode) == 'numeric'
  tbl.main[numeric_columns] <-  round(tbl.main[numeric_columns], digits)

  if(ncol(tbl.exponentials) > 0){
     tbl.exponentials <- apply(tbl.exponentials, 2, function(col) as.numeric(formatC(col, format = "e", digits = 2)))
     }

  tbl.out <- cbind(tbl.main, tbl.exponentials)[, colnames(tbl)]
  tbl.out

} # roundNumericColumns
#------------------------------------------------------------------------------------------------------------------------


