#' Non-overlapping temporal aggregation of a time series
#'
#' Non-overlapping temporal aggregation of a time series according
#' to a specific aggregation order.
#'
#' @param agg_order Aggregation order to consider.
#' @param x Univariate time series: a vector or a \code{ts} object.
#' @param align Specifies whether the aggregates should be aligned with the start ()
#' or end of the series.
#' @param rm_na logical. Should missing values be removed?
#'
#' @return A vector or \code{ts} object
#'
#' @family utilities
#' @importFrom stats is.ts na.omit ts tsp
#'
#' @examples
#' data(FoReco_data)
#' annual_ts <- agg_ts(12, FoReco_data$obs$k1[,1]) # == FoReco_data$obs$k12[,1]
#'
#' @export
agg_ts <- function(agg_order, x, align = "end", rm_na = FALSE){
  if(is.ts(x)){
    tspx <- tsp(x)
    x <- as.matrix(x)
    isvec <- NCOL(x)==1
  }else{
    if(is.vector(x)){
      x <- cbind(x)
    }
    tspx <- NULL
  }

  align <- match.arg(align, c("end","start"))
  n <- NROW(x)
  agg_order <- as.integer(agg_order)

  if(align=='end'){
    start <- n%%agg_order + 1L
  }else{
    start <- 1L
  }

  nk <- trunc(n/agg_order)
  out <- apply(x, 2, function(col){
    tmp <- matrix(col[start - 1L + seq_len(agg_order*nk)], ncol=nk)
    colSums(tmp)
  }, simplify = FALSE)
  out <- do.call(cbind, out)

  if(align=='end' & n%%agg_order != 0){
    out <- rbind(NA, out)
  }else if(align=='start' & n%%agg_order != 0){
    out <- rbind(out, NA)
  }

  if(NCOL(out)==1){
    out <- out[,]
  }


  if(!is.null(tspx)){
    out <- ts(out, frequency=tspx[3]/agg_order,
              start=tspx[1]-1+tspx[3]/agg_order)
  }

  if(rm_na){
    out <- na.omit(out)
  }

  return(out)
}
