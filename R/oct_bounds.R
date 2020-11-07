# Title
#
# @param hts_bounds
# @param thf_bounds
# @param m
# @param C
# @param Ut
# @param nb
#
# @return
# @export
#
# @examples
# hts_bounds <- matrix(c(rnorm(3), rnorm(3, 100)), nrow = 3)
# hts_bounds[1,2] <- Inf
# hts_bounds[2,1] <- -Inf
# thf_bounds <- matrix(c(rnorm(7), rnorm(7, 100)), nrow = 7)
# obj <- oct_bounds(hts_bounds = hts_bounds, m = c(4,1), C = matrix(1, ncol = 2))
oct_bounds <- function(hts_bounds, thf_bounds, m, C, Ut, nb){
  if(missing(thf_bounds) & missing(hts_bounds)){
    stop("ciao")
  }else if(missing(thf_bounds)){
    hts_nrow <- NROW(hts_bounds)
    thf_nrow <- thf_tools(m = m)$kt

    if(!missing(C)){
      n <- hts_tools(C = C)$n
      if(hts_nrow != n){
        stop("hts_bounds must be a matrix (", n, " x 2)", call. = FALSE)
      }
    }else if(!(missing(Ut) & missing(nb))){
      n <- hts_tools(Ut = Ut, nb = nb)$n
      if(hts_nrow != n){
        stop("hts_bounds must be a matrix (", n, " x 2)", call. = FALSE)
      }
    }

    bhts <- t(simplify2array(rep(split(hts_bounds, 1:hts_nrow), each = thf_nrow)))
    dimnames(bhts) <- NULL
    colnames(bhts) <- c("lower", "upper")
    return(bhts)
  }else if(missing(hts_bounds)){
    if(!missing(C)){
      hts_nrow <- hts_tools(C = C)$n
    }else if(!(missing(Ut) & missing(nb))){
      hts_nrow <- hts_tools(Ut = Ut, nb = nb)$n
    }else{
      stop("the argument C (or Ut AND nb) is not specified", call. = FALSE)
    }

    thf_nrow <- NROW(thf_bounds)

    if(!missing(m)){
      kt <- thf_tools(m = m)$kt
      if(thf_nrow != kt){
        stop("thf_bounds must be a matrix (", kt, " x 2)", call. = FALSE)
      }
    }

    bthf <- t(simplify2array(rep(split(thf_bounds, 1:thf_nrow), hts_nrow)))
    dimnames(bthf) <- NULL
    colnames(bthf) <- c("lower", "upper")
    return(bthf)
  }else{
    kt <- thf_tools(m = m)$kt

    if(!missing(C)){
      n <- hts_tools(C = C)$n
    }else if(!(missing(Ut) & missing(nb))){
      n <- hts_tools(Ut = Ut, nb = nb)$n
    }else{
      stop("the argument C (or Ut AND nb) is not specified", call. = FALSE)
    }
    thf_nrow <- NROW(thf_bounds)
    hts_nrow <- NROW(hts_bounds)

    if(hts_nrow != n){
      stop("hts_bounds must be a matrix (", n, " x 2)", call. = FALSE)
    }

    if(thf_nrow != kt){
      stop("thf_bounds must be a matrix (", kt, " x 2)", call. = FALSE)
    }

    bthf <- t(simplify2array(rep(split(thf_bounds, 1:thf_nrow), hts_nrow)))
    bhts <- t(simplify2array(rep(split(hts_bounds, 1:hts_nrow), each = thf_nrow)))
    out <- cbind(apply(cbind(bhts[,1], bthf[,1]), 1, max),
                 apply(cbind(bhts[,2], bthf[,2]), 1, min))
    dimnames(out) <- NULL
    colnames(out) <- c("lower", "upper")
    return(out)
  }
}

