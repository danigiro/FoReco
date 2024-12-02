# Print pretty number
fnumber <- function(x, m = 30){
  sapply(x, function(z){
    fulln <- formatC(z, format = "f", digits = 2)
    format <- ifelse((z<1 & z !=0) | nchar(fulln)>m, "e", "f")
    formatC(z, width = m,
            format = format, digits = 2)
  })
}

# Split cross-temporal matrix in a temporal list
mat2list <- function(mat, kset){
  m <- max(kset)
  h <- NCOL(mat)/sum(kset)
  kid <- rep(kset, h*m/kset)
  split.data.frame(t(mat), kid)[as.character(kset)]
}

# Trasform the cross-temporal [n x (h*kt)] matrix into a [h x (n*kt)] matrix
# See also: hmat2mat()
mat2hmat <- function(mat, h, kset, n){
  m <- max(kset)
  i <- rep(rep(rep(1:h, length(kset)), rep(m/kset, each = h)), n)
  vec <- as.vector(t(mat))
  matrix(vec[order(i)], nrow = h, byrow = T)
}

# Trasform the [h x (n*kt)] matrix into a cross-temporal [n x (h*kt)] matrix
# See also: mat2hmat()
hmat2mat <- function(hmat, h, kset, n){
  m <- max(kset)
  i <- rep(1:sum(m/kset), h*n)
  it <- rep(rep(m/kset, m/kset), h*n)
  ih <- rep(1:h, each = n*sum(m/kset))
  out <- matrix(as.vector(t(hmat))[order(it, ih, i)], nrow = n)
  colnames(out) <- namesTE(kset = kset, h = h)
  out
}

# Trasform the temporal [(h*kt) x 1] vector into a [h x kt] matrix
# See also: hmat2vec()
vec2hmat <- function(vec, h, kset){
  m <- max(kset)
  i <- rep(rep(1:h, length(kset)), rep(m/kset, each = h))
  matrix(vec[order(i)], nrow = h, byrow = T)
}

# Trasform the [h x kt] matrix into a temporal [(h*kt) x 1] vector
# See also: mat2hmat()
hmat2vec <- function(hmat, h, kset){
  m <- max(kset)
  i <- rep(1:sum(m/kset), h)
  it <- rep(rep(m/kset, m/kset), h)
  ih <- rep(1:h, each = sum(m/kset))
  out <- as.vector(t(hmat))[order(it, ih, i)]
  names_vec <- namesTE(kset = kset, h = h)
  setNames(out, names_vec)
}

# Build a named vector to specify k and h position
namesTE <- function(kset, h){
  m <- max(kset)
  seqk <- h * (m/kset)
  paste0("k-", rep(kset, seqk),
         " h-", Reduce("c", sapply(seqk, seq_len)))
}

# Build a cross-sectional name
namesCS <- function(n, names_vec = NULL, names_list = NULL){
  if(!is.null(names_vec)){
    return(names_vec)
  }else if(length(names_list)==2 && !any(sapply(names_list, is.null))){
    return(unlist(names_list))
  }else{
    return(paste0("s-", 1:n))
  }
}

# x is a int number
# return: all factors of x
all_factors <- function(x){
  x <- as.integer(x)
  div <- seq_len(abs(x))
  factors <- div[x %% div == 0L]
  return(factors)
}

# Solve a System of Equations: Robust function
lin_sys <- function(msx, mdx){
  if(NCOL(msx)>100){
    if(!is(msx, "symmetricMatrix")){
      #msx <- methods::as(methods::as(forceSymmetric(msx), "unpackedMatrix"), "symmetricMatrix")
      msx <- forceSymmetric(msx)
      mdx <- methods::as(mdx, "CsparseMatrix")
    }
  }

  out <- tryCatch(solve(msx, mdx), error = function(cond){
    tryCatch(solve(qr(msx), mdx), error = function(cond){
      backsolve(chol(msx), mdx)
    })
  })
  out[is.na(out)] <- 0
  return(out)
}

# Fast cov2cor
covcor <- function(V){
  p <- (d <- dim(V))[1L]
  if(length(d) != 2L || p != d[2L])
    cli_abort("{.arg V} is not a square numeric matrix", call = NULL)
  Is <- sqrt(1/diag(V))
  if(any(!is.finite(Is)))
    cli_warn(c("!"="diag(.) had 0 or NA entries; non-finite result is doubtful"), call = NULL)
  r <- V / tcrossprod(diag(V) ^ 0.5)
  r[cbind(1L:p, 1L:p)] <- 1
  r
}

# Sample covariance matrix
sample_estim <- function(x, mse = TRUE){
  if(mse){
    if(any(is.na(x))){
      x <- remove_na(x)
    }
    if(any(is.na(x))){
      if(is.vector(x)){
        n <- sum(!is.na(x))
      }else{
        n <- colSums(!is.na(x))
      }
      n <- tcrossprod(n, rep(1, length(n)))
      x[is.na(x)] <- 0
      crossprod(x) / pmin(t(n), n)
    }else{
      crossprod(x) / NROW(x)
    }
  }else{
    stats::var(x, na.rm = TRUE, use = "na.or.complete")
  }
}

# find the bottom time series given strc_mat
find_bts <- function(strc_mat){
  strc_mat <- Matrix(strc_mat, sparse = TRUE)
  strc_mat@i[strc_mat@p[-1]] + 1
}

# Sparse matrix to dense
sparse2dense <- function(input, sparse = TRUE){
  if(!sparse){
    class_check <- "Matrix"
    class_out <- "matrix"
  }else{
    class_check <- "matrix"
    class_out <- "Matrix"
  }
    if(is.list(input)){
      output <- lapply(input, function(x){
        if(is(x, class_check)){
          as(x, class_out)
        }else{
          x
        }
      })
      return(output)
    }else{
      if(is(input, class_check)){
        return(as(input, class_out))
      }else{
        return(input)
      }
    }
}

# Re-arrange the strc_mat according to the new bts
transform_strc_mat <- function(strc_mat, bts){
  if (length(bts) != NCOL(strc_mat)){
    stop(simpleError(sprintf('length of basis set should be %d', NCOL(strc_mat))))
  }
  S1 <- strc_mat[bts,]
  S2 <- strc_mat[-bts,]
  transitionMat <- solve(S1, Diagonal(NCOL(strc_mat)))
  strc_mat[-bts,] <- S2 %*% transitionMat
  strc_mat[bts,] <- Diagonal(NCOL(strc_mat))
  return(strc_mat)
}

# Remove NA values (row and columns)
remove_na <- function(x){
  inax <- is.na(x)
  if(any(inax)){
    out <- stats::na.omit(x)
    if(NROW(out) == 0){
      x <- x[, !(colSums(!inax) == 0)]
      inax <- is.na(x)
      #out <- stats::na.omit(x)
    }
    row_na <- rowSums(inax)
    x <- x[!(rowSums(inax) == NCOL(x)), , drop = FALSE]
  }
  return(x)
}



