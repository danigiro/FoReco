reco <- function(approach, base, immutable = NULL, nn = NULL, ...){
  # Fri Feb  9 2024
  tsp(base) <- NULL # Remove ts

  if(any(approach %in% c("proj_osqp", "strc_osqp",
                         "proj_immutable",
                         "proj_immutable2", "strc_immutable") | is.null(immutable))){
    class_base <- approach
  }else{
    class_base <- paste0(approach, "_immutable")
  }

  # Set class of 'base' to include 'approach' and reconcile
  class(approach) <- c(class(approach), class_base)
  rmat <- .reco(approach = approach, base = base, nn = nn, immutable = immutable, ...)

  # Check if 'nn' is provided and adjust 'rmat' accordingly
  if(!is.null(nn)){
    if(nn == "osqp"){
      nn <- paste(approach, nn, sep = "_")
    }

    if(!all(rmat >= -sqrt(.Machine$double.eps))){
      class(approach)[length(class(approach))] <- nn
      rmat <- .reco(approach = approach, base = base, nn = nn, reco = rmat,
                    immutable = immutable, ...)
    } else if(!all(rmat >= 0)){
      rmat[rmat < 0] <- 0
    }
  }
  return(rmat)
}

.reco <- function(approach, ...){
  UseMethod("reco", approach)
}

reco.default <- function(approach, ...){
  cli_abort(c("Please provide a valid approach.",
    "i"= "{.strong Optimal}: {.code proj}, {.code strc}, {.code proj_osqp}, {.code strc_osqp}",
    "i"= "{.strong Non-negative}: {.code sntz}, {.code proj_osqp}, {.code strc_osqp}",
    "i"= "{.strong Immutable}: {.code proj}, {.code strc}, {.code proj_osqp}, {.code strc_osqp}"))
}

reco.proj <- function(base, cons_mat, cov_mat, ...){

  # check input
  if(missing(base) | missing(cons_mat) | missing(cov_mat)){
    cli_abort("Mandatory arguments: {.arg base}, {.arg cons_mat} and {.arg cov_mat}.",
              call = NULL)
  }

  if(NCOL(cons_mat) != NROW(cov_mat) | NCOL(base) != NROW(cov_mat)){
    cli_abort("The size of the matrices does not match.", call = NULL)
  }

  # Point reconciled forecasts
  lm_dx <- methods::as(Matrix::tcrossprod(cons_mat, base), "CsparseMatrix")
  lm_sx <- methods::as(Matrix::tcrossprod(cons_mat %*% cov_mat, cons_mat), "CsparseMatrix")
  reco <- base - t(cov_mat %*% Matrix::crossprod(cons_mat, lin_sys(lm_sx, lm_dx)))
  return(as.matrix(reco))
}

reco.strc <- function(base, strc_mat, cov_mat, ...){
  # check input
  if(missing(base) | missing(strc_mat) | missing(cov_mat)){
    cli_abort("Mandatory arguments: {.arg base}, {.arg strc_mat} and {.arg cov_mat}.",
              call = NULL)
  }

  if(is.null(strc_mat)){
    cli_abort("Please provide a valid {.arg agg_mat} for the structural approach.",
              call = NULL)
  }

  if(NROW(strc_mat) != NROW(cov_mat) | NCOL(base) != NROW(cov_mat)){
    cli_abort("The size of the matrices does not match.", call = NULL)
  }

  # Point reconciled forecasts
  if(isDiagonal(cov_mat)){
    cov_mat_inv <- .sparseDiagonal(x = diag(cov_mat)^(-1))
    StWm <- Matrix::crossprod(strc_mat, cov_mat_inv)
    lm_sx1 <- methods::as(StWm %*% strc_mat, "CsparseMatrix")
    lm_dx1 <- methods::as(Matrix::tcrossprod(StWm, base), "CsparseMatrix")
    reco <- t(strc_mat %*% lin_sys(lm_sx1, lm_dx1))
    return(as.matrix(reco))
  } else {
    Q <- lin_sys(cov_mat, strc_mat)
    lm_sx1 <- methods::as(t(strc_mat) %*% Q, "CsparseMatrix")
    lm_dx1 <- methods::as(t(base %*% Q), "CsparseMatrix")
    reco <- t(strc_mat %*% lin_sys(lm_sx1, lm_dx1))
    return(as.matrix(reco))
  }
}

reco.proj_osqp <- function(base, cons_mat, cov_mat,
                           nn = NULL, id_nn = NULL, bounds = NULL,
                           reco = NULL, settings = NULL, immutable = NULL, ...){
  # check input
  if(missing(base) | missing(cons_mat) | missing(cov_mat)){
    cli_abort("Mandatory arguments: {.arg base}, {.arg cons_mat} and {.arg cov_mat}.",
              call = NULL)
  }

  if(NCOL(cons_mat) != NROW(cov_mat) | NCOL(base) != NROW(cov_mat)){
    cli_abort("The size of the matrices does not match.", call = NULL)
  }

  if(is.null(id_nn)){
    id_nn <- rep(1, NCOL(base))
  }

  if(!is.null(nn) & !is.null(reco)){
    id <- which(rowSums(reco < (-sqrt(.Machine$double.eps))) != 0)
    if(!is.null(bounds)){
      id_b <- which(apply(reco, 1, function(x) all(bounds[,1] <= x) & all(bounds[,2] >= x)))
      if(length(id_b) > 0){
        id <- sort(unique(c(id, id_b)))
      }
    }

    if(length(id) == 0){
      reco[reco < 0] <- 0
      return(reco)
    }
  } else {
    id <- 1:NROW(base)
    reco <- base
  }

  c <- ncol(cons_mat)
  r <- nrow(cons_mat)
  # Linear constrains H = 0
  l <- rep(0, r)
  u <- rep(0, r)
  A <- cons_mat

  # P matrix
  if(isDiagonal(cov_mat)){
    P <- Diagonal(x = diag(cov_mat)^(-1))
    #P <- .sparseDiagonal(x = diag(cov_mat)^(-1))
  } else {
    R <- chol(cov_mat)
    P <- as(chol2inv(R), "CsparseMatrix")
  }

  # nn constraints (only on the building block variables)
  if(!is.null(nn)){
    if(!(nn %in% c("osqp", TRUE, "proj_osqp"))){
      cli_warn("Non-negative reconciled forecasts obtained with osqp.", call = NULL)
    }
    A <- rbind(A, .sparseDiagonal(c)[id_nn == 1, ])
    l <- c(l, rep(0, sum(id_nn)))
    u <- c(u, rep(Inf, sum(id_nn)))
  }

  # other constraints
  if(!is.null(bounds)){
    bounds_rows <- rowSums(abs(bounds) == Inf) < 2
    A <- rbind(A, Diagonal(c)[bounds_rows, ])
    l <- c(l, bounds[bounds_rows, 1, drop = TRUE])
    u <- c(u, bounds[bounds_rows, 2, drop = TRUE])
  }

  if(is.null(settings)){
    settings <- osqpSettings(
      verbose = FALSE,
      eps_abs = 1e-5,
      eps_rel = 1e-5,
      polish_refine_iter = 100,
      polish = TRUE
    )
  }

  # OSQP
  osqp_step <- apply(base[id, , drop = FALSE], 1, function(x){
    if(!is.null(immutable)){
      A <- rbind(A, .sparseDiagonal(c)[immutable, , drop = FALSE])
      l <- c(l, x[immutable])
      u <- c(u, x[immutable])
    }
    q <- (-1) * t(P) %*% as.vector(x)
    rec <- solve_osqp(P, q, A, l, u, settings)

    # Fix a problem of osqp
    if(rec$info$status_val == -4){
      u[u == Inf] <- max(x)*100
      rec <- solve_osqp(P, q, A, l, u, settings)
    }

    out <- list()
    out$reco <- rec$x

    if(rec$info$status_val != 1){
      cli_warn(c("x"="OSQP failed: check the results.",
                 "i"="OSQP flag = {rec$info$status_val}",
                 "i"="OSQP pri_res = {rec$info$pri_res}"), call = NULL)
    }

    out$info <- c(rec$info$obj_val, rec$info$run_time, rec$info$iter,
                  rec$info$pri_res, rec$info$status_val, rec$info$status_polish)

    return(out)
  })
  osqp_step <- do.call("rbind", osqp_step)

  # Point reconciled forecasts
  reco[id, ] <- do.call("rbind", osqp_step[, "reco"])
  if(!is.null(nn)){
    reco[which(reco <= sqrt(.Machine$double.eps))] <- 0
  }

  class(reco) <- setdiff(class(reco), "proj_osqp")

  info <- do.call("rbind", osqp_step[, "info"])
  colnames(info) <- c(
    "obj_val", "run_time", "iter", "pri_res",
    "status", "status_polish"
  )
  rownames(info) <- id
  attr(reco, "info") <- info
  return(reco)
}

reco.strc_osqp <- function(base, strc_mat, cov_mat,
                           nn = NULL, id_nn = NULL, bounds = NULL,
                           reco = NULL, settings = NULL, immutable = NULL, ...){
  # check input
  if(missing(base) | missing(strc_mat) | missing(cov_mat)){
    cli_abort("Mandatory arguments: {.arg base}, {.arg strc_mat} and {.arg cov_mat}.",
              call = NULL)
  }

  if(is.null(strc_mat)){
    cli_abort("Please provide a valid {.arg agg_mat} for the structural approach.",
              call = NULL)
  }

  if(NROW(strc_mat) != NROW(cov_mat) | NCOL(base) != NROW(cov_mat)){
    cli_abort("The size of the matrices does not match.", call = NULL)
  }

  if(is.null(id_nn)){
    bts <- find_bts(strc_mat)
    id_nn <- rep(0, NCOL(base))
    id_nn[bts] <- 1
  }

  if(!is.null(nn) & !is.null(reco)){
    id <- which(rowSums(reco < (-sqrt(.Machine$double.eps))) != 0)
    if(!is.null(bounds)){
      id_b <- which(apply(reco, 1, function(x) all(bounds[,1] <= x) & all(bounds[,2] >= x)))
      if(length(id_b) > 0){
        id <- sort(unique(c(id, id_b)))
      }
    }
    if(length(id) == 0){
      reco[reco < 0] <- 0
      return(reco)
    }
  } else {
    id <- 1:NROW(base)
    reco <- base
  }

  r <- NROW(strc_mat)
  c <- NCOL(strc_mat)
  A <- NULL
  l <- NULL
  u <- NULL

  # P matrix and q1 vector
  if(isDiagonal(cov_mat)){
    Q <- Diagonal(x = diag(cov_mat)^(-1))
    P <- t(strc_mat) %*% Q %*% strc_mat
    q1 <- (-1) * t(Q %*% strc_mat)
  } else {
    Q <- lin_sys(cov_mat, strc_mat)
    P <- t(strc_mat) %*% Q
    q1 <- (-1) * t(Q)
  }

  # nn constraints (only on the building block variables - bottom variables)
  if(!is.null(nn)){
    if(!(nn %in% c("osqp", TRUE, "strc_osqp"))){
      cli_warn("Non-negative reconciled forecasts obtained with osqp.", call = NULL)
    }
    A <- .sparseDiagonal(c)
    l <- rep(0, sum(c))
    u <- rep(Inf, sum(c))
  }

  # other constraints
  if(!is.null(bounds)){
    bounds_rows <- rowSums(abs(bounds) == Inf) < 2
    A <- rbind(A, strc_mat[bounds_rows, ,drop = FALSE])
    l <- c(l, bounds[bounds_rows, 1, drop = TRUE])
    u <- c(u, bounds[bounds_rows, 2, drop = TRUE])
  }

  if(is.null(settings)){
    settings <- osqpSettings(
      verbose = FALSE,
      eps_abs = 1e-5,
      eps_rel = 1e-5,
      polish_refine_iter = 100,
      polish = TRUE
    )
  }

  # OSQP
  osqp_step <- apply(base[id, , drop = FALSE], 1, function(x){
    if(!is.null(immutable)){
      A <- rbind(A, strc_mat[immutable, , drop = FALSE])
      l <- c(l, x[immutable])
      u <- c(u, x[immutable])
    }
    q <- q1 %*% as.vector(x)
    rec <- solve_osqp(P, q, A, l, u, settings)

    # Fix a problem of osqp
    if(rec$info$status_val == -4){
      u[u == Inf] <- max(x)*100
      rec <- solve_osqp(P, q, A, l, u, settings)
    }

    out <- list()
    out$reco <- as.numeric(strc_mat %*% rec$x)

    if(rec$info$status_val != 1){
      cli_warn(c("x"="OSQP failed: check the results.",
                 "i"="OSQP flag = {rec$info$status_val}",
                 "i"="OSQP pri_res = {rec$info$pri_res}"), call = NULL)
    }

    out$info <- c(
      rec$info$obj_val, rec$info$run_time, rec$info$iter, rec$info$pri_res,
      rec$info$status_val, rec$info$status_polish
    )

    return(out)
  })
  osqp_step <- do.call("rbind", osqp_step)

  # Point reconciled forecasts
  reco[id, ] <- do.call("rbind", osqp_step[, "reco"])
  if(!is.null(nn)){
    reco[which(reco <= sqrt(.Machine$double.eps))] <- 0
  }

  class(reco) <- setdiff(class(reco), "strc_osqp")

  info <- do.call("rbind", osqp_step[, "info"])
  colnames(info) <- c(
    "obj_val", "run_time", "iter", "pri_res",
    "status", "status_polish"
  )
  rownames(info) <- id
  attr(reco, "info") <- info
  return(reco)
}


reco.sntz <- function(base, reco, strc_mat, cov_mat, id_nn = NULL, settings = NULL, ...){
  # Check input
  if(missing(strc_mat) | missing(cov_mat)){
    cli_abort("Mandatory arguments: {.arg strc_mat} and {.arg cov_mat}.",
              call = NULL)
  }

  if(missing(reco)){
    reco <- base
  }

  if(is.null(strc_mat)){
    cli_abort(c("Argument {.arg agg_mat} is missing. The {.strong sntz} approach
                is available only for hierarchical/groupped time series."), call = NULL)
  }

  if(is.null(id_nn)){
    bts <- find_bts(strc_mat)
    id_nn <- rep(0, NCOL(reco))
    id_nn[bts] <- 1
  }

  bts <- reco[, id_nn == 1, drop = FALSE]
  if(is.null(settings$type)){
    sntz_type <- "bu"
  }else{
    sntz_type <- settings$type
  }
  tol <- sqrt(.Machine$double.eps)
  switch(sntz_type,
         bu = {
           bts[bts<tol] <- 0
         },
         tdp = {
           bts <- t(apply(bts, 1, function(x){
             d <- sum(x[x<tol])
             Ip <- (x>tol)
             w <- rep(0, length(x))
             w[Ip] <- x[Ip]/sum(x[Ip])
             while(any(w[Ip] > x[Ip]/abs(d))){
               Ip[Ip][w[Ip] > x[Ip]/abs(d)] <- FALSE
               d <- sum(x[!Ip])
               w <- rep(0, length(x))
               w[Ip] <- x[Ip]/sum(x[Ip])
             }
             x[!Ip] <- 0
             x + w*d
           }, simplify = TRUE))
         },
         tdsp = {
           bts <- t(apply(bts, 1, function(x){
             d <- sum(x[x<tol])
             Ip <- (x>tol)
             w <- rep(0, length(x))
             w[Ip] <- (x[Ip]^2)/sum(x[Ip]^2)
             while(any(w[Ip] > x[Ip]/abs(d))){
               Ip[Ip][w[Ip] > x[Ip]/abs(d)] <- FALSE
               d <- sum(x[!Ip])
               w <- rep(0, length(x))
               w[Ip] <- (x[Ip]^2)/sum(x[Ip]^2)
             }
             x[!Ip] <- 0
             x + w*d
           }, simplify = TRUE))
         },
         tdvw = {
           sigma2 <- diag(cov_mat)[id_nn == 1]
           bts <- t(apply(bts, 1, function(x){
             Ip <- (x>tol)
             d <- sum(x[!Ip])
             w <- rep(0, length(x))
             w[Ip] <- sigma2[Ip]/sum(sigma2[Ip])
             while(any(w[Ip] > x[Ip]/abs(d))){
               Ip[Ip][w[Ip] > x[Ip]/abs(d)] <- FALSE
               d <- sum(x[!Ip])
               w <- rep(0, length(x))
               w[Ip] <- sigma2[Ip]/sum(sigma2[Ip])
             }
             x[!Ip] <- 0
             x + w*d
           }, simplify = TRUE))
         })

  as.matrix(bts %*% t(strc_mat))
}

reco.proj_immutable <- function(base, cons_mat, cov_mat, immutable = NULL, ...){

  # Check input
  if(missing(base) | missing(cons_mat) | missing(cov_mat)){
    cli_abort("Mandatory arguments: {.arg base}, {.arg cons_mat} and {.arg cov_mat}.",
              call = NULL)
  }

  if(NCOL(cons_mat) != NROW(cov_mat) | NCOL(base) != NROW(cov_mat)){
    cli_abort("The size of the matrices does not match.", call = NULL)
  }

  if(is.null(immutable)){
    cli_warn("No immutable forecasts!")
    reco <- reco.proj(base = base, cons_mat = cons_mat, cov_mat = cov_mat)
    return(reco)
  } else if(max(immutable) > NCOL(base)){
    cli_abort("{.code max(immutable)} must be less or equal to {NCOL(base)}", call = NULL)
  }

  # Complete constraints matrix (immutable forecasts + linear constraints)
  imm_cons_mat <- .sparseDiagonal(NCOL(base))[immutable, , drop = FALSE]
  imm_cons_vec <- base[, immutable, drop = FALSE]
  compl_cons_mat <- rbind(cons_mat, imm_cons_mat)
  compl_cons_vec <- cbind(
    Matrix(0, nrow = NROW(imm_cons_vec), ncol = NROW(cons_mat)),
    imm_cons_vec
  )

  # check immutable feasibility
  # TODO: can proj_immutable2 be more stable than proj_immutable?
  # Answer issue: https://github.com/danigiro/FoReco/issues/6#issue-2397642027 (@AngelPone)
  if(rankMatrix(cons_mat) + length(immutable) != rankMatrix(compl_cons_mat)){
    cli_abort("There is no solution with this {.arg immutable} set.",  call = NULL)
  }

  # Point reconciled forecasts
  lm_dx <- t(compl_cons_vec) - Matrix::tcrossprod(compl_cons_mat, base)
  lm_sx <- methods::as(Matrix::tcrossprod(compl_cons_mat %*% cov_mat,
                                          compl_cons_mat), "CsparseMatrix")
  reco <- base + t(cov_mat %*% Matrix::crossprod(compl_cons_mat, lin_sys(lm_sx, lm_dx)))
  return(as.matrix(reco))
}

reco.proj_immutable2 <- function(base, cons_mat, cov_mat, immutable, ...){
  # Check input
  if(missing(base) | missing(cons_mat) | missing(cov_mat)){
    cli_abort("Mandatory arguments: {.arg base}, {.arg cons_mat} and {.arg cov_mat}.",
              call = NULL)
  }

  if(NCOL(cons_mat) != NROW(cov_mat) | NCOL(base) != NROW(cov_mat)){
    cli_abort("The size of the matrices does not match.", call = NULL)
  }

  if(is.null(immutable)){
    cli_warn("No immutable forecasts!", call = NULL)
    reco <- reco.proj(base = base, cons_mat = cons_mat, cov_mat = cov_mat)
    return(reco)
  } else if(max(immutable) > NCOL(base)){
    cli_abort("{.code max(immutable)} must be less or equal to {NCOL(base)}", call = NULL)
  }

  cons_mat_red <- cons_mat[ , -immutable, drop = FALSE]
  cons_vec <- apply(-cons_mat[ , immutable, drop = FALSE], 1, function(w)
    rowSums(base[, immutable, drop = FALSE]%*%w))
  cov_mat_red <- cov_mat[-immutable , -immutable, drop = FALSE]
  base_red <- base[, -immutable, drop = FALSE]

  if(length(immutable)>2){
    check <- which(rowSums(cons_mat_red != 0) == 0)
    if(length(check) > 0 && any(cons_vec[, check] > sqrt(.Machine$double.eps))){
      cli_abort("There is no solution with this {.arg immutable} set.",  call = NULL)
    }
  }

  # Point reconciled forecasts
  lm_dx <- t(cons_vec) - Matrix::tcrossprod(cons_mat_red, base_red)
  lm_sx <- methods::as(Matrix::tcrossprod(cons_mat_red %*% cov_mat_red,
                                          cons_mat_red), "CsparseMatrix")
  reco <- Matrix(base)
  reco[ , -immutable] <- (base_red + t(cov_mat_red %*% Matrix::crossprod(cons_mat_red,
                                                                         lin_sys(lm_sx, lm_dx))))
  if(any(is.nan(reco))){
    cli_abort("There is no solution with this {.arg immutable} set.",  call = NULL)
  }
  return(as.matrix(reco))
}

reco.strc_immutable <- function(base, strc_mat, cov_mat, immutable = NULL, ...){
  # Check input
  if(missing(base) | missing(strc_mat) | missing(cov_mat)){
    cli_abort("Mandatory arguments: {.arg base}, {.arg strc_mat} and {.arg cov_mat}.",
              call = NULL)
  }

  if(is.null(strc_mat)){
    cli_abort("Please provide a valid {.arg agg_mat} for the structural approach.",
              call = NULL)
  }

  if(NROW(strc_mat) != NROW(cov_mat) | NCOL(base) != NROW(cov_mat)){
    cli_abort("The size of the matrices does not match.", call = NULL)
  }

  if(is.null(immutable)){
    cli_warn("No immutable forecasts!", call = NULL)
    reco <- reco.strc(base = base, strc_mat = strc_mat, cov_mat = cov_mat)
    return(reco)
  } else if(max(immutable) > NCOL(base)){
    cli_abort("{.code max(immutable)} must be less or equal to {NCOL(base)}", call = NULL)
  }

  # Code idea: https://github.com/AngelPone/chf
  bts <- find_bts(strc_mat)
  immutable <- sort(immutable)
  candidate <- setdiff(bts, immutable)
  determined <- setdiff(1:NROW(strc_mat), bts)
  mutable <- candidate
  if(any(immutable %in% determined)){
    i <- max(which(immutable %in% determined))
    while(i > 0){
      corr_leaves <- bts[which(strc_mat[immutable[i], ] != 0)]
      free_leaves <- setdiff(corr_leaves, c(immutable, determined))
      if(length(free_leaves) == 0){
        if(all(corr_leaves %in% immutable)){
          cli_warn("All children of {immutable[i]}th series are immutable, it is removed from the condition.", call = NULL)
          immutable <- immutable[immutable != immutable[i]]
          i <- i - 1
          next
        }else{
          cli_abort("There is no solution with this {.arg immutable} set.",  call = NULL)
        }
      }
      determined <- determined[determined != immutable[i]]
      determined <- c(determined, free_leaves[1])
      mutable <- mutable[mutable != free_leaves[1]]
      i <- i - 1
    }
  }
  new_basis <- sort(c(sort(mutable), immutable))
  snew <- transform_strc_mat(strc_mat, new_basis)
  S1 <- snew[-immutable, -which(new_basis %in% immutable), drop = FALSE]
  S2 <- snew[-new_basis, which(new_basis %in% immutable), drop = FALSE]
  S2u <- base[, immutable, drop = FALSE]%*%t(S2)
  base2 <- Matrix(base)
  base2[, -new_basis] <- (base[, -new_basis, drop = FALSE] - S2u)
  reco_bts <- base2[, new_basis, drop = FALSE]
  cov_mat_red <- cov_mat[-immutable, -immutable, drop = FALSE]
  tmp <- reco.strc(base2[,-immutable, drop = FALSE], strc_mat = S1,
                   cov_mat = cov_mat_red)[, find_bts(S1), drop = FALSE]
  reco_bts[, !(new_basis %in% immutable)] <- tmp
  reco <- reco_bts %*% t(snew)
  return(as.matrix(reco))
}

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

reco.kann <- function(base, cons_mat, cov_mat, nn = NULL,
                      reco = NULL, settings = NULL, immutable = NULL, ...){
  # Check input
  if(missing(base) | missing(cons_mat) | missing(cov_mat)){
    cli_abort("Mandatory arguments: {.arg base}, {.arg cons_mat} and {.arg cov_mat}.",
              call = NULL)
  }

  tol <- settings$tol
  if(is.null(tol)){
    tol <- sqrt(.Machine$double.eps)
  }

  itmax <- settings$itmax
  if(is.null(itmax)){
    itmax <- 100
  }

  if(is.null(reco)){
    if(is.null(immutable)){
      reco <- reco.proj(base = base, cons_mat = cons_mat, cov_mat = cov_mat)
    }else{
      reco <- reco.proj_immutable(base = base, cons_mat = cons_mat,
                                  cov_mat = cov_mat, immutable = immutable)
    }
  }

  if(is.null(nn)){
    return(reco)
  }

  if(all(reco>tol)){
    reco[reco<0] <- 0
    return(reco)
  }

  rowid <- which(rowSums(reco < (-sqrt(.Machine$double.eps))) != 0)
  kann_step <- apply(reco[rowid, , drop = FALSE], 1, function(x){
    start <- Sys.time()
    for(i in 1:itmax){
      x[x < sqrt(.Machine$double.eps)] <- 0
      if(is.null(immutable)){
        x <- reco.proj(base = rbind(x), cons_mat = cons_mat, cov_mat = cov_mat)
      }else{
        x <- reco.proj_immutable(base = rbind(x), cons_mat = cons_mat, cov_mat = cov_mat,
                                 immutable = immutable)
      }
      x <- as.numeric(x)

      if(all(x >= (-tol))){
        flag <- 1
        break
      }else{
        flag <- -2
      }
    }

    if(flag == 1){
      x[x <= sqrt(.Machine$double.eps)] <- 0
      if(i == itmax){
        flag <- 2
      }
    }
    end <- Sys.time()
    out <- list()
    out$reco <- x
    out$info <- c(difftime(end,start, units = "secs"), i, flag)
    if(flag %in% c(-2, 2)){
      cli_warn(c("x"="KANN failed: check the results.",
                 "i"="Flag = {flag},  tol = {tol}, itmax = {itmax}"), call = NULL)
    }
    out
  })

  kann_step <- do.call("rbind", kann_step)

  # Point reconciled forecasts
  reco[rowid, ] <- do.call("rbind", kann_step[, "reco"])

  info <- do.call("rbind", kann_step[, "info"])
  colnames(info) <- c("run_time", "iter", "status")
  rownames(info) <- rowid
  attr(reco, "info") <- info
  return(reco)
}

reco.gauss <- function(...){
  # TODO
  return(NULL)
}

reco.fbpp <- function(base, cons_mat, cov_mat, id_nn = NULL, nn = NULL,
                       reco = NULL, settings = NULL, immutable = NULL, ...){
  # Check input
  if(missing(base) | missing(cons_mat) | missing(cov_mat)){
    cli_abort("Mandatory arguments: {.arg base}, {.arg cons_mat} and {.arg cov_mat}.",
              call = NULL)
  }

  tol <- settings$tol
  if(is.null(tol)){
    tol <- sqrt(.Machine$double.eps)
  }

  itmax <- settings$itmax
  if(is.null(itmax)){
    itmax <- 100
  }

  if(is.null(reco)){
    if(is.null(immutable)){
      reco <- reco.proj(base = base, cons_mat = cons_mat, cov_mat = cov_mat)
    }else{
      reco <- reco.proj_immutable(base = base, cons_mat = cons_mat, cov_mat = cov_mat,
                                  immutable = immutable)
    }
  }

  if(is.null(nn)){
    return(reco)
  }

  if(all(reco>-tol)){
    reco[reco<=sqrt(.Machine$double.eps)] <- 0
    return(reco)
  }

  if(is.null(id_nn)){
    qrtmp <- base::qr(cons_mat)
    id_nn <- rep(1, NCOL(base))
    id_nn[qrtmp$pivot[1:qrtmp$rank]] <- 0
  }

  rowid <- which(rowSums(reco < (-sqrt(.Machine$double.eps))) != 0)
  fbpp_step <- apply(reco[rowid, , drop = FALSE], 1, function(x){
    start <- Sys.time()
    idx <- NULL
    for(i in 1:itmax){
      idx <- c(idx, which(x < -tol))
      idx <- idx[idx %in% which(id_nn == 1)]
      block <- sparseMatrix(i = 1:length(idx),
                            j = idx,
                            x = 1,
                            dims = c(length(idx), NCOL(cons_mat)))
      cons_matx <- rbind(cons_mat, block)
      if(is.null(immutable)){
        x <- reco.proj(base = rbind(x), cons_mat = cons_matx, cov_mat = cov_mat)
      }else{
        x <- reco.proj_immutable(base = rbind(x), cons_mat = cons_matx, cov_mat = cov_mat,
                                 immutable = immutable)
      }
      x <- as.numeric(x)

      if(all(x >= (-tol))){
        flag <- 1
        break
      }else{
        flag <- -2
      }
    }

    if(flag == 1){
      x[x <= sqrt(.Machine$double.eps)] <- 0
      if(i == itmax){
        flag <- 2
      }
    }
    end <- Sys.time()
    out <- list()
    out$reco <- x
    out$info <- c(difftime(end,start, units = "secs"), i, flag)
    out$idx <- idx
    if(flag %in% c(-2, 2)){
      cli_warn(c("x"="FBPP failed: check the results.",
                 "i"="Flag = {flag},  tol = {tol}, itmax = {itmax}"), call = NULL)
    }
    out
  })

  fbpp_step <- do.call("rbind", fbpp_step)

  # Point reconciled forecasts
  reco[rowid, ] <- do.call("rbind", fbpp_step[, "reco"])

  info <- do.call("rbind", fbpp_step[, "info"])
  colnames(info) <- c("run_time", "iter", "status")
  rownames(info) <- rowid
  attr(reco, "info") <- info
  return(reco)
}
