reco <- function(approach, base, immutable = NULL, nn = NULL, bounds = NULL, ...){
  # Fri Feb  9 2024
  tsp(base) <- NULL # Remove ts

  if(approach %in% c("proj_osqp", "strc_osqp", "proj_immutable",
                     "proj_immutable2", "strc_immutable") | is.null(immutable)){
    class_base <- approach
  }else{
    class_base <- paste0(approach, "_immutable")
  }

  # Set class of 'base' to include 'approach' and reconcile
  class(approach) <- c(class(approach), class_base)
  rmat <- .reco(approach = approach, base = base, nn = nn,
                immutable = immutable, bounds = bounds, ...)

  # Check if 'nn' is provided and adjust 'rmat' accordingly
  if(!is.null(nn)){
    if(nn %in% c("osqp", TRUE)){
      nn <- paste(approach, "osqp", sep = "_")
    }

    if(!all(rmat >= -sqrt(.Machine$double.eps))){
      class(approach)[length(class(approach))] <- nn
      rmat <- .reco(approach = approach, base = base, nn = nn, reco = rmat,
                    immutable = immutable, ...)
    }else if(!all(rmat >= 0)){
      rmat[rmat < 0] <- 0
    }
  }

  if(!is.null(bounds)){
    nbid <- bounds[,1,drop = TRUE]

    checkb <- apply(rmat, 1, function(x){
      idl <- any(x[nbid]<bounds[,2,drop = TRUE] - sqrt(.Machine$double.eps))
      idb <- any(x[nbid]>bounds[, 3, drop = TRUE] + sqrt(.Machine$double.eps))

      idl0 <- any(x[nbid]<bounds[,2,drop = TRUE])
      idb0 <- any(x[nbid]>bounds[, 3, drop = TRUE])
      c(any(c(idl, idb)), any(c(idl0, idb0)))
    })

    if(any(checkb[1,])){
      if(is.null(attr(bounds, "approach")) || attr(bounds, "approach") == "osqp"){
        attr(bounds, "approach") <- paste(approach, "osqp", sep = "_")
      }

      class(approach)[length(class(approach))] <- attr(bounds, "approach")
      rmat <- .reco(approach = approach, base = base, bounds = bounds, reco = rmat,
                    immutable = immutable, ...)
    }else if(any(checkb[2,])){
      rmat <- t(apply(rmat, 1, function(x){
        id <- x[nbid]<=bounds[,2,drop = TRUE]+sqrt(.Machine$double.eps)
        x[nbid][id] <- bounds[,2,drop = TRUE][id]

        id <- x[nbid]>=bounds[, 3, drop = TRUE]-sqrt(.Machine$double.eps)
        x[nbid][id] <- bounds[, 3, drop = TRUE][id]
        x
      }))
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
      id_b <- which(apply(reco[, bounds[,1], drop = FALSE], 1,
                          function(x) any(x <= bounds[,2]) | any(x >= bounds[,3])))
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
    A <- rbind(A, Diagonal(c)[bounds[,1,drop = TRUE], ])
    l <- c(l, bounds[,2,drop = TRUE])
    u <- c(u, bounds[,3,drop = TRUE])
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

    if(!is.null(bounds)){
      nbid <- bounds[,1,drop = TRUE]
      id <- out$reco[nbid]<=bounds[,2,drop = TRUE]+sqrt(.Machine$double.eps)
      out$reco[nbid][id] <- bounds[,2,drop = TRUE][id]

      id <- out$reco[nbid]>=bounds[, 3, drop = TRUE]-sqrt(.Machine$double.eps)
      out$reco[nbid][id] <- bounds[, 3, drop = TRUE][id]
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
    A <- rbind(A, strc_mat[bounds[,1,drop = TRUE], ,drop = FALSE])
    l <- c(l, bounds[,2,drop = TRUE])
    u <- c(u, bounds[,3,drop = TRUE])
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

    if(!is.null(bounds)){
      nbid <- bounds[,1,drop = TRUE]
      id <- out$reco[nbid]<=bounds[,2,drop = TRUE]+sqrt(.Machine$double.eps)
      out$reco[nbid][id] <- bounds[,2,drop = TRUE][id]

      id <- out$reco[nbid]>=bounds[, 3, drop = TRUE]-sqrt(.Machine$double.eps)
      out$reco[nbid][id] <- bounds[, 3, drop = TRUE][id]
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
  rank_cm <- rankMatrix(cons_mat, method = "qr", warn.t = FALSE)
  rank_ccm <-  rankMatrix(compl_cons_mat, method = "qr", warn.t = FALSE)
  if(rank_cm + length(immutable) != rank_ccm){
    cli_warn(c(
      "There may be redundant constraints.",
      "i" = "A solution may exist, but numerical issues could occur.",
      "x" = "Alternatively, there is no solution with this {.arg immutable} set.",
      "i" = "Please check the {.arg immutable} parameter and carefully inspect the resulting solution."
    ), call = NULL)
    #cli_abort("There is no solution with this {.arg immutable} set.",  call = NULL)
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
  if(is.vector(cons_vec)){
    cons_vec <- unname(rbind(cons_vec))
  }
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
          cli_warn("All children of {immutable[i]}th series are immutable, it is removed.",
                   call = NULL)
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

reco.bpv <- function(base, strc_mat, cov_mat, cons_mat, id_nn = NULL, nn = NULL,
                     reco = NULL, settings = NULL, immutable = NULL, approach = "proj", ...){
  if(!is.null(immutable)){
    cli_warn("immutable not supported")
  }

  if(is.null(strc_mat)){
    cli_abort("Please provide a valid {.arg agg_mat} for bpv.",
              call = NULL)
  }

  if(is.null(reco)){
    if(approach == "proj"){
      reco <- reco.proj(base = base, cons_mat = cons_mat, cov_mat = cov_mat)
    }else{
      reco <- reco.strc(base = base, strc_mat = strc_mat, cov_mat = cov_mat)
    }
  }

  if(is.null(nn) | all(reco>=0)){
    return(reco)
  }

  if(is.null(id_nn)){
    bts <- find_bts(strc_mat)
    id_nn <- rep(0, NCOL(reco))
    id_nn[bts] <- 1
  }

  rowid <- which(rowSums(reco<0)>0)
  base0 <- base[rowid, , drop = FALSE]
  reco0 <- reco[rowid, , drop = FALSE]

  #controls
  nb <- NCOL(strc_mat)

  tol <- ifelse(is.null(settings$tol), sqrt(.Machine$double.eps), settings$tol)
  pbar <- ifelse(is.null(settings$pbar), 10, settings$pbar)
  ptype <- ifelse(is.null(settings$ptype), "fixed", settings$ptype)
  gtol <- ifelse(is.null(settings$gtol), sqrt(.Machine$double.eps), settings$gtol)
  itmax <- ifelse(is.null(settings$itmax), 100, settings$itmax)
  ninf <- nb + 1
  if(ptype == "fixed") {
    alpha <- 1:nb
  }else if(ptype == "random") {
    alpha <- sample(1:nb, nb, replace = FALSE)
  }
  maxp <- pbar

  z <- t(solve(cov_mat, strc_mat))
  grad_all <- -t(strc_mat) %*% solve(cov_mat, t(base0))
  grad0 <- z %*% t(reco0) + grad_all

  bpv_step <- lapply(1:length(rowid), function(j){
    start <- Sys.time()
    baseh <- base0[j, ]
    grad <- grad0[,j]

    # Starting set
    Fset <- 1:nb
    Gset <- numeric(0)

    rbts <- reco0[j, id_nn == 1]
    gradg <- grad[Gset]

    # To avoid nondegenerate problem (as done by Jason Cantarella)
    rbts[abs(rbts) < tol] <- 0L
    gradg[abs(gradg) < tol] <- 0L
    verb <- TRUE
    for(i in 1:itmax){
      i1 <- Fset[which(rbts < -tol)]
      i2 <- Gset[which(gradg < -tol)]
      ivec <- union(i1, i2)

      if(length(ivec) < ninf) {
        ninf <- length(ivec)
        maxp <- pbar
      }else if(maxp >= 1) {
        maxp <- maxp - 1
      }else{
        if(ptype == "fixed") {
          if(verb){
            cat("Slow zone: it might take some time to converge! \n")
            verb <- FALSE
          }

          r <- max(ivec)
        }else if(ptype == "random") {
          if(verb){
            cat("Slow zone: it might take some time to converge! \n")
            verb <- FALSE
          }
          r <- alpha[max(which(alpha %in% ivec))]
        }
        if(is.element(r, i1)) {
          i1 <- r
          i2 <- numeric(0)
        }else{
          i1 <- numeric(0)
          i2 <- r
        }
      }

      # updating f and g
      Fset <- union(Fset[!Fset %in% i1], i2)
      Gset <- union(Gset[!Gset %in% i2], i1)

      if(length(Gset) == 0) {
        rfc <- reco0[j,]
      }else{
        strc_mat_0 <- strc_mat[, -Gset, drop = FALSE]
        if(approach == "proj"){
          id_c <- which(id_nn==1)[Gset]
          idb0 <- which(id_nn==1)[-Gset]
          ida0 <- (1:NROW(strc_mat_0))[-idb0]
          id0 <- c(ida0, idb0)
          agg_mat0 <- strc_mat[-idb0, -Gset, drop = FALSE]
          cons_mat_0 <- cbind(Diagonal(NROW(agg_mat0)), -agg_mat0)

          tmp <- reco.proj(base = t(baseh[id0]), cons_mat = cons_mat_0,
                                     cov_mat = cov_mat[id0, id0])
          tmp <- tmp[sort.int(id0, index.return = T)$ix]
        }else{
          tmp2 <- reco.strc(base = t(baseh), strc_mat = strc_mat_0, cov_mat = cov_mat)
        }
        rfc <- as.numeric(tmp)
      }

      grad <- as.matrix(z %*% rfc) + grad_all[, j]
      rbts <-  rfc[id_nn == 1][Fset]
      gradg <- grad[Gset]

      # To avoid nondegenerate problem (as done by Jason Cantarella)
      rbts[abs(rbts) < tol] <- 0L
      gradg[abs(gradg) < tol] <- 0L
      if(all(rbts > -gtol) & all(gradg > -gtol)){
        flag <- 1
        break
      }else{
        flag <- -2
      }
    }

    if(flag == 1){
      rfc[rfc <= sqrt(.Machine$double.eps)] <- 0
      if(i == itmax){
        flag <- 2
      }
    }
    end <- Sys.time()
    out <- list()
    out$reco <- rfc
    out$info <- c(difftime(end,start, units = "secs"), i, flag)
    out$idx <- Gset
    out
  })

  bpv_step <- do.call("rbind", bpv_step)

  # Point reconciled forecasts
  reco[rowid, ] <- do.call("rbind", bpv_step[, "reco"])

  info <- do.call("rbind", bpv_step[, "info"])
  colnames(info) <- c("run_time", "iter", "status")
  rownames(info) <- rowid
  attr(reco, "info") <- info
  return(reco)
}

reco.sftb <- function(base, reco, strc_mat, id_nn = NULL, bounds = NULL, ...){
  # Check input
  if(missing(strc_mat)){
    cli_abort("Mandatory arguments: {.arg strc_mat} and {.arg cov_mat}.",
              call = NULL)
  }

  if(missing(reco)){
    reco <- base
  }

  if(is.null(bounds)){
    return(reco)
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

  nbid <- bounds[,1,drop = TRUE]
  reco <- t(apply(reco, 1, function(x){
    id <- x[nbid]<=bounds[,2,drop = TRUE]+sqrt(.Machine$double.eps)
    x[nbid][id] <- bounds[,2,drop = TRUE][id]

    id <- x[nbid]>=bounds[, 3, drop = TRUE]-sqrt(.Machine$double.eps)
    x[nbid][id] <- bounds[, 3, drop = TRUE][id]
    x
  }))

  bts <- reco[, id_nn == 1, drop = FALSE]
  as.matrix(bts %*% t(strc_mat))
}
