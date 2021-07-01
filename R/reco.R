# basef: base forecasts (h x col)
#    W: covariance
#    H: Ht in ctr, Zt in thf e Ut in hts
recoM <- function(basef, W, Ht, sol = "direct", nn = FALSE, nn_type = "osqp", settings, b_pos = NULL,
                  keep = "list", bounds = NULL, S) {

  sol <- match.arg(sol, c("direct", "osqp"))
  nn_type <- match.arg(nn_type, c("osqp", "fbpp", "KAnn", "sntz"))

  if(!is.null(bounds)){
    if(!is.matrix(bounds) | NCOL(bounds) != 2 | NROW(bounds) != NCOL(basef)){
      stop("bounds must be a matrix (", NCOL(basef), "x2)", call. = FALSE)
    } else {
      sol <- "osqp"
    }
  }

  switch(sol,
         direct = {
           out <- list()

           lm_dx <- methods::as(Ht %*% t(basef), "CsparseMatrix")
           lm_sx <- Matrix::Matrix(Ht %*% W %*% t(Ht), sparse = TRUE, forceCheck = TRUE)
           out$recf <- t(.sparseDiagonal(NCOL(W)) %*% t(basef) - W %*% t(Ht) %*% solveLin(lm_sx, lm_dx))

           if(nn & any(out$recf < (-sqrt(.Machine$double.eps)))){
             switch(nn_type,
                    osqp = {
                      if(isDiagonal(W)){
                        P <- .sparseDiagonal(x = diag(W)^(-1))
                      }else{
                        R <- chol(W)
                        P <- chol2inv(R)
                      }
                      id <- which(rowSums(out$recf<(-sqrt(.Machine$double.eps)))!=0)
                      rec <- apply(basef[id,,drop = FALSE], 1, function(x){
                        M_osqp(y = x, Ht = Ht, P = P, nn = nn,
                               bounds = bounds, settings = settings, b_pos = b_pos)
                      })
                      rec <- do.call("rbind",rec)
                      out$recf[id,] <- do.call("rbind", rec[,"recf"])
                      out$info <- do.call("rbind", rec[,"info"])
                      colnames(out$info) <- c("obj_val", "run_time", "iter", "pri_res",
                                              "status", "status_polish")
                      rownames(out$info) <- id

                      out$recf[-id,,drop=FALSE] <- out$recf[-id,,drop=FALSE] * (out$recf[-id,,drop=FALSE] > 0)

                      if (keep == "list") {
                        out$W <- as.matrix(W)
                      }
                    },
                    fbpp = {
                      param <- list(tol = 1e-6,
                                    itmax = 10)
                      param[names(settings)] <- settings

                      rec <- apply(out$recf, 1, function(x) {
                        M_fbpp(y = x, Ht = Ht, W = W, b_pos = b_pos,
                               param = param, keep = keep)
                      })

                      rec <- do.call("rbind",rec)
                      out$recf <- do.call("rbind", rec[,"recf"])
                      out$idx <- as.vector(stats::na.omit(unlist(rec[,"idx"])))
                      out$info <- do.call("rbind", rec[,"info"])
                      colnames(out$info) <- c("run_time (s)", "iter", "status", "tol", "itmax", "lidx")
                      rownames(out$info) <- 1:NROW(out$recf)
                      out$info <- out$info[out$info[,2]!=0, , drop=FALSE]

                      if(keep == "list"){
                        out$varf <- do.call("rbind", rec[,"varf"])
                        out$M <- unlist(rec[,"M"])
                        if(length(out$M)){
                          out$M <- out$M[[1]]
                        }
                        out$W <- as.matrix(W)
                      }
                    },
                    KAnn = {
                      param <- list(tol = 1e-6,
                                    itmax = 10)
                      param[names(settings)] <- settings

                      rec <- apply(out$recf, 1, function(x) {
                        M_KAnn(y = x, Ht = Ht, W = W, b_pos = b_pos,
                               param = param, keep = keep)
                      })

                      rec <- do.call("rbind",rec)
                      out$recf <- do.call("rbind", rec[,"recf"])
                      out$info <- do.call("rbind", rec[,"info"])
                      colnames(out$info) <- c("run_time (s)", "iter", "status", "tol", "itmax")
                      rownames(out$info) <- 1:NROW(out$recf)
                      out$info <- out$info[out$info[,2]!=0, , drop=FALSE]

                      if(keep == "list"){
                        out$varf <- do.call("rbind", rec[,"varf"])
                        out$M <- unlist(rec[,"M"])
                        if(length(out$M)){
                          out$M <- out$M[[1]]
                        }
                        out$W <- as.matrix(W)
                      }
                    },
                    sntz = {
                      na <- NROW(S)-NCOL(S)
                      bottom <- out$recf[,-c(1:na),drop = FALSE]
                      bottom <- bottom*(bottom > 0)
                      out$recf <- bottom %*% t(S)

                      if (keep == "list") {
                        out$W <- as.matrix(W)
                      }
                    })
           }else{
             if(nn){
               out$recf <- out$recf * (out$recf > 0)
             }

             if(keep == "list"){
               M <- unname(.sparseDiagonal(NCOL(W)) - W %*% t(Ht) %*% solveLin(Ht %*% W %*% t(Ht), Ht))
               out$varf <- diag(M %*% W)
               out$M <- as.matrix(M)
               out$W <- as.matrix(W)
             }
           }
         },
         osqp = {
           if(isDiagonal(W)){
             P <- .sparseDiagonal(x = diag(W)^(-1))
           }else{
             R <- chol(W)
             P <- chol2inv(R)
           }

           rec <- apply(basef, 1, function(x){
             M_osqp(y = x, Ht = Ht, P = P, nn = nn, bounds = bounds,
                    settings = settings, b_pos = b_pos)
           })
           rec <- do.call("rbind",rec)
           out <- list()
           out$recf <- do.call("rbind", rec[,"recf"])
           out$info <- do.call("rbind", rec[,"info"])
           colnames(out$info) <- c("obj_val", "run_time", "iter", "pri_res",
                                   "status", "status_polish")
           rownames(out$info) <- 1:NROW(out$recf)

           if(keep == "list"){
             out$W <- as.matrix(W)
           }
         }
  )
  return(out)
}

recoS <- function(basef, W, S, sol = "direct", nn = FALSE, settings, b_pos = NULL,
                  keep = "list", bounds = NULL, nn_type = "osqp") {
  sol <- match.arg(sol, c("direct", "osqp"))


  nn_type <- match.arg(nn_type, c("osqp", "sntz"))

  if(!is.null(bounds)){
    if(!is.matrix(bounds) | NCOL(bounds) != 2 | NROW(bounds) != NCOL(basef)){
      stop("bounds must be a matrix (", NCOL(basef), "x2)", call. = FALSE)
    } else {
      sol <- "osqp"
    }
  }

  switch(sol,
         direct = {
           out <- list()

           if(isDiagonal(W)){
             Wm1 <- .sparseDiagonal(x = diag(W)^(-1))
             lm_dx1 <- methods::as(t(S) %*% Wm1 %*% t(basef), "CsparseMatrix")
             lm_sx1 <- t(S) %*% Wm1 %*% S
             out$recf <- t(S %*% solveLin(lm_sx1, lm_dx1))

             if(nn & any(out$recf < (-sqrt(.Machine$double.eps)))){
               switch(nn_type,
                      osqp = {
                        P <- t(S) %*% Wm1 %*% S
                        q <- (-1) * t(Wm1 %*% S)

                        id <- which(rowSums(out$recf<(-sqrt(.Machine$double.eps)))!=0)
                        rec <- apply(basef[id,,drop = FALSE], 1, function(x) {
                          S_osqp(y = x, S = S, P = P, q = q, nn = nn, bounds = bounds,
                                 settings = settings, b_pos = b_pos)
                        })
                        rec <- do.call("rbind",rec)
                        out$recf[id,] <- do.call("rbind", rec[,"recf"])
                        out$info <- do.call("rbind", rec[,"info"])
                        colnames(out$info) <- c("obj_val", "run_time", "iter", "pri_res",
                                                "status", "status_polish")
                        rownames(out$info) <- id

                        out$recf[-id,,drop=FALSE] <- out$recf[-id,,drop=FALSE] * (out$recf[-id,,drop=FALSE] > 0)
                      },
                      sntz = {
                        na <- NROW(S)-NCOL(S)
                        bottom <- out$recf[,-c(1:na),drop = FALSE]
                        bottom <- bottom*(bottom > 0)
                        out$recf <- bottom %*% t(S)
                      })

               if (keep == "list") {
                 out$W <- as.matrix(W)
               }
             }else{
               if(nn){
                 out$recf <- out$recf * (out$recf > 0)
               }
               if(keep == "list"){
                 lm_dx2 <- methods::as(t(S) %*% Wm1, "CsparseMatrix")
                 G <- unname(solveLin(lm_sx1, lm_dx2))
                 M <- S %*% G
                 out$varf <- diag(M %*% W)
                 out$G <- as.matrix(G)
                 out$M <- as.matrix(M)
                 out$W <- W
               }
             }
           }else{
             Q <- solveLin(W, S)

             lm_dx1 <- methods::as(t(Q) %*% t(basef), "CsparseMatrix")
             lm_sx1 <- t(S) %*% Q
             out$recf <- t(S %*% solveLin(lm_sx1, lm_dx1))

             if(nn & any(out$recf < (-sqrt(.Machine$double.eps)))){
               switch(nn_type,
                      osqp = {
                        P <- t(S) %*% Q
                        q <- (-1) * t(Q)
                        id <- which(rowSums(out$recf<(-sqrt(.Machine$double.eps)))!=0)
                        rec <- apply(basef[id, , drop = FALSE], 1, function(x) {
                          S_osqp(y = x, S = S, P = P, q = q, nn = nn, bounds = bounds,
                                 settings = settings, b_pos = b_pos)
                        })
                        rec <- do.call("rbind",rec)
                        out$recf[id,] <- do.call("rbind", rec[,"recf"])
                        out$info <- do.call("rbind", rec[,"info"])
                        colnames(out$info) <- c("obj_val", "run_time", "iter", "pri_res",
                                                "status", "status_polish")
                        rownames(out$info) <- id
                      },
                      sntz = {
                        na <- NROW(S)-NCOL(S)
                        bottom <- out$recf[,-c(1:na),drop = FALSE]
                        bottom <- bottom*(bottom > 0)
                        out$recf <- bottom %*% t(S)
                      })

               if (keep == "list") {
                 out$W <- as.matrix(W)
               }
             }else{
               if(nn){
                 out$recf <- out$recf * (out$recf > 0)
               }
               if(keep == "list"){
                 G <- unname(solveLin(lm_sx1, t(Q)))
                 M <- S %*% G
                 out$varf <- diag(M %*% W)
                 out$G <- as.matrix(G)
                 out$M <- as.matrix(M)
                 out$W <- W
               }
             }
           }
         },
         osqp = {
           if(isDiagonal(W)){
             Q <- .sparseDiagonal(x = diag(W)^(-1))
             P <- t(S) %*% Q %*% S
             q <- (-1) * t(Q %*% S)
           }else{
             Q <- solveLin(W, S)
             P <- t(S) %*% Q
             q <- (-1) * t(Q)
           }

           rec <- apply(basef, 1, function(x) {
             S_osqp(y = x, S = S, q = q, P = P, nn = nn, bounds = bounds,
                    settings = settings, b_pos = b_pos)
           })
           rec <- do.call("rbind",rec)
           out <- list()
           out$recf <- do.call("rbind", rec[,"recf"])
           out$info <- do.call("rbind", rec[,"info"])
           colnames(out$info) <- c("obj_val", "run_time", "iter", "pri_res",
                                   "status", "status_polish")
           rownames(out$info) <- 1:NROW(out$recf)

           if(keep == "list"){
             out$W <- W
           }
         }
  )

  return(out)
}

# y = riga di basef
M_osqp <- function(y, P = NULL, Ht = NULL, nn = NULL, settings = NULL,
                   b_pos = NULL, bounds = NULL) {
  c <- ncol(Ht)
  r <- nrow(Ht)
  # Linear constrains H = 0
  l <- rep(0, r)
  u <- rep(0, r)
  A <- Ht

  q <- (-1) * t(P) %*% as.vector(y)

  if(nn){
    # bts >= 0
    A <- rbind(A, .sparseDiagonal(c)[b_pos == 1, ])
    l <- c(l, rep(0, sum(b_pos)))
    u <- c(u, rep(Inf, sum(b_pos)))
  }

  if(!is.null(bounds)){
    bounds_rows <- rowSums(abs(bounds) == Inf) < 2
    A <- rbind(A, .sparseDiagonal(c)[bounds_rows, ])
    l <- c(l, bounds[bounds_rows, 1, drop = TRUE])
    u <- c(u, bounds[bounds_rows, 2, drop = TRUE])
  }

  if(length(settings)==0){
    settings = osqpSettings(verbose = FALSE,
                            eps_abs = 1e-5,
                            eps_rel = 1e-5,
                            polish_refine_iter = 100,
                            polish=TRUE)
  }

  rec <- solve_osqp(P, q, A, l, u, settings)

  out <- list()
  out$recf <- rec$x

  if(rec$info$status_val != 1){
    warning(paste("OSQP flag", rec$info$status_val, "OSQP pri_res", rec$info$pri_res), call. = FALSE)
  }

  if(nn){
    out$recf[which(out$recf < 0)] <- 0
  }

  out$info <- c(rec$info$obj_val, rec$info$run_time, rec$info$iter, rec$info$pri_res,
                rec$info$status_val, rec$info$status_polish)

  return(out)
}

#extract <- function(l, num) {
#  do.call("rbind", l)[, num]
#}

# S_sol <- function(basef, W, S){
#   out <- list()
#   c <- ncol(W)
#
#   if(isDiagonal(W)){
#     Wm1 <- .sparseDiagonal(x = diag(W)^(-1))
#     G <- S%*%solve(t(S)%*%Wm1%*%S)%*%t(S)%*%Wm1
#   }else{
#     R <- chol(W)
#     Q <- t(solve(R))
#
#     Sstar <- Q%*%S
#     P <- t(Sstar)%*%Sstar
#
#     G <- S%*%solve(t(Sstar)%*%Sstar)%*%t(Sstar)%*%Q
#   }
#
#   out$recf <- basef%*%t(G)
#   out$G <- G
#   out$W <- Matrix(W)
#   return(out)
# }

# y = riga di basef
S_osqp <- function(y, q = NULL, P = NULL, S = NULL, nn = NULL,
                   settings = NULL, b_pos = NULL, bounds = NULL) {
  q <- q %*% y
  r <- nrow(S)
  c <- sum(b_pos)
  A <- NULL
  l <- NULL
  u <- NULL

  if(nn){
    A <- .sparseDiagonal(c)
    l <- rep(0, sum(b_pos))
    u <- rep(Inf, sum(b_pos))
  }

  if(!is.null(bounds)){
    bounds_rows <- rowSums(abs(bounds) == Inf) < 2
    A <- .sparseDiagonal(c)[bounds_rows, ]
    l <- bounds[bounds_rows, 1, drop = TRUE]
    u <- bounds[bounds_rows, 2, drop = TRUE]
  }

  if(length(settings)==0){
    settings <- osqpSettings(verbose = FALSE,
                             eps_abs = 1e-5,
                             eps_rel = 1e-5,
                             polish_refine_iter = 100,
                             polish=TRUE)
  }

  rec <- solve_osqp(P, q, A, l, u, settings)

  out <- list()
  out$recf <- as.vector(S %*% rec$x)

  if(rec$info$status_val != 1){
    warning(paste("OSQP flag", rec$info$status_val, "OSQP pri_res", rec$info$pri_res), call. = FALSE)
  }

  if(nn){
    out$recf[which(out$recf < 0)] <- 0
  }

  out$info <- c(rec$info$obj_val, rec$info$run_time, rec$info$iter,
                rec$info$pri_res, rec$info$status_val, rec$info$status_polish)
  return(out)
}

M_fbpp <- function(y, Ht, W, b_pos, param, keep){
  tol <- param$tol
  itmax <- param$itmax
  if(all(y>(-tol))){
    out <- list()
    out$recf <- y*(y>0)
    out$info <- c(0, 0, tol, 0, itmax, 0)
    out$idx <- NA

    if(keep == "list"){
      out$M <- .sparseDiagonal(NCOL(W)) - W %*% t(Ht) %*% solveLin(Ht %*% W %*% t(Ht), Ht)
      out$varf <- diag(out$M%*%W)
    }
    return(out)
  }
  start <- Sys.time()
  bid <- Position(function(x) x == 1, b_pos)
  idx <- c()
  id <- list()
  Htnn <- Ht
  if(keep == "list"){
    #M <- list()
    #M[[1]] <- .sparseDiagonal(NCOL(W)) - W %*% t(Ht) %*% solveLin(Ht %*% W %*% t(Ht), Ht)
    M <- .sparseDiagonal(NCOL(W)) - W %*% t(Ht) %*% solveLin(Ht %*% W %*% t(Ht), Ht)
  }

  for(i in 1:itmax){
    id[[i]] <- which(y < -tol)
    idx <- unique(c(idx, id[[i]][id[[i]]>=bid]))
    block <- sparseMatrix(i = 1:length(idx),
                          j = idx,
                          x = 1,
                          dims = c(length(idx), NCOL(Htnn)))
    Htnn <- rbind(Ht, block)
    lm_dx <- methods::as(Htnn %*% y, "CsparseMatrix")
    lm_sx <- Matrix::Matrix(Htnn %*% W %*% t(Htnn), sparse = TRUE, forceCheck = TRUE)
    y <- as.vector(.sparseDiagonal(NCOL(W)) %*% y - W %*% t(Htnn) %*% solveLin(lm_sx, lm_dx, verb = FALSE))

    if(keep == "list"){
      M2 <- .sparseDiagonal(NCOL(W)) - W %*% t(Htnn) %*% solveLin(lm_sx, Htnn, verb = FALSE)
      M <- M2%*%M
      #M[[i+1]] <- M2
    }

    if(all(y >= (-tol))){
      flag <- 1
      break
    }else{
      flag <- -2
    }
  }

  if(flag == 1){
    y <- (y > 0)*y + (y < 0) * 0
    if(i == itmax){
      flag <- 2
    }
  }

  end <- Sys.time()
  out <- list()
  out$recf <- y
  out$info <- c(difftime(end,start, units = "secs"), i, flag, tol, itmax, length(idx))
  out$idx <- id

  if(keep == "list"){
    out$M <- M
    out$varf <- diag(M%*%W)
  }
  return(out)
}


M_KAnn <- function(y, Ht, W, b_pos, param, keep){
  tol <- param$tol
  itmax <- param$itmax
  if(all(y>(-tol))){
    out <- list()
    out$recf <- y*(y>0)
    out$info <- c(0, 0, tol, 0, itmax)

    if(keep == "list"){
      out$M <- .sparseDiagonal(NCOL(W)) - W %*% t(Ht) %*% solveLin(Ht %*% W %*% t(Ht), Ht)
      out$varf <- diag(out$M%*%W)
    }
    return(out)
  }
  start <- Sys.time()
  lm_sx <- Matrix::Matrix(Ht %*% W %*% t(Ht), sparse = TRUE, forceCheck = TRUE)
  if(keep == "list"){
    M <- .sparseDiagonal(NCOL(W)) - W %*% t(Ht) %*% solveLin(lm_sx, Ht)
  }
  id <- list()
  for(i in 1:itmax){
    id[[i]] <- (y > -tol)
    yc <- y * id[[i]]

    lm_dx <- methods::as(Ht %*% yc, "CsparseMatrix")
    y <- as.vector(.sparseDiagonal(NCOL(W)) %*% yc - W %*% t(Ht) %*% solveLin(lm_sx, lm_dx, verb = FALSE))

    if(keep == "list"){
      M1 <- M%*%diag(id[[i]])
      M <- M1%*%M
    }

    if(all(y >= (-tol))){
      flag <- 1
      break
    }else{
      flag <- -2
    }
  }

  if(flag == 1){
    y <- (y > 0)*y + (y < 0) * 0
    if(i == itmax){
      flag <- 2
    }
  }

  end <- Sys.time()
  out <- list()
  out$recf <- y
  out$info <- c(difftime(end,start, units = "secs"), i, flag, tol, itmax)

  if(keep == "list"){
    out$M <- M
    out$varf <- diag(M%*%W)
  }
  return(out)
}

solveLin <- function(msx, mdx, verb = FALSE) {
  tryCatch(solve(msx, mdx), error = function(cond){
    if(verb){
      warning(
        "An error in LU decomposition occurred, with the following message:\n",
        cond$message, "\n Trying QR decomposition instead...", call. = FALSE)
    }
    tryCatch(solve(qr(msx), mdx), error = function(cond){
      if(verb){
        warning(
          "An error in QR decomposition occurred, with the following message:\n",
          cond$message, "\n Trying chol decomposition instead...", call. = FALSE)
      }
      backsolve(chol(msx), mdx)
    })
  })
}
