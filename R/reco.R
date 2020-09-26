# basef: base forecasts (h x col)
#    W: covariance
#    H: Ht in ctr, Zt in thf e Ut in hts
recoM <- function(basef, W, H, sol = "direct", nn = FALSE, settings, b_pos = NULL, keep = "list") {
  sol <- match.arg(sol, c("direct", "osqp"))

  switch(sol,
    direct = {
      out <- list()
      c_W <- ncol(W)

      lm_dx <- methods::as(H %*% t(basef), "CsparseMatrix")
      lm_sx <- Matrix::Matrix(H %*% W %*% t(H), sparse = TRUE, forceCheck = TRUE)
      out$recf <- t(.sparseDiagonal(NCOL(W)) %*% t(basef) - W %*% t(H) %*% solveLin(lm_sx, lm_dx))

      if (nn) {
        if (any(out$recf < (-1e-6))) {
          if (isDiagonal(W)) {
            P <- .sparseDiagonal(x = diag(W)^(-1))
          } else {
            R <- chol(W)
            P <- chol2inv(R)
          }

          rec <- apply(basef, 1, function(x) {
            M_osqp(
              y = x, H = H, P = P, nn = nn,
              settings = settings, b_pos = b_pos
            )
          })

          out <- list()
          out$recf <- do.call("rbind", extract(rec, 1))
          out$info <- do.call("rbind", extract(rec, 2))
          colnames(out$info) <- c("obj_val", "run_time", "iter", "pri_res", "status", "status_polish")

          if (keep == "list") {
            out$W <- as.matrix(W)
          }
        } else {
          out$recf <- out$recf * (out$recf > 0)

          if (keep == "list") {
            M <- .sparseDiagonal(c_W) - W %*% t(H) %*% solveLin(H %*% W %*% t(H), H)
            out$varf <- diag(M %*% W)
            out$M <- as.matrix(M)
            out$W <- as.matrix(W)
          }
        }
      } else if (keep == "list") {
        M <- .sparseDiagonal(c_W) - W %*% t(H) %*% solveLin(H %*% W %*% t(H), H)
        out$varf <- diag(M %*% W)
        out$M <- as.matrix(M)
        out$W <- as.matrix(W)
      }
    },
    osqp = {
      if (isDiagonal(W)) {
        P <- .sparseDiagonal(x = diag(W)^(-1))
      } else {
        R <- chol(W)
        P <- chol2inv(R)
      }

      rec <- apply(basef, 1, function(x) {
        M_osqp(
          y = x, H = H, P = P, nn = nn,
          settings = settings, b_pos = b_pos
        )
      })

      out <- list()
      out$recf <- do.call("rbind", extract(rec, 1))
      out$info <- do.call("rbind", extract(rec, 2))
      colnames(out$info) <- c("obj_val", "run_time", "iter", "pri_res", "status", "status_polish")

      if (keep == "list") {
        out$W <- as.matrix(W)
      }
    }
  )
  return(out)
}

recoS <- function(basef, W, S, sol = "direct", nn = FALSE, settings, b_pos = NULL, keep = "list") {
  sol <- match.arg(sol, c("direct", "osqp"))

  switch(sol,
    direct = {
      out <- list()

      if (isDiagonal(W)) {
        Wm1 <- .sparseDiagonal(x = diag(W)^(-1))
        lm_dx1 <- methods::as(t(S) %*% Wm1 %*% t(basef), "CsparseMatrix")
        lm_sx1 <- t(S) %*% Wm1 %*% S
        out$recf <- t(S %*% solveLin(lm_sx1, lm_dx1))

        if (nn) {
          if (any(out$recf < (-1e-6))) {
            P <- t(S) %*% Wm1 %*% S
            q <- (-1) * t(Wm1 %*% S)
            rec <- apply(basef, 1, function(x) {
              S_osqp(
                y = x, S = S, P = P, q = q, nn = nn,
                settings = settings, b_pos = b_pos
              )
            })

            out <- list()
            out$recf <- do.call("rbind", extract(rec, 1))
            out$info <- do.call("rbind", extract(rec, 2))
            colnames(out$info) <- c("obj_val", "run_time", "iter", "pri_res", "status", "status_polish")

            if (keep == "list") {
              out$W <- W
            }
          } else {
            out$recf <- out$recf * (out$recf > 0)
            if (keep == "list") {
              lm_dx2 <- methods::as(t(S) %*% Wm1, "CsparseMatrix")
              G <- S %*% solveLin(lm_sx1, lm_dx2)
              out$varf <- diag(G %*% W)
              out$G <- as.matrix(G)
              out$W <- W
            }
          }
        } else if (keep == "list") {
          lm_dx2 <- methods::as(t(S) %*% Wm1, "CsparseMatrix")
          G <- S %*% solveLin(lm_sx1, lm_dx2)
          out$varf <- diag(G %*% W)
          out$G <- as.matrix(G)
          out$W <- W
        }
      } else {
        Q <- solveLin(W, S)

        lm_dx1 <- methods::as(t(Q) %*% t(basef), "CsparseMatrix")
        lm_sx1 <- t(S) %*% Q
        out$recf2 <- t(S %*% solveLin(lm_sx1, lm_dx1))

        if (nn) {
          if (any(out$recf < (-1e-6))) {
            P <- t(S) %*% Q
            q <- (-1) * t(Q)
            rec <- apply(basef, 1, function(x) {
              S_osqp(
                y = x, S = S, P = P, q = q, nn = nn,
                settings = settings, b_pos = b_pos
              )
            })

            out <- list()
            out$recf <- do.call("rbind", extract(rec, 1))
            out$info <- do.call("rbind", extract(rec, 2))
            colnames(out$info) <- c("obj_val", "run_time", "iter", "pri_res", "status", "status_polish")

            if (keep == "list") {
              out$W <- W
            }
          } else {
            out$recf <- out$recf * (out$recf > 0)
            if (keep == "list") {
              G <- S %*% solveLin(lm_sx1, lm_dx2)
              out$varf <- diag(G %*% W)
              out$G <- as.matrix(G)
              out$W <- W
            }
          }
        } else if (keep == "list") {
          G <- S %*% solveLin(lm_sx1, lm_dx2)
          out$varf <- diag(G %*% W)
          out$G <- as.matrix(G)
          out$W <- W
        }
      }
    },
    osqp = {
      if (isDiagonal(W)) {
        Q <- .sparseDiagonal(x = diag(W)^(-1))
        P <- t(S) %*% Q %*% S
        q <- (-1) * t(Q %*% S)
      } else {
        Q <- solveLin(W, S)
        P <- t(S) %*% Q
        q <- (-1) * t(Q)
      }

      rec <- apply(basef, 1, function(x) {
        S_osqp(
          y = x, S = S, q = q, P = P, nn = nn,
          settings = settings, b_pos = b_pos
        )
      })

      out <- list()
      out$recf <- do.call("rbind", extract(rec, 1))
      out$info <- do.call("rbind", extract(rec, 2))
      colnames(out$info) <- c("obj_val", "run_time", "iter", "pri_res", "status", "status_polish")

      if (keep == "list") {
        out$W <- W
      }
    }
  )

  return(out)
}

# y = riga di basef
M_osqp <- function(y, P, H, nn, settings, b_pos) {
  c <- ncol(H)
  r <- nrow(H)
  l <- rep(0, r)
  u <- rep(0, r)
  A <- H

  q <- (-1) * t(P) %*% as.vector(y)

  if (nn) {
    A <- rbind(A, Diagonal(c)[b_pos == 1, ])
    l <- c(l, rep(0, sum(b_pos)))
    u <- c(u, rep(Inf, sum(b_pos)))
  }

  rec <- solve_osqp(P, q, A, l, u, settings)

  out <- list()
  out$recf <- rec$x

  if (rec$info$status_val != 1 | rec$info$pri_res > 1e-6) {
    warning(paste("OSQP flag", rec$info$status_val, "OSQP pri_res", rec$info$pri_res), call. = FALSE)
  }

  if (nn) {
    out$recf[which(out$recf < 0 & out$recf > -1e-6)] <- 0
  }

  out$info <- c(
    rec$info$obj_val, rec$info$run_time, rec$info$iter, rec$info$pri_res,
    rec$info$status_val, rec$info$status_polish
  )

  return(out)
}

extract <- function(l, num) {
  do.call("rbind", l)[, num]
}

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
S_osqp <- function(y, q, P, S, nn, settings, b_pos) {
  q <- q %*% y
  r <- nrow(S)
  c <- sum(b_pos)
  A <- NULL
  l <- NULL
  u <- NULL

  if (nn) {
    A <- .sparseDiagonal(c)
    l <- rep(0, sum(b_pos))
    u <- rep(Inf, sum(b_pos))
  }

  rec <- solve_osqp(P, q, A, l, u, settings)

  out <- list()
  out$recf <- as.vector(S %*% rec$x)

  if (rec$info$status_val != 1 | rec$info$pri_res > 1e-6) {
    warning(paste("OSQP flag", rec$info$status_val, "OSQP pri_res", rec$info$pri_res), call. = FALSE)
  }

  if (nn) {
    out$recf[which(out$recf < 0 & out$recf > -1e-6)] <- 0
  }

  out$info <- c(
    rec$info$obj_val, rec$info$run_time, rec$info$iter,
    rec$info$pri_res, rec$info$status_val, rec$info$status_polish
  )
  return(out)
}

solveLin <- function(msx, mdx) {
  tryCatch(solve(msx, mdx), error = function(cond) {

    # browser()
    warning(
      "An error in LU decomposition occurred, with the following message:\n",
      cond$message, "\n Trying QR decomposition instead..."
    )
    tryCatch(solve(qr(msx), mdx), error = function(cond) {

      # browser()
      warning(
        "An error in QR decomposition occurred, with the following message:\n",
        cond$message, "\n Trying chol decomposition instead..."
      )
      backsolve(chol(msx), mdx)
    })
  })
}
