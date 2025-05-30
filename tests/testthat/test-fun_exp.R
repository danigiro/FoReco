# test fun_export.R
if(require(testthat)){
  test_that("Commutation matrix", {
    Y <- matrix(rnorm(30), 5, 6)
    P <- commat(5, 6)
    expect_equal(as.vector(P %*% as.vector(Y)), as.vector(t(Y)))
  })

  test_that("Balanced and unbalanced hierarchy", {
    A <- matrix(c(1,1,1,1,
                  1,1,0,0), 2, byrow = TRUE)
    obj <- balance_hierarchy(A, nodes = "auto", sparse = FALSE)
    expect_equal(obj,
                 balance_hierarchy(A, nodes = c(1, 1), sparse = FALSE))

    expect_equal(as.matrix(unbalance_hierarchy(obj$bam)), A)


    A <- matrix(c(1,1,1,1,1,1,
                  1,1,0,0,0,0,
                  0,0,0,1,1,0,
                  1,1,1,0,0,0), 4, byrow = TRUE)
    obj <- balance_hierarchy(A, nodes = "auto", sparse = FALSE)
    expect_equal(unbalance_hierarchy(obj$bam, sparse = FALSE), A)
  })


  test_that("LCMAT", {
    C <- rbind(c(1,-1,-1, 0, 0, 0, 0),
               c(0, 1, 0,-1,-1, 0, 0),
               c(0, 0, 1, 0, 0,-1,-1))
    A <- matrix(c(1,1,1,1,
                  1,1,0,0,
                  0,0,1,1), 3, byrow = TRUE)

    r1 <- lcmat(C, method = "qr", sparse = FALSE)
    r2 <- lcmat(C, method = "rref", sparse = FALSE)

    expect_equal(r1, r2)
    expect_equal(r1$agg_mat, A)
    expect_equal(r1$pivot, 1:NCOL(C))
  })

  test_that("df2aggmat", {
    data_bts <- data.frame(X1 = c("A", "A", "B", "B", "B"),
                           X2 = c("A", "B", "A", "B", "C"),
                           stringsAsFactors = FALSE)
    A <- matrix(c(1,1,1,1,1,
                  1,1,0,0,0,
                  0,0,1,1,1), 3, byrow = TRUE)
    agg_mat <- df2aggmat(~ X1 / X2, data_bts, sep = "",
                         verbose = FALSE, sparse = FALSE)
    expect_equal(agg_mat, A, ignore_attr = TRUE)
  })

  test_that("shr", {
    # Original function: https://github.com/earowang/hts/blob/master/R/MinT.R
    shrink.estim <- function(x, tar){
      if (is.matrix(x) == TRUE && is.numeric(x) == FALSE)
        stop("The data matrix must be numeric!", call. = FALSE)
      p <- ncol(x)
      n <- nrow(x)
      covm <- crossprod(x) / n
      corm <- cov2cor(covm)
      xs <- scale(x, center = FALSE, scale = sqrt(diag(covm)))
      v <- (1/(n * (n - 1))) * (crossprod(xs^2) - 1/n * (crossprod(xs))^2)
      diag(v) <- 0
      corapn <- cov2cor(tar)
      d <- (corm - corapn)^2
      lambda <- sum(v)/sum(d)
      lambda <- max(min(lambda, 1), 0)
      shrink.cov <- lambda * tar + (1 - lambda) * covm
      return(list(shrink.cov, lambda))
    }

    N <- 500
    n <- 2
    mean <- c(0,0)
    sigma <- matrix(c(1,1,1,2), nrow=2)
    chol_sigma <- chol(sigma)
    x <- t(mean + t(matrix(stats::rnorm(n*N), nrow=N,ncol=n)%*%chol_sigma))

    original_hts <- shrink.estim(x, diag(diag(crossprod(x) / NROW(x))))
    FoReco_shr <- FoReco::shrink_estim(x, mse = TRUE)
    expect_equal(original_hts[[1]], as.matrix(FoReco_shr))
    expect_equal(original_hts[[2]], attr(FoReco_shr, "lambda"))
  })
}
