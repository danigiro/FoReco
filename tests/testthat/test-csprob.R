# test cross-sectional reconciliation
if(require(testthat)){
  A <- matrix(c(1,1,1,1,
                1,1,0,0,
                0,0,1,1), 3, byrow = TRUE)
  colnames(A) <- LETTERS[4:7]
  rownames(A) <- LETTERS[1:3]
  set.seed(123)
  res <- matrix(rnorm(100*sum(dim(A))), 100, sum(dim(A)))
  base <- t(rnorm(sum(dim(A)), 1))
  C <- (cbind(diag(NROW(A)), -A))
  colnames(C) <- LETTERS[1:7]
  comb <- "sam"
  test_that("Cross-sectional gaussian, full cov", {
    rp <- csrec(base = base, agg_mat = A, comb = comb, res = res)
    rd1 <- csgauss(base = base, agg_mat = A, comb = comb, res = res)
    rd2 <- csgauss(base = base, agg_mat = A, comb = comb,
                   comb_base = comb, res = res)
    rd3 <- csgauss(base = base, agg_mat = A, comb = comb,
                   comb_base = "ols", res = res)
    rd4 <- csgauss(base = base, agg_mat = A, comb = comb,
                   res = res, reduce_form = TRUE)

    rd5 <- csgauss(base = base, agg_mat = A, comb = comb,
                   res = res, approach = "strc")
    rd6 <- csgauss(base = base, agg_mat = A, comb = comb, approach = "strc",
                   res = res, reduce_form = TRUE)

    expect_equal(rp, mean(rd1), ignore_attr = TRUE)
    expect_equal(rp, mean(rd3), ignore_attr = TRUE)
    expect_equal(rp, csbu(mean(rd4), A), ignore_attr = TRUE)
    expect_equal(rp, csbu(mean(rd6), A), ignore_attr = TRUE)
    expect_equal(rd1, rd2, ignore_attr = TRUE)
    expect_equal(rd1, rd5, ignore_attr = TRUE)
    expect_equal(rd4, rd6, ignore_attr = TRUE)
  })

  comb <- "wls"
  test_that("Cross-sectional gaussian, diag cov", {
    rp <- csrec(base = base, agg_mat = A, comb = comb, res = res)
    rd1 <- csgauss(base = base, agg_mat = A, comb = comb, res = res)
    rd2 <- csgauss(base = base, agg_mat = A, comb = comb,
                   comb_base = comb, res = res)
    rd3 <- csgauss(base = base, agg_mat = A, comb = comb,
                   comb_base = "ols", res = res)
    rd4 <- csgauss(base = base, agg_mat = A, comb = comb,
                   res = res, reduce_form = TRUE)

    rd5 <- csgauss(base = base, agg_mat = A, comb = comb,
                   res = res, approach = "strc")
    rd6 <- csgauss(base = base, agg_mat = A, comb = comb, approach = "strc",
                   res = res, reduce_form = TRUE)

    expect_equal(rp, mean(rd1), ignore_attr = TRUE)
    expect_equal(rp, mean(rd3), ignore_attr = TRUE)
    expect_equal(rp, csbu(mean(rd4), A), ignore_attr = TRUE)
    expect_equal(rp, csbu(mean(rd6), A), ignore_attr = TRUE)
    expect_equal(rd1, rd2, ignore_attr = TRUE)
    expect_equal(rd1, rd5, ignore_attr = TRUE)
    expect_equal(rd4, rd6, ignore_attr = TRUE)
  })

  set.seed(123)
  A <- t(c(1,1)) # Aggregation matrix for Z = X + Y

  # (100 x 3) base forecasts sample (simulated) for h = 1
  base_h1 <- matrix(rnorm(100*3, mean = c(20, 10, 10)), 100, byrow = TRUE)
  # (100 x 3) base forecasts sample (simulated) for h = 2
  base_h2 <- matrix(rnorm(100*3, mean = c(20, 10, 10)), 100, byrow = TRUE)
  # (2 x 3 x 100) base forecasts sample array with
  # 2 forecast horizons, 3 time series and 100 sample
  base_sample <- aperm(simplify2array(list(base_h1, base_h2)), c(3,2,1))
  res <- t(matrix(rnorm(n = 300), nrow = 3))
  comb <- "shr"
  test_that("Cross-sectional sample", {
    # Optimal cross-sectional probabilistic reconciliation
    r1 <- cssample(base_sample, agg_mat = A, comb = comb, res = res)
    r11 <- apply(base_sample, 3, csrec, agg_mat = A, comb = comb, res = res)
    expect_equal(rowMeans(r11), as.vector(mean(r1)))

    # Bottom-up probabilistic reconciliation
    r2 <- cssample(base_sample[,-1,], agg_mat = A, fun = csbu)
    r21 <- apply(base_sample[,-1,], 3, csbu, agg_mat = A)
    expect_equal(rowMeans(r21), as.vector(mean(r2)))

    # Level conditional coherent probabilistic reconciliation
    r3 <- cssample(base_sample, agg_mat = A, fun = cslcc, comb = comb,
                   res = res)
    r31 <- apply(base_sample, 3, cslcc, agg_mat = A, comb = comb, res = res)
    expect_equal(rowMeans(r31), as.vector(mean(r3)))
  })

}
