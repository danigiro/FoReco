# test temporal reconciliation
if(require(testthat)){
  agg_order <- 4
  set <- c(4, 2, 1)
  set.seed(123)
  res <- rnorm(100*sum(set))
  base <- rnorm(sum(set), 1)

  test_that("Temporal tools", {
    expect_equal(tetools(agg_order = 4, sparse = FALSE),
                 list(dim = c(m = 4, p = 3, ks = 3, kt = 7),
                      set = c(4, 2, 1),
                      agg_mat = matrix(c(1,1,1,1,
                                         1,1,0,0,
                                         0,0,1,1), 3, 4, byrow = TRUE),
                      strc_mat = matrix(c(1,1,1,1,
                                          1,1,0,0,
                                          0,0,1,1,
                                          1,0,0,0,
                                          0,1,0,0,
                                          0,0,1,0,
                                          0,0,0,1), 7, 4, byrow = TRUE),
                      cons_mat = matrix(c(1,0,0,-1,-1,-1,-1,
                                          0,1,0,-1,-1,0,0,
                                          0,0,1,0,0,-1,-1), 3, 7, byrow = TRUE)))
  })

  C <- tetools(agg_order, sparse = FALSE)$cons_mat
  comb <- "wlsv"

  test_that("Optimal temporal reconciliation", {
    r1 <- terec(base = base, agg_order = agg_order, comb = comb,
                res = res, approach = "strc")
    r2 <- terec(base = base, agg_order = agg_order, comb = comb,
                res = res, approach = "proj")
    r3 <- terec(base = base, agg_order = agg_order, comb = comb,
                res = res, approach = "strc_osqp")
    r4 <- terec(base = base, agg_order = agg_order, comb = comb,
                res = res, approach = "proj_osqp")

    expect_equal(r1, r2, ignore_attr = TRUE)
    expect_equal(r1, r3, ignore_attr = TRUE)
    expect_equal(r1, r4, ignore_attr = TRUE)
    expect_equal(max(abs(C%*%r1)), 0)
  })

  base[length(base)] <- -10
  test_that("Optimal nonegative temporal reconciliation", {
    r1 <- terec(base = base, agg_order = agg_order, comb = comb, res = res,
                approach = "proj", nn = "strc_osqp")
    r2 <- terec(base = base, agg_order = agg_order, comb = comb, res = res,
                approach = "proj", nn = "proj_osqp")
    r3 <- terec(base = base, agg_order = agg_order, comb = comb, res = res,
                approach = "proj", nn = "sntz")


    expect_equal(r1, r2, ignore_attr = TRUE)
    expect_equal(max(abs(C%*%r1)), 0)
    expect_equal(max(abs(C%*%r3)), 0)
  })

  test_that("Optimal immutable temporal reconciliation", {
    r1 <- terec(base = base, agg_order = agg_order, comb = comb, res = res,
                approach = "strc", immutable = c(agg_order, 1))
    r2 <- terec(base = base, agg_order = agg_order, comb = comb, res = res,
                approach = "proj", immutable = c(agg_order, 1))
    r3 <- terec(base = base, agg_order = agg_order, comb = comb, res = res,
                approach = "strc_osqp", immutable = c(agg_order, 1))
    r4 <- terec(base = base, agg_order = agg_order, comb = comb, res = res,
                approach = "proj_osqp", immutable = c(agg_order, 1))
    r5 <- terec(base = base, agg_order = agg_order, comb = comb, res = res,
                approach = "strc_osqp", immutable = c(agg_order, 1), nn = "osqp")
    r6 <- terec(base = base, agg_order = agg_order, comb = comb, res = res,
                approach = "proj_osqp", immutable = c(agg_order, 1), nn = "osqp")

    fix_r <- c(r1[1], r2[1], r3[1], r4[1], r5[1], r6[1])

    expect_equal(r1, r2, ignore_attr = TRUE)
    expect_equal(r1, r3, ignore_attr = TRUE)
    expect_equal(r1, r4, ignore_attr = TRUE)
    expect_equal(r5, r6, ignore_attr = TRUE)
    expect_equal(max(abs(fix_r - base[1])), 0)
    expect_equal(max(abs(C%*%r1)), 0)
    expect_equal(max(abs(C%*%r5)), 0)
  })

  test_that("telcc and BU", {
    r1 <- telcc(base = base, agg_order = agg_order, comb = comb, res = res)
    r2 <- tebu(base[-c(1:3)], agg_order = agg_order)

    fix <- unlist(mapply(function(z, y) z[y], y = list(1, c(2,3), c(4:7)),
                         z = recoinfo(r1, verbose = FALSE)$lcc))

    expect_equal(max(abs(fix - base)), 0)
    expect_equal(recoinfo(r1, verbose = FALSE)$lcc[[3]], r2)
    expect_equal(max(abs(C%*%r1)), 0)
    expect_equal(max(abs(C%*%r2)), 0)
  })

  test_that("Top-down and Middle-out", {
    set.seed(123)
    topf <- rnorm(2, 10)

    fix_weights <- runif(4)
    r1 <- tetd(base = topf, agg_order = 4, weights = fix_weights)

    h_weights <- c(fix_weights, fix_weights)
    r2 <- tetd(base = topf, agg_order = 4, weights = h_weights)

    # Normalization check
    r3 <- tetd(base = topf, agg_order = 4, weights = fix_weights/sum(fix_weights))

    # Middle-out
    r4 <- temo(base = topf, agg_order = 4, weights = fix_weights)
    r5 <- temo(base = topf, agg_order = 4, weights = h_weights)
    r6 <- temo(base = topf, order = 2, agg_order = 4, weights = fix_weights)

    expect_equal(max(abs(r1[1:2] - topf)), 0)
    expect_equal(max(abs(r6[2:3] - topf)), 0)
    expect_equal(max(abs(r6[1] - sum(topf))), 0)
    expect_equal(r1, r2, ignore_attr = TRUE)
    expect_equal(r1, r3, ignore_attr = TRUE)
    expect_equal(r1, r4, ignore_attr = TRUE)
    expect_equal(r1, r5, ignore_attr = TRUE)
    expect_equal(max(abs(C%*%r6)), 0)
  })

  test_that("Temporal tools", {
    M <- teprojmat(agg_order = 4, comb = "wlsv", res = res)
    G <- teprojmat(agg_order = 4, comb = "wlsv", res = res, mat = "G")
    S <- tetools(agg_order = 4)$strc_mat

    expect_equal(M, S%*%G, ignore_attr = TRUE)
  })
}
