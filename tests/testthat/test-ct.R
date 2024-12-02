# test cross-sectional reconciliation
if(require(testthat)){
  A <- matrix(c(1,1,1,1,
                1,1,0,0,
                0,0,1,1), 3, byrow = TRUE)
  agg_order <- 4
  set.seed(123)
  base <- rbind(rnorm(7, rep(c(40, 20, 10), c(1, 2, 4))),
                rnorm(7, rep(c(20, 10, 5), c(1, 2, 4))),
                rnorm(7, rep(c(20, 10, 5), c(1, 2, 4))),
                rnorm(7, rep(c(10, 5, 2.5), c(1, 2, 4))),
                rnorm(7, rep(c(10, 5, 2.5), c(1, 2, 4))),
                rnorm(7, rep(c(10, 5, 2.5), c(1, 2, 4))),
                rnorm(7, rep(c(10, 5, 2.5), c(1, 2, 4))))
  res <- matrix(rnorm(70*7), nrow = 7)
  Ccs <- cbind(diag(NROW(A)), -A)
  comb <- "bdshr"

  test_that("Cross-temporal tools", {
    expect_equal(cttools(agg_mat = A,
                         agg_order = agg_order,
                         sparse = FALSE)[c("dim", "set")],
                 list(dim = c(n = 7, na = 3, nb = 4,
                              m = 4, p = 3, ks = 3, kt = 7),
                      set = c(4, 2, 1)))
  })

  Cte <- tetools(agg_order, sparse = FALSE)$cons_mat

  test_that("Optimal cross-temporal reconciliation", {
    r1 <- ctrec(base = base, agg_mat = A, agg_order = agg_order, comb = comb,
                res = res, approach = "strc")
    r2 <- ctrec(base = base, agg_mat = A, agg_order = agg_order, comb = comb,
                res = res, approach = "proj")
    r3 <- ctrec(base = base, agg_mat = A, agg_order = agg_order, comb = comb,
                res = res, approach = "strc_osqp")
    r4 <- ctrec(base = base, agg_mat = A, agg_order = agg_order, comb = comb,
                res = res, approach = "proj_osqp")
    r5 <- ctrec(base = base, cons_mat = Ccs, agg_order = agg_order, comb = comb,
                res = res, approach = "strc")
    r6 <- ctrec(base = base, cons_mat = Ccs, agg_order = agg_order, comb = comb,
                res = res, approach = "proj")

    expect_equal(r1, r2, ignore_attr = TRUE)
    expect_equal(r1, r3, ignore_attr = TRUE)
    expect_equal(r1, r4, ignore_attr = TRUE)
    expect_equal(r1, r5, ignore_attr = TRUE)
    expect_equal(r1, r6, ignore_attr = TRUE)
    expect_equal(max(abs(Ccs%*%t(r1))), 0)
    expect_equal(max(abs(Cte%*%r1)), 0)
  })

  test_that("Cross-temporal reconciliation with agg_order subset", {
    r1 <- ctrec(base = base[, -c(2:3)], agg_mat = A, agg_order = c(4, 1), comb = "wlsv",
                res = res[, -c(11:30)])
    r2 <- iterec(base = base[, -c(2:3)], cslist = list(agg_mat = A, comb = "wls"),
                 telist = list(agg_order = c(4, 1), comb = "wlsv"),
                 res = res[, -c(11:30)], verbose = FALSE, tol = 1e-8)


    r3 <- ctrec(base = base[, -c(2:3)], agg_mat = A, agg_order = c(4, 1), comb = "ols",
                res = res[, -c(11:30)])
    r4 <- iterec(base = base[, -c(2:3)], cslist = list(agg_mat = A, comb = "ols"),
                 telist = list(agg_order = c(4, 1), comb = "ols"),
                 res = res[, -c(11:30)], verbose = FALSE)

    r5 <- ctrec(base = base[, -c(2:3)], agg_mat = A, agg_order = c(4, 1), comb = "str",
                res = res[, -c(11:30)])
    r6 <- iterec(base = base[, -c(2:3)], cslist = list(agg_mat = A, comb = "str"),
                 telist = list(agg_order = c(4, 1), comb = "str"),
                 res = res[, -c(11:30)], verbose = FALSE)
    r6k <- tcsrec(base = base[, -c(2:3)], cslist = list(agg_mat = A, comb = "str"),
                  telist = list(agg_order = c(4, 1), comb = "str"),
                  res = res[, -c(11:30)])

    r7 <- ctrec(base = base[, -c(2:3)], agg_mat = A, agg_order = c(4, 1), comb = "csstr",
                res = res[, -c(11:30)])
    r8 <- iterec(base = base[, -c(2:3)], cslist = list(agg_mat = A, comb = "str"),
                 telist = list(agg_order = c(4, 1), comb = "ols"),
                 res = res[, -c(11:30)], verbose = FALSE)
    r8k <- tcsrec(base = base[, -c(2:3)], cslist = list(agg_mat = A, comb = "str"),
                  telist = list(agg_order = c(4, 1), comb = "ols"),
                  res = res[, -c(11:30)])

    r9 <- ctrec(base = base[, -c(2:3)], agg_mat = A, agg_order = c(4, 1), comb = "testr",
                res = res[, -c(11:30)])
    r10 <- iterec(base = base[, -c(2:3)], cslist = list(agg_mat = A, comb = "ols"),
                  telist = list(agg_order = c(4, 1), comb = "str"),
                  res = res[, -c(11:30)], verbose = FALSE)
    r10k <- tcsrec(base = base[, -c(2:3)], cslist = list(agg_mat = A, comb = "ols"),
                   telist = list(agg_order = c(4, 1), comb = "str"),
                   res = res[, -c(11:30)])

    expect_true(max(abs(r1-r2)) < 1e-6)
    expect_equal(r3, r4, ignore_attr = TRUE)
    expect_equal(r5, r6, ignore_attr = TRUE)
    expect_equal(r5, r6k, ignore_attr = TRUE)
    expect_equal(r7, r8, ignore_attr = TRUE)
    expect_equal(r7, r8k, ignore_attr = TRUE)
    expect_equal(r9, r10, ignore_attr = TRUE)
    expect_equal(r9, r10k, ignore_attr = TRUE)
  })

  base[4:5,NCOL(base)] <- -10
  test_that("Optimal nonegative cross-sectional reconciliation", {
    r1 <- ctrec(base = base, agg_mat = A, agg_order = agg_order, comb = comb,
                res = res, approach = "proj", nn = "strc_osqp")
    r2 <- ctrec(base = base, agg_mat = A, agg_order = agg_order, comb = comb,
                res = res, approach = "proj", nn = "proj_osqp")
    r3 <- ctrec(base = base, agg_mat = A, agg_order = agg_order, comb = comb,
                res = res, approach = "proj", nn = "sntz")
    r4 <- ctrec(base = base, agg_mat = A, agg_order = agg_order, comb = comb,
                res = res, approach = "proj", nn = "bpv")

    expect_equal(r1, r2, ignore_attr = TRUE)
    expect_equal(r1, r4, ignore_attr = TRUE)
    expect_equal(max(abs(Ccs%*%t(r1))), 0)
    expect_equal(max(abs(Ccs%*%t(r3))), 0)
    expect_equal(max(abs(Cte%*%r1)), 0)
    expect_equal(max(abs(Cte%*%r3)), 0)
  })

  test_that("Optimal immutable cross-sectional reconciliation", {
    r1 <- ctrec(base = base, agg_mat = A, agg_order = agg_order, comb = comb, res = res,
                approach = "strc", immutable = c(1, 4, 1))
    r2 <- ctrec(base = base, agg_mat = A, agg_order = agg_order, comb = comb, res = res,
                approach = "proj", immutable = c(1, 4, 1))
    r3 <- ctrec(base = base, agg_mat = A, agg_order = agg_order, comb = comb, res = res,
                approach = "strc_osqp", immutable = c(1, 4, 1))
    r4 <- ctrec(base = base, agg_mat = A, agg_order = agg_order, comb = comb, res = res,
                approach = "proj_osqp", immutable = c(1, 4, 1))
    r5 <- ctrec(base = base, agg_mat = A, agg_order = agg_order, comb = comb, res = res,
                approach = "strc_osqp", immutable = c(1, 4, 1), nn = "osqp")
    r6 <- ctrec(base = base, agg_mat = A, agg_order = agg_order, comb = comb, res = res,
                approach = "proj_osqp", immutable = c(1, 4, 1), nn = "osqp")

    fix_r <- c(r1[1,1], r2[1,1], r3[1,1], r4[1,1], r5[1,1], r6[1,1])

    expect_equal(r1, r2, ignore_attr = TRUE)
    expect_equal(r1, r3, ignore_attr = TRUE)
    expect_equal(r1, r4, ignore_attr = TRUE)
    expect_equal(r5, r6, ignore_attr = TRUE)
    expect_equal(max(abs(fix_r - base[1,1])), 0)
    expect_equal(max(abs(Ccs%*%t(r1))), 0)
    expect_equal(max(abs(Cte%*%r1)), 0)
    expect_equal(max(abs(Ccs%*%t(r5))), 0)
    expect_equal(max(abs(Cte%*%r5)), 0)
  })

  test_that("cslcc and BU", {
    r0 <- ctlcc(base = base, agg_mat = A, agg_order = agg_order, comb = comb, res = res)
    r1 <- ctlcc(base = base, agg_mat = A, agg_order = agg_order, nodes = c(1, 2), comb = comb, res = res)
    r2 <- ctbu(base[-c(1:3),-c(1:3)], agg_mat = A, agg_order = agg_order)

    expect_equal(r1, r0, ignore_attr = TRUE)
    expect_equal(recoinfo(r1, verbose = FALSE)$lcc[[9]], r2, ignore_attr = TRUE)
    expect_equal(max(abs(Ccs%*%t(r1))), 0)
    expect_equal(max(abs(Cte%*%r1)), 0)
    expect_equal(max(abs(Ccs%*%t(r2))), 0)
    expect_equal(max(abs(Cte%*%r2)), 0)
  })

  test_that("Top-down and Middle-out", {
    topf <- rnorm(2, 10)

    fix_weights <- matrix(runif(4*4), 4, 4)
    r1 <- cttd(base = topf, agg_mat = A, agg_order = agg_order, weights = fix_weights)

    h_weights <- cbind(fix_weights, fix_weights)
    r2 <- cttd(base = topf, agg_mat = A, agg_order = agg_order, weights = h_weights)

    # Normalization check
    r3 <- cttd(base = topf, agg_mat = A, agg_order = agg_order, weights = fix_weights/sum(fix_weights))

    # Middle-out
    r4 <- ctmo(base = topf, agg_mat = A, agg_order = agg_order, weights = fix_weights)
    r5 <- ctmo(base = topf, agg_mat = A, agg_order = agg_order, weights = h_weights)
    r6 <- ctmo(base = cbind(topf), agg_mat = A, agg_order = agg_order,
               weights = fix_weights, id_rows = 2:3)
    r7 <- ctmo(base = rbind(topf), agg_mat = A, agg_order = agg_order,
               weights = fix_weights, order = 2)

    expect_equal(max(abs(r1[1,1:2] - topf)), 0)

    expect_equal(max(abs(r7[1,2:3] - topf)), 0)
    expect_equal(max(abs(r6[2:3, 1] - topf)), 0)
    expect_equal(r1, r2, ignore_attr = TRUE)
    expect_equal(r1, r3, ignore_attr = TRUE)
    expect_equal(r1, r4, ignore_attr = TRUE)
    expect_equal(r1, r5, ignore_attr = TRUE)
    expect_equal(max(abs(Ccs%*%t(r6))), 0)
    expect_equal(max(abs(Ccs%*%t(r7))), 0)
    expect_equal(max(abs(Cte%*%r1)), 0)
    expect_equal(max(abs(Cte%*%r6)), 0)
    expect_equal(max(abs(Cte%*%r7)), 0)
  })

  test_that("Cross-sectional tools", {
    M <- ctprojmat(cons_mat = Ccs, agg_order = agg_order, comb = comb, res = res)
    G <- ctprojmat(agg_mat = A, agg_order = agg_order, comb = comb, res = res, mat = "G")
    S <- cttools(agg_mat = A, agg_order = agg_order)$strc_mat

    expect_equal(M, S%*%G, ignore_attr = TRUE)
  })

  test_that("Covariance", {
    for(i in c("acov", "bdsam", "bdshr", "bsam", "bshr", "csstr", "hbsam", "hbshr",
               "hsam", "hshr", "ols", "sam", "shr", "Ssam", "Sshr", "str", "testr", "wlsh",
               "wlsv")){
      expect_no_error(ctrec(base = base, agg_mat = A, agg_order = agg_order, comb = i,
                            res = res))
    }
  })
}
