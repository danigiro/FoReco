# test cross-sectional reconciliation
if (require(testthat)) {
  A <- matrix(c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1), 3, byrow = TRUE)
  colnames(A) <- LETTERS[4:7]
  rownames(A) <- LETTERS[1:3]
  set.seed(123)
  res <- matrix(rnorm(100 * sum(dim(A))), 100, sum(dim(A)))
  base <- t(rnorm(sum(dim(A)), rep(c(40, 20, 10), c(1, 2, 4))))
  C <- (cbind(diag(NROW(A)), -A))
  colnames(C) <- LETTERS[1:7]
  comb <- "shr"
  test_that("Cross-sectional tools", {
    expect_equal(
      cstools(agg_mat = A, sparse = FALSE),
      #fmt: skip
      list(dim = c(n = 7, na = 3, nb = 4),
           agg_mat = matrix(c(1,1,1,1,
                              1,1,0,0,
                              0,0,1,1), 3, 4, byrow = TRUE,
                            dimnames = list(LETTERS[1:3], LETTERS[4:7])),
           strc_mat = matrix(c(1,1,1,1,
                               1,1,0,0,
                               0,0,1,1,
                               1,0,0,0,
                               0,1,0,0,
                               0,0,1,0,
                               0,0,0,1), 7, 4, byrow = TRUE,
                             dimnames = list(LETTERS[1:7], LETTERS[4:7])),
           cons_mat = matrix(c(1,0,0,-1,-1,-1,-1,
                               0,1,0,-1,-1,0,0,
                               0,0,1,0,0,-1,-1), 3, 7, byrow = TRUE,
                             dimnames = list(LETTERS[1:3], LETTERS[1:7])))
    )
    expect_equal(
      cstools(cons_mat = C[-2, ], sparse = FALSE),
      list(
        dim = c(n = 7),
        cons_mat = matrix(
          c(1, 0, 0, -1, -1, -1, -1, 0, 0, 1, 0, 0, -1, -1),
          2,
          7,
          byrow = TRUE,
          dimnames = list(LETTERS[c(1, 3)], LETTERS[1:7])
        )
      )
    )
  })

  test_that("Optimal cross-sectional reconciliation", {
    r1 <- csrec(
      base = base,
      agg_mat = A,
      comb = comb,
      res = res,
      approach = "strc"
    )
    r2 <- csrec(
      base = base,
      agg_mat = A,
      comb = comb,
      res = res,
      approach = "proj"
    )
    r3 <- csrec(
      base = base,
      agg_mat = A,
      comb = comb,
      res = res,
      approach = "strc_osqp"
    )
    r4 <- csrec(
      base = base,
      agg_mat = A,
      comb = comb,
      res = res,
      approach = "proj_osqp"
    )
    r5 <- csrec(
      base = base,
      cons_mat = C,
      comb = comb,
      res = res,
      approach = "strc"
    )
    r6 <- csrec(
      base = base,
      cons_mat = C,
      comb = comb,
      res = res,
      approach = "proj"
    )

    r7 <- csrec(base = base[1, ], agg_mat = A, comb = comb, res = res)

    expect_equal(r1, r2, ignore_attr = TRUE)
    expect_equal(r1, r3, ignore_attr = TRUE)
    expect_equal(r1, r4, ignore_attr = TRUE)
    expect_equal(r1, r5, ignore_attr = TRUE)
    expect_equal(r1, r6, ignore_attr = TRUE)
    expect_equal(r1, r7, ignore_attr = TRUE)
    expect_equal(FoReco2matrix(r1)$`k-1`, r1, ignore_attr = TRUE)
    expect_equal(max(abs(C %*% t(r1))), 0)
  })

  test_that("Bounds", {
    tmp <- set_bounds(n = c(4:7), lb = 10)
    r1 <- csrec(
      base = base,
      agg_mat = A,
      comb = comb,
      res = res,
      approach = "strc",
      bounds = tmp
    )
    r2 <- csrec(
      base = base,
      agg_mat = A,
      comb = comb,
      res = res,
      approach = "proj",
      bounds = tmp
    )
    expect_equal(r1[1, ], c(40, 20, 20, 10, 10, 10, 10), ignore_attr = TRUE)
    expect_equal(r1, r2, ignore_attr = TRUE)

    tmp <- set_bounds(n = c(4:7), lb = 10, approach = "sftb")
    r3 <- csrec(
      base = base,
      agg_mat = A,
      comb = comb,
      res = res,
      approach = "strc",
      bounds = tmp
    )
    expect_true(all(r3 >= 10))

    tmp <- set_bounds(n = 8, lb = 3)
    expect_warning(csrec(
      base = base,
      agg_mat = A,
      comb = comb,
      res = res,
      bounds = tmp
    ))
  })

  base[1, NCOL(base)] <- -10
  test_that("Optimal nonegative cross-sectional reconciliation", {
    r1 <- csrec(
      base = base,
      agg_mat = A,
      comb = comb,
      res = res,
      approach = "proj",
      nn = "strc_osqp"
    )
    r2 <- csrec(
      base = base,
      agg_mat = A,
      comb = comb,
      res = res,
      approach = "proj",
      nn = "proj_osqp"
    )
    r3 <- csrec(
      base = base,
      agg_mat = A,
      comb = comb,
      res = res,
      approach = "proj",
      nn = "sntz"
    )
    r4 <- csrec(
      base = base,
      agg_mat = A,
      comb = comb,
      res = res,
      approach = "proj",
      nn = "bpv"
    )
    r5 <- csrec(
      base = base,
      agg_mat = A,
      comb = comb,
      res = res,
      approach = "proj",
      nn = "nnic"
    )
    r6 <- csrec(
      base = base,
      agg_mat = A,
      comb = comb,
      res = res,
      approach = "proj",
      nn = "nfca"
    )
    r7 <- csrec(
      base = base,
      agg_mat = A,
      comb = comb,
      res = res,
      approach = "proj",
      nn = "sntz",
      settings = list(type = "tdp")
    )
    r8 <- csrec(
      base = base,
      agg_mat = A,
      comb = comb,
      res = res,
      approach = "proj",
      nn = "sntz",
      settings = list(type = "tdsp")
    )
    r9 <- csrec(
      base = base,
      agg_mat = A,
      comb = comb,
      res = res,
      approach = "proj",
      nn = "sntz",
      settings = list(type = "tdvw")
    )

    rf <- csrec(
      base = base,
      agg_mat = A,
      comb = comb,
      res = res,
      approach = "proj"
    )
    rbu <- csbu(rf[, -c(1:3)], agg_mat = unname(A), sntz = TRUE)

    expect_true(all(r1 >= 0))
    expect_true(all(r3 >= 0))
    expect_true(all(r6 >= 0))
    expect_true(all(r7 >= 0))
    expect_true(all(r8 >= 0))
    expect_true(all(r9 >= 0))
    expect_equal(r1, r2, ignore_attr = TRUE)
    expect_equal(r1, r4, ignore_attr = TRUE)
    expect_equal(r1, r5, ignore_attr = TRUE)
    expect_equal(r3, rbu, ignore_attr = TRUE)
    expect_equal(r7[1], rf[1])
    expect_equal(r8[1], rf[1])
    expect_equal(r9[1], rf[1])
    expect_equal(max(abs(C %*% t(r1))), 0)
    expect_equal(max(abs(C %*% t(r3))), 0)
    expect_equal(max(abs(C %*% t(r6))), 0)
    expect_equal(max(abs(C %*% t(r7))), 0)
    expect_equal(max(abs(C %*% t(r8))), 0)
    expect_equal(max(abs(C %*% t(r9))), 0)
  })

  test_that("sntz-variants", {
    yhat <- c(40, 35, -5, 10)
    agg_mat <- t(c(1, 1, 1))
    r1 <- csrec(yhat, agg_mat = agg_mat, nn = "sntz")
    r2 <- csrec(
      yhat,
      agg_mat = agg_mat,
      nn = "sntz",
      settings = list(type = "bu")
    )
    r3 <- csrec(
      yhat,
      agg_mat = agg_mat,
      nn = "sntz",
      settings = list(type = "tdp")
    )
    r4 <- csrec(
      yhat,
      agg_mat = agg_mat,
      nn = "sntz",
      settings = list(type = "tdsp")
    )
    r5 <- csrec(
      yhat,
      agg_mat = agg_mat,
      nn = "sntz",
      settings = list(type = "tdvw"),
      comb = diag(c(1, 64, 1, 16))
    )
    expect_equal(r1, c(45, 35, 0, 10), ignore_attr = TRUE)
    expect_equal(r2, c(45, 35, 0, 10), ignore_attr = TRUE)
    expect_equal(round(r3, 2), c(40, 31.11, 0, 8.89), ignore_attr = TRUE)
    expect_equal(round(r4, 2), c(40, 30.38, 0, 9.62), ignore_attr = TRUE)
    expect_equal(r5, c(40, 31, 0, 9), ignore_attr = TRUE)
  })

  test_that("Optimal immutable cross-sectional reconciliation", {
    r1 <- csrec(
      base = base,
      agg_mat = A,
      comb = comb,
      res = res,
      approach = "strc",
      immutable = 1
    )
    r2 <- csrec(
      base = base,
      agg_mat = A,
      comb = comb,
      res = res,
      approach = "proj",
      immutable = 1
    )
    r3 <- csrec(
      base = base,
      agg_mat = A,
      comb = comb,
      res = res,
      approach = "strc_osqp",
      immutable = 1
    )
    r4 <- csrec(
      base = base,
      agg_mat = A,
      comb = comb,
      res = res,
      approach = "proj_osqp",
      immutable = 1
    )
    r5 <- csrec(
      base = base,
      agg_mat = A,
      comb = comb,
      res = res,
      approach = "strc_osqp",
      immutable = 1,
      nn = "osqp"
    )
    r6 <- csrec(
      base = base,
      agg_mat = A,
      comb = comb,
      res = res,
      approach = "proj_osqp",
      immutable = 1,
      nn = "osqp"
    )
    expect_no_error(suppressMessages(recoinfo(r1, verbose = TRUE)))

    fix_r <- c(r1[1, 1], r2[1, 1], r3[1, 1], r4[1, 1], r5[1, 1], r6[1, 1])

    expect_equal(r1, r2, ignore_attr = TRUE)
    expect_equal(r1, r3, ignore_attr = TRUE)
    expect_equal(r1, r4, ignore_attr = TRUE)
    expect_equal(r5, r6, ignore_attr = TRUE)
    expect_equal(max(abs(fix_r - base[1, 1])), 0)
    expect_equal(max(abs(C %*% t(r1))), 0)
    expect_equal(max(abs(C %*% t(r5))), 0)
  })

  test_that("cslcc and BU", {
    r0 <- cslcc(base = base, agg_mat = A, comb = comb, res = res)
    r1 <- cslcc(
      base = base,
      agg_mat = A,
      nodes = c(1, 2),
      comb = comb,
      res = res
    )
    r2 <- csbu(base[, -c(1:3)], agg_mat = A)
    r3 <- csbu(base[, -c(1:3)], agg_mat = A, round = TRUE)

    expect_error(csbu(base, agg_mat = A))
    expect_error(csbu(base))
    expect_error(csbu(agg_mat = A))

    fix <- unlist(mapply(
      function(z, y) z[, y],
      y = list(1, c(2, 3), c(4:7)),
      z = recoinfo(r1, verbose = FALSE)$lcc
    ))

    expect_no_error(suppressMessages(recoinfo(r1, verbose = TRUE)))
    expect_no_error(suppressMessages(recoinfo(r2, verbose = TRUE)))
    expect_equal(max(abs(fix - base)), 0)
    expect_equal(r1, r0, ignore_attr = TRUE)
    expect_equal(recoinfo(r1, verbose = FALSE)$lcc[[3]], r2)
    expect_equal(max(abs(C %*% t(r1))), 0)
    expect_equal(max(abs(C %*% t(r2))), 0)
    expect_equal(max(abs(C %*% t(r3))), 0)
    expect_equal(r3, round(r3), ignore_attr = TRUE)
  })

  test_that("Top-down and Middle-out", {
    topf <- rnorm(2, 10)
    fix_weights <- runif(4)
    r1 <- cstd(base = topf, agg_mat = A, weights = fix_weights)

    h_weights <- rbind(fix_weights, fix_weights)
    r2 <- cstd(base = topf, agg_mat = A, weights = h_weights)

    # Normalization check
    r3 <- cstd(
      base = topf,
      agg_mat = A,
      weights = fix_weights / sum(fix_weights)
    )

    # Middle-out
    r4 <- csmo(base = cbind(topf), agg_mat = A, weights = fix_weights)
    r5 <- csmo(base = cbind(topf), agg_mat = A, weights = h_weights)
    r6 <- csmo(
      base = rbind(topf),
      agg_mat = A,
      weights = fix_weights,
      id_rows = 2:3
    )
    expect_no_error(suppressMessages(recoinfo(r1, verbose = TRUE)))
    expect_no_error(suppressMessages(recoinfo(r4, verbose = TRUE)))

    expect_equal(max(abs(r1[, 1] - topf)), 0)
    expect_equal(max(abs(r6[, 2:3] - topf)), 0)
    expect_equal(max(abs(r6[, 1] - sum(topf))), 0)
    expect_equal(r1, r2, ignore_attr = TRUE)
    expect_equal(r1, r3, ignore_attr = TRUE)
    expect_equal(r1, r4, ignore_attr = TRUE)
    expect_equal(r1, r5, ignore_attr = TRUE)
    expect_equal(max(abs(C %*% t(r1))), 0)
    expect_equal(max(abs(C %*% t(r6))), 0)
  })

  test_that("Cross-sectional tools", {
    M <- csprojmat(cons_mat = C, comb = "shr", res = res)
    G <- csprojmat(agg_mat = A, comb = "shr", res = res, mat = "G")
    S <- cstools(agg_mat = A)$strc_mat

    expect_equal(M, unname(S %*% G), ignore_attr = TRUE)
  })

  test_that("Covariance", {
    for (i in c("wls", "shr", "sam", "oasd")) {
      expect_no_error(csrec(base = base, agg_mat = A, comb = i, res = res))
    }

    for (i in c("ols", "str")) {
      expect_no_error(csrec(base = base, agg_mat = A, comb = i))
    }
  })

  test_that("Errors", {
    expect_error(csrec(base = base, comb = comb, res = res))
    expect_error(csrec(agg_mat = A, comb = comb, res = res))
    expect_error(csrec(base = base[, 1:2], agg_mat = A, comb = comb, res = res))
    expect_error(cstools())
    expect_error(csrec(
      base = base,
      agg_mat = A,
      comb = comb,
      res = res,
      immutable = cbind(1, 1)
    ))
    expect_error(csrec(
      base = base,
      agg_mat = A,
      comb = comb,
      res = res,
      immutable = c(1:7)
    ))
  })
}
