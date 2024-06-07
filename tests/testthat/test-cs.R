# test cross-sectional reconciliation
if(require(testthat)){
  A <- matrix(c(1,1,1,1,
                1,1,0,0,
                0,0,1,1), 3, byrow = TRUE)
  set.seed(123)
  res <- matrix(rnorm(100*sum(dim(A))), 100, sum(dim(A)))
  base <- t(rnorm(sum(dim(A)), 1))
  C <- cbind(diag(NROW(A)), -A)
  comb <- "shr"

  test_that("Cross-sectional tools", {
    expect_equal(cstools(agg_mat = A,
                         sparse = FALSE),
                 list(dim = c(n = 7, na = 3, nb = 4),
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

  test_that("Optimal cross-sectional reconciliation", {
    r1 <- csrec(base = base, agg_mat = A, comb = comb,
                res = res, approach = "strc")
    r2 <- csrec(base = base, agg_mat = A, comb = comb,
                res = res, approach = "proj")
    r3 <- csrec(base = base, agg_mat = A, comb = comb,
                res = res, approach = "strc_osqp")
    r4 <- csrec(base = base, agg_mat = A, comb = comb,
                res = res, approach = "proj_osqp")
    r5 <- csrec(base = base, cons_mat = C, comb = comb,
                res = res, approach = "strc")
    r6 <- csrec(base = base, cons_mat = C, comb = comb,
                res = res, approach = "proj")

    expect_equal(r1, r2, ignore_attr = TRUE)
    expect_equal(r1, r3, ignore_attr = TRUE)
    expect_equal(r1, r4, ignore_attr = TRUE)
    expect_equal(r1, r5, ignore_attr = TRUE)
    expect_equal(r1, r6, ignore_attr = TRUE)
    expect_equal(max(abs(C%*%t(r1))), 0)
  })

  base[1,NCOL(base)] <- -10
  test_that("Optimal nonegative cross-sectional reconciliation", {
    r1 <- csrec(base = base, agg_mat = A, comb = comb, res = res,
                approach = "proj", nn = "strc_osqp")
    r2 <- csrec(base = base, agg_mat = A, comb = comb, res = res,
                approach = "proj", nn = "proj_osqp")
    r3 <- csrec(base = base, agg_mat = A, comb = comb, res = res,
                approach = "proj", nn = "sntz")


    expect_equal(r1, r2, ignore_attr = TRUE)
    expect_equal(max(abs(C%*%t(r1))), 0)
    expect_equal(max(abs(C%*%t(r3))), 0)
  })

  test_that("Optimal immutable cross-sectional reconciliation", {
    r1 <- csrec(base = base, agg_mat = A, comb = comb, res = res,
                approach = "strc", immutable = 1)
    r2 <- csrec(base = base, agg_mat = A, comb = comb, res = res,
                approach = "proj", immutable = 1)
    r3 <- csrec(base = base, agg_mat = A, comb = comb, res = res,
                approach = "strc_osqp", immutable = 1)
    r4 <- csrec(base = base, agg_mat = A, comb = comb, res = res,
                approach = "proj_osqp", immutable = 1)
    r5 <- csrec(base = base, agg_mat = A, comb = comb, res = res,
                approach = "strc_osqp", immutable = 1, nn = "osqp")
    r6 <- csrec(base = base, agg_mat = A, comb = comb, res = res,
                approach = "proj_osqp", immutable = 1, nn = "osqp")

    fix_r <- c(r1[1,1], r2[1,1], r3[1,1], r4[1,1], r5[1,1], r6[1,1])

    expect_equal(r1, r2, ignore_attr = TRUE)
    expect_equal(r1, r3, ignore_attr = TRUE)
    expect_equal(r1, r4, ignore_attr = TRUE)
    expect_equal(r5, r6, ignore_attr = TRUE)
    expect_equal(max(abs(fix_r - base[1,1])), 0)
    expect_equal(max(abs(C%*%t(r1))), 0)
    expect_equal(max(abs(C%*%t(r5))), 0)
  })

  test_that("cslcc and BU", {
    r0 <- cslcc(base = base, agg_mat = A, comb = comb, res = res)
    r1 <- cslcc(base = base, agg_mat = A,  nodes = c(1, 2), comb = comb, res = res)
    r2 <- csbu(base[,-c(1:3)], agg_mat = A)

    fix <- unlist(mapply(function(z, y) z[,y], y = list(1, c(2,3), c(4:7)),
                         z = recoinfo(r1, verbose = FALSE)$lcc))

    expect_equal(max(abs(fix - base)), 0)
    expect_equal(r1, r0, ignore_attr = TRUE)
    expect_equal(recoinfo(r1, verbose = FALSE)$lcc[[3]], r2)
    expect_equal(max(abs(C%*%t(r1))), 0)
    expect_equal(max(abs(C%*%t(r2))), 0)
  })

  test_that("Top-down and Middle-out", {
    topf <- rnorm(2, 10)
    fix_weights <- runif(4)
    r1 <- cstd(base = topf, agg_mat = A, weights = fix_weights)

    h_weights <- rbind(fix_weights, fix_weights)
    r2 <- cstd(base = topf, agg_mat = A, weights = h_weights)

    # Normalization check
    r3 <- cstd(base = topf, agg_mat = A, weights = fix_weights/sum(fix_weights))

    # Middle-out
    r4 <- csmo(base = cbind(topf), agg_mat = A, weights = fix_weights)
    r5 <- csmo(base = cbind(topf), agg_mat = A, weights = h_weights)
    r6 <- csmo(base = rbind(topf), agg_mat = A, weights = fix_weights, id_rows = 2:3)

    expect_equal(max(abs(r1[,1] - topf)), 0)
    expect_equal(max(abs(r6[,2:3] - topf)), 0)
    expect_equal(max(abs(r6[,1] - sum(topf))), 0)
    expect_equal(r1, r2, ignore_attr = TRUE)
    expect_equal(r1, r3, ignore_attr = TRUE)
    expect_equal(r1, r4, ignore_attr = TRUE)
    expect_equal(r1, r5, ignore_attr = TRUE)
    expect_equal(max(abs(C%*%t(r1))), 0)
    expect_equal(max(abs(C%*%t(r6))), 0)
  })

  test_that("Cross-sectional tools", {
    M <- csprojmat(cons_mat = C, comb = "shr", res = res)
    G <- csprojmat(agg_mat = A, comb = "shr", res = res, mat = "G")
    S <- cstools(agg_mat = A)$strc_mat

    expect_equal(M, S%*%G, ignore_attr = TRUE)
  })
}
