# test fun_export.R
if(require(testthat)){
  test_that("Layouts", {
    h <- 10
    agg_order <- 4
    tmp <- tetools(agg_order)
    kt <- tmp$dim["kt"]
    vec <- rnorm(kt*h)

    out <- as_hstack_telayout(vec, agg_order)
    index <- c(0, cumsum(agg_order/tmp$set))
    for(i in 1:length(index[-1])){
      expect_equal(as.vector(t(out[,(index[i]+1):index[i+1]])),
                   vec[(h*index[i]+1):(h*index[i+1])])
    }
    expect_error(as_hstack_telayout(vec, agg_order-1))
    expect_equal(vec, unname(as_tevector(out, agg_order)))

    n <- 3
    mat <- rbind(rnorm(kt*h), rnorm(kt*h), rnorm(kt*h))
    index <- c(0, cumsum(agg_order/tmp$set))
    out <- as_hstack_ctlayout(mat, agg_order) # matrix h x (3*kt)
    for(i in 1:length(index[-1])){
      for(j in 1:n){
        expect_equal(as.vector(t(out[,((index[i]+1):index[i+1])+(kt*(j-1))])),
                     mat[j, (h*index[i]+1):(h*index[i+1])])
      }
    }
    expect_error(as_hstack_ctlayout(mat, agg_order-1))
    expect_equal(mat, unname(as_ctmatrix(out, agg_order, n)))
    expect_equal(mat, unname(as_ctmatrix(out, agg_order, n,
                                         row_names = LETTERS[1:n])))
  })

  test_that("arrange_hres", {
    # Input: 4 (forecast horizons) vectors with 4*10 elements
    input <-  list(rnorm(4*10), rnorm(4*10), rnorm(4*10), rnorm(4*10))
    # Output: 1 vector with 4*10 elements
    out <- arrange_hres(input)
    for(i in 1:4){
      expect_equal(out[rep(1:4, 10)==i], input[[i]][rep(1:4, 10)==i])
    }

    # Matrix version
    input <-  list(matrix(rnorm(4*10*3), 4*10), matrix(rnorm(4*10*3), 4*10),
                   matrix(rnorm(4*10*3), 4*10), matrix(rnorm(4*10*3), 4*10))
    out <- arrange_hres(input)
    for(i in 1:4){
      expect_equal(out[rep(1:4, 10)==i, ], input[[i]][rep(1:4, 10)==i, ])
    }

    expect_error(arrange_hres(input[[1]]))
    expect_equal(arrange_hres(input[1]), input[[1]])
  })
}
