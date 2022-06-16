#' @title Level conditional coherent forecast reconciliation for genuine hierarchical/grouped time series
#'
#' @description
#' \loadmathjax
#' Forecast reconciliation procedure built on and extending the original
#' proposal by Hollyman et al. (2021). Level conditional coherent
#' reconciled forecasts may be computed in cross-sectional, temporal, and
#' cross-temporal frameworks. The reconciled forecasts are conditional to (i.e.,
#' constrained by) the base forecasts of a specific upper level of the hierarchy
#' (exogenous constraints). The linear constraints linking the variables may be
#' dealt with endogenously as well (Di Fonzo and Girolimetto, 2021).
#' \emph{Combined Conditional Coherent} (CCC)
#' forecasts may be calculated as simple averages of LCC and bottom-up
#' reconciled forecasts, with either endogenous or exogenous constraints.
#'
#' @usage lccrec(basef, m, C, nl, weights, bnaive = NULL, const = "exogenous",
#'        CCC = TRUE, nn = FALSE, nn_type = "osqp", settings = list(), ...)
#'
#' @param basef matrix/vector of base forecasts to be reconciled:
#' (\mjseqn{h \times n}) matrix in the cross-sectional framework;
#' (\mjseqn{h(k^\ast + m) \times 1}) vector in the temporal framework;
#' (\mjseqn{n \times h(k^\ast+m)}) matrix in the cross-temporal framework.
#' \mjseqn{n} is the total number of variables, \mjseqn{m} is the highest time
#' frequency, \mjseqn{k^\ast} is the sum of (a subset of) (\mjseqn{p-1}) factors
#' of \mjseqn{m}, excluding \mjseqn{m}, and \mjseqn{h} is the forecast horizon.
#' @param C (\mjseqn{n_a \times n_b}) cross-sectional (contemporaneous) matrix
#' mapping the bottom level series into the higher level ones (or a list
#' of matrices forming \mjseqn{\mathbf{C} = [\mathbf{C}_1' \; \mathbf{C}_2' \; ... \; \mathbf{C}_L']'},
#' \mjseqn{1, ..., L} being the number of the cross-sectional upper levels.
#' @param nl (\mjseqn{L \times 1}) vector containing the number of time series
#' in each level of the hierarchy (\code{nl[1] = 1}).
#' @param m Highest available sampling frequency per seasonal cycle (max. order
#' of temporal aggregation, \mjseqn{m}), or a subset of the \mjseqn{p} factors
#' of \mjseqn{m}.
#' @param weights covariance matrix or a vector (weights used in the reconciliation:
#' either (\mjseqn{n_b \times 1}) for exogenous or (\mjseqn{n \times 1}) for
#' endogenous constraints).
#' @param bnaive matrix/vector of naive bts base forecasts
#' (e.g., seasonal averages, as in Hollyman et al., 2021):
#' (\mjseqn{h \times n_b}) matrix in the cross-sectional framework;
#' (\mjseqn{h m \times 1}) vector in the temporal framework;
#' (\mjseqn{n_b \times h m}) matrix in the cross-temporal framework.
#' Ignore it, if only basef are to be used as base forecasts.
#' @param const \strong{exo}genous (\emph{default}) or \strong{endo}genous
#' constraints
#' @param CCC Option to return Combined Conditional Coherent reconciled
#' forecasts (\emph{default} is TRUE).
#' @param nn Logical value: \code{TRUE} if non-negative reconciled forecasts
#' are wished.
#' @param nn_type Non-negative method: "osqp" (\emph{default}) or "sntz"
#' (\emph{set-negative-to-zero}, only if \code{CCC = TRUE}) with exogenous
#' constraints (\code{const = "exo"}); "osqp" (\emph{default}), "KAnn"
#' (only \code{type == "M"}) or "sntz" with endogenous constraints
#' (\code{const = "endo"}).
#' @param settings Settings for \pkg{osqp} (object \code{\link[osqp]{osqpSettings}}).
#' The default options are: \code{verbose = FALSE}, \code{eps_abs = 1e-5},
#' \code{eps_rel = 1e-5}, \code{polish_refine_iter = 100} and \code{polish = TRUE}.
#' For details, see the \href{https://osqp.org/}{\pkg{osqp} documentation}
#' (Stellato et al., 2019).
#' @param ... Additional functional arguments passed to \link[FoReco]{htsrec} of
#' FoReco.
#'
#' @details
#' \strong{Cross-sectional hierarchies}
#'
#' To be as simple as possible, we fix the forecast horizon equal to 1.
#' Let the base forecasts be a vector \mjseqn{\widehat{\mathbf{y}} =
#' \left[\widehat{\mathbf{a}}' \; \widehat{\mathbf{b}}'\right]'}, where
#' vector \mjseqn{\widehat{\mathbf{a}}} consists of the sub-vectors forming each
#' of the upper \mjseqn{L} levels of the hierarchy/grouping:
#' \mjsdeqn{\widehat{\mathbf{a}} = \left[\begin{array}{c}
#'     \widehat{a}_1 \cr \widehat{\mathbf{a}}_2 \cr \vdots \cr \widehat{\mathbf{a}}_L
#'     \end{array}\right],}
#' where \mjseqn{\widehat{\mathbf{a}}_l}, \mjseqn{l=1,\ldots, L}, has dimension
#' \mjseqn{(n_l \times 1)} and \mjseqn{\sum_{l=1}^{L} n_l = n_a}. Denote
#' \mjseqn{\mathbf{C}_l} the \mjseqn{(n_l \times n_b)} matrix mapping the
#' bts into the level-\mjseqn{l} uts, then the aggregation matrix
#' \mjseqn{\mathbf{C}} may be written as
#' \mjsdeqn{\mathbf{C} = \left[\begin{array}{c}
#'   \mathbf{C}_1 \cr\mathbf{C}_2 \cr \vdots \cr \mathbf{C}_L
#'   \end{array}\right],}
#' where the generic matrix \mjseqn{\mathbf{C}_L} is (\mjseqn{n_L \times n_b}),
#' \mjseqn{l=1, \ldots, L}. Given a generic level \mjseqn{l}, the reconciled
#' forecasts coherent with the base forecasts of that level are the solution to
#' a quadratic minimization problem with linear
#' exogenous constraints (\code{const = "exo"}):
#' \mjsdeqn{\begin{array}{c}\widetilde{\mathbf{b}}_{l} = \arg\min_{\mathbf{b}}
#' \left(\mathbf{b} - \widehat{\mathbf{b}}\right)'\mathbf{W}_b^{-1}
#' \left(\mathbf{b} - \widehat{\mathbf{b}}\right) \quad \mbox{ s.t. }
#' \mathbf{C}_l\mathbf{b} = \widehat{\mathbf{a}}_l, \quad l=1, \ldots, L ,\cr
#' \Downarrow \cr
#' \widetilde{\mathbf{b}}_{l} = \widehat{\mathbf{b}} +
#' \mathbf{W}_b\mathbf{C}_l'\left(\mathbf{C}_l\mathbf{W}_b\mathbf{C}_l'
#' \right)^{-1}\left(\widehat{\mathbf{a}}_l -\mathbf{C}_l\widehat{\mathbf{b}}
#' \right), \quad l=1,\ldots,L,\end{array}}
#' where \mjseqn{\mathbf{W}_b} is a (\mjseqn{n_b \times n_b}) p.d. matrix
#' (in Hollyman et al., 2021, \mjseqn{\mathbf{W}_b} is diagonal).
#' If endogenous constraints (\code{const = "endo"}) are considered,
#' denote \mjseqn{\widehat{\mathbf{y}}_l =
#' \left[\widehat{\mathbf{a}}_l' \; \widehat{\mathbf{b}}'\right]'} and
#' \mjseqn{\mathbf{U}_l' = \left[\mathbf{I}_{n_l}'  \; \mathbf{C}_l'\right]'},
#' then
#' \mjsdeqn{\begin{array}{c}\widetilde{\mathbf{y}}_{l} = \arg\min_{\mathbf{y}_l}
#' \left(\mathbf{y}_l - \widehat{\mathbf{y}}_l\right)'\mathbf{W}_l^{-1}
#' \left(\mathbf{y}_l - \widehat{\mathbf{y}}_l\right) \quad \mbox{ s.t. }
#' \mathbf{U}_l'\mathbf{y}_l = \mathbf{0}, \quad l=1, \ldots, L ,\cr
#' \Downarrow \cr
#' \widetilde{\mathbf{y}}_{l} = \left(\mathbf{I}_{n_b+n_l} -
#' \mathbf{W}_l\mathbf{U}_l\left(\mathbf{U}_l'\mathbf{W}_l
#' \mathbf{U}_l\right)^{-1}\mathbf{U}_l'\right)\widehat{\mathbf{y}}_{l},
#' \quad l=1,...,L,
#' \end{array}}
#' where \mjseqn{\mathbf{W}_l} is a (\mjseqn{n_l + n_b \times n_l + n_b})
#' p.d. matrix.
#' Combined Conditional Coherent (CCC) forecasts are obtained by taking
#' the simple average of the \mjseqn{L} level conditional, and of the bottom-up
#' reconciled forecasts, respectively (Di Fonzo and Girolimetto, 2021):
#' \mjsdeqn{\widetilde{\mathbf{y}}_{CCC}=\frac{1}{L+1}\sum_{l=1}^{L+1} \mathbf{S}\widetilde{\mathbf{b}}_l,}
#' with \mjsdeqn{\widetilde{\mathbf{b}}_{L+1} = \widehat{\mathbf{b}}.}
#'
#' Non-negative reconciled forecasts may be obtained by setting \code{nn_type}
#' alternatively as:
#' \itemize{
#'   \item to apply non-negative constraints to each level:
#'   \itemize{
#'     \item \code{nn_type = "KAnn"} (only \code{const = "endo"})
#'     \item \code{nn_type = "osqp"}
#'   }
#'   \item to apply non-negative constraints only to the CCC forecasts:
#'   \itemize{
#'     \item \code{nn_type = "sntz"} ("set-negative-to-zero")
#'   }
#' }
#'
#' \strong{Temporal hierarchies}
#'
#' The extension to the case of \strong{temporal hierarchies} is quite simple.
#' Using the same notation as in \code{\link[FoReco]{thfrec}()}, the
#' following `equivalences' hold between the symbols used
#' for the level conditional cross-sectional reconciliation and the ones
#' used in temporal reconciliation:
#' \mjsdeqn{L \equiv p-1, \quad (n_a, n_b, n) \equiv (k^*, m, k^*+m),}
#' and
#' \mjsdeqn{\mathbf{C} \equiv \mathbf{K} , \;
#' \mathbf{C}_1 \equiv \mathbf{K}_1 = \mathbf{1}_m', \;
#' \mathbf{C}_2 \equiv \mathbf{K}_2 = \mathbf{I}_{\frac{m}{k_{p-1}}},\; \ldots, \;
#' \mathbf{C}_L \equiv \mathbf{K}_{p-1} = \mathbf{I}_{\frac{m}{k_{2}}} \otimes
#' \mathbf{1}_{k_{2}}',\; \mathbf{S} \equiv \mathbf{R}.}
#'
#' The description of the \strong{cross-temporal extension} is currently under progress.
#'
#' @return The function returns the level reconciled forecasts
#' in case of an elementary hierarchy with one level.
#' Otherwise it returns a list with
#' \item{\code{recf}}{Level Conditional Coherent (\code{CCC = FALSE}) forecasts
#' matrix/vector calculated as simple averages of upper level conditional
#' reconciled forecasts, with either endogenous or exogenous constraints.
#' If \code{CCC = TRUE} then it is the Combined Conditional Coherent matrix/vector,
#' as weighted averages of LCC and bottom-up reconciled forecasts.}
#' \item{\code{levrecf}}{list of level conditional reconciled forecasts (+ BU).}
#' \item{\code{info} (\code{type="osqp"})}{matrix with some useful indicators (columns)
#' for each forecast horizon \mjseqn{h} (rows): run time (\code{run_time}), number of iteration,
#' norm of primal residual (\code{pri_res}), status of osqp's solution (\code{status}) and
#' polish's status (\code{status_polish}).}
#'
#' @references
#' Di Fonzo, T., and Girolimetto, D. (2021), Cross-temporal forecast reconciliation:
#' Optimal combination method and heuristic alternatives, \emph{International Journal
#' of Forecasting}, in press.
#'
#' Di Fonzo, T., Girolimetto, D. (2021), Forecast combination based forecast reconciliation:
#' insights and extensions, \href{https://arxiv.org/abs/2106.05653}{arXiv:2006.08570}.
#'
#' Hollyman, R., Petropoulos, F., Tipping, M.E. (2021), Understanding Forecast Reconciliation,
#' \emph{European Journal of Operational Research} (in press).
#'
#' @family reconciliation procedures
#'
#' @examples
#' data(FoReco_data)
#' ### LCC and CCC CROSS-SECTIONAL FORECAST RECONCILIATION
#' # Cross sectional aggregation matrix
#' C <- rbind(FoReco_data$C, c(0,0,0,0,1))
#' # monthly base forecasts
#' id <- which(simplify2array(strsplit(colnames(FoReco_data$base), split = "_"))[1, ] == "k1")
#' mbase <- t(FoReco_data$base[, id])[,c("T", "A","B","C","AA","AB","BA","BB","C")]
#' # residuals
#' id <- which(simplify2array(strsplit(colnames(FoReco_data$res), split = "_"))[1, ] == "k1")
#' mres <- t(FoReco_data$res[, id])[,c("T", "A","B","C","AA","AB","BA","BB","C")]
#' # covariance matrix of all the base forecasts, computed using the in-sample residuals
#' Wres <- cov(mres)
#' # covariance matrix of the bts base forecasts, computed using the in-sample residuals
#' Wb <- Wres[c("AA","AB", "BA", "BB", "C"),
#'            c("AA","AB", "BA", "BB", "C")]
#' # bts seasonal averages
#' obs_1 <- FoReco_data$obs$k1
#' bts_sm <- apply(obs_1, 2, function(x) tapply(x[1:168],rep(1:12, 14), mean))
#' bts_sm <- bts_sm[,c("AA", "AB", "BA", "BB", "C")]
#'
#' ## EXOGENOUS CONSTRAINTS AND DIAGONAL COVARIANCE MATRIX (Hollyman et al., 2021)
#' # Forecast reconciliation for an elementary hierarchy:
#' # 1 top-level series + 5 bottom-level series (Level 2 base forecasts are not considered).
#' # The input is given by the base forecasts of the top and bottom levels series,
#' # along with a vector of positive weights for the bts base forecasts
#' exo_EHd <- lccrec(basef = mbase[, c("T","AA","AB", "BA", "BB", "C")],
#'                  weights = diag(Wb), bnaive = bts_sm)
#'
#' # Level conditional reconciled forecasts
#' # recf/L1: Level 1 reconciled forecasts for the whole hierarchy
#' # L2: Middle-Out reconciled forecasts hinging on Level 2 conditional reconciled forecasts
#' # L3: Bottom-Up reconciled forecasts
#' exo_LCd <- lccrec(basef = mbase, C = C, nl = c(1,3), weights = diag(Wb),
#'                  CCC = FALSE, bnaive = bts_sm)
#'
#' # Combined Conditional Coherent (CCC) reconciled forecasts
#' # recf: CCC reconciled forecasts matrix
#' # L1: Level 1 conditional reconciled forecasts for the whole hierarchy
#' # L2: Middle-Out reconciled forecasts hinging on Level 2 conditional reconciled forecasts
#' # L3: Bottom-Up reconciled forecasts
#' exo_CCCd <- lccrec(basef = mbase, C = C, nl = c(1,3), weights = diag(Wb))
#'
#' ## EXOGENOUS CONSTRAINTS AND FULL COVARIANCE MATRIX
#' # Simply substitute weights=diag(Wb) with weights=Wb
#' exo_EHf <- lccrec(basef = mbase[, c("T","AA","AB", "BA", "BB", "C")], weights = Wb, bnaive = bts_sm)
#' exo_LCf <- lccrec(basef = mbase, C = C, nl = c(1,3), weights = Wb, CCC = FALSE, bnaive = bts_sm)
#' exo_CCCf <- lccrec(basef = mbase, C = C, nl = c(1,3), weights = Wb, bnaive = bts_sm)
#'
#' ## ENDOGENOUS CONSTRAINTS AND DIAGONAL COVARIANCE MATRIX
#' # parameters of function htsrec(), like "type" and "nn_type" may be used
#'
#' # Forecast reconciliation for an elementary hierarchy with endogenous constraints
#' W1 <- Wres[c("T","AA","AB", "BA", "BB", "C"),
#'            c("T","AA","AB", "BA", "BB", "C")]
#' endo_EHd <- lccrec(basef = mbase[, c("T","AA","AB", "BA", "BB", "C")],
#'                  weights = diag(W1), const = "endo", CCC = FALSE)
#'
#' # Level conditional reconciled forecasts with endogenous constraints
#' endo_LCd <- lccrec(basef = mbase, C = C, nl = c(1,3), weights = diag(Wres),
#'                   const = "endo", CCC = FALSE)
#'
#' # Combined Conditional Coherent (CCC) reconciled forecasts with endogenous constraints
#' endo_CCCd <- lccrec(basef = mbase, C = C, nl = c(1,3),
#'                     weights = diag(Wres), const = "endo")
#'
#' ## ENDOGENOUS CONSTRAINTS AND FULL COVARIANCE MATRIX
#' # Simply substitute weights=diag(Wres) with weights=Wres
#' endo_EHf <- lccrec(basef = mbase[, c("T","AA","AB", "BA", "BB", "C")],
#'                    weights = W1,
#'                    const = "endo")
#' endo_LCf <- lccrec(basef = mbase, C = C, nl = c(1,3),
#'                    weights = Wres, const = "endo", CCC = FALSE)
#' endo_CCCf <- lccrec(basef = mbase-40, C = C, nl = c(1,3),
#'                    weights = Wres, const = "endo")
#'
#' ### LCC and CCC TEMPORAL FORECAST RECONCILIATION
#' # top ts base forecasts ([lowest_freq' ...  highest_freq']')
#' topbase <- FoReco_data$base[1, ]
#' # top ts residuals ([lowest_freq' ...  highest_freq']')
#' topres <- FoReco_data$res[1, ]
#' Om_bt <- cov(matrix(topres[which(simplify2array(strsplit(names(topres),
#'             "_"))[1,]=="k1")], ncol = 12, byrow = TRUE))
#' t_exo_LCd <- lccrec(basef = topbase, m = 12, weights = diag(Om_bt), CCC = FALSE)
#' t_exo_CCCd <- lccrec(basef = topbase, m = 12, weights = diag(Om_bt))
#'
#' ### LCC and CCC CROSS-TEMPORAL FORECAST RECONCILIATION
#' idr <- which(simplify2array(strsplit(colnames(FoReco_data$res), "_"))[1,]=="k1")
#' bres <- FoReco_data$res[-c(1:3), idr]
#' bres <- lapply(1:5, function(x) matrix(bres[x,], nrow=14, byrow = TRUE))
#' bres <- do.call(cbind, bres)
#' ctbase <- FoReco_data$base[c("T", "A","B","C","AA","AB","BA","BB","C"),]
#' ct_exo_LCf <- lccrec(basef = ctbase, m = 12, C = C, nl = c(1,3),
#'                     weights = diag(cov(bres)), CCC = FALSE)
#' ct_exo_CCCf <- lccrec(basef = ctbase, m = 12, C = C, nl = c(1,3),
#'                      weights = diag(cov(bres)), CCC = TRUE)
#'
#' @export
#' @import Matrix
#'
lccrec <- function(basef, m, C, nl, weights, bnaive = NULL,
                   const = "exogenous", CCC = TRUE, nn = FALSE,
                   nn_type = "osqp", settings = list(), ...){

  const <- match.arg(const, c("exogenous", "endogenous"))

  if(missing(C) | missing(m)){
    if(missing(m)){
      # lccrec cross-sectional
      if(missing(C)){
        C <- NULL
      }
      out <- lccrec_hts(basef = basef, C = C, nl = nl, weights = weights,
                        bnaive = bnaive, nn = nn, settings = settings,
                        const = const, CCC = CCC, nn_type = nn_type, ...)
    }else{
      # lccrec temporal
      out <- lccrec_thf(basef = basef, m = m, weights = weights,
                        bnaive = bnaive, nn = nn, settings = settings,
                        const = const, CCC = CCC, nn_type = nn_type, ...)
    }
  }else{
    # lccrec cross-temporal
    out <- lccrec_ctf(basef = basef, m = m, C = C, nl = nl, weights = weights,
                      bnaive = bnaive, nn = nn, settings = settings,
                      const = const, CCC = CCC, nn_type = nn_type, ...)
  }
  return(out)
}

lccrec_thf <- function(basef, m, weights, bnaive = NULL,
                       const = "exogenous", CCC = TRUE, nn = FALSE,
                       nn_type = "osqp", settings = list(), ...){
  # m condition
  if (missing(m)) {
    stop("The argument m is not specified", call. = FALSE)
  }
  tools <- thf_tools(m)
  kset <- tools$kset
  m <- tools$m
  p <- tools$p
  kt <- tools$kt
  ks <- tools$ks

  # matrix
  K <- tools$K
  R <- tools$R
  Zt <- tools$Zt

  # base forecasts condition
  if (missing(basef)) {
    stop("The argument basef is not specified", call. = FALSE)
  }

  if (NCOL(basef) != 1) {
    stop("basef must be a vector", call. = FALSE)
  }

  # Base Forecasts matrix
  h <- length(basef) / kt
  Dh <- Dmat(h = h, m = kset, n = 1)
  basef <- matrix(Dh %*% basef, nrow = h, byrow = TRUE)

  idr <- rep(1:length(kset[-1]), rev(kset[-1]))
  Kl <- lapply(unique(idr), function(i) K[idr==i, , drop = FALSE])

  check_bil <- any(sapply(Kl, function(x) any(colSums(x)!=1)))
  if(check_bil){
    warning("The hierarchy is not balanced. \nThis could produce results that are inconsistent with the reconciliation procedure.", call. = FALSE)
  }

  idc <- rep(1:length(kset), rev(kset))
  lal <- lapply(unique(idc), function(i) basef[,idc==i, drop = FALSE])
  b <- lal[[max(idc)]]
  lal <- lal[-max(idc)]

  if(!is.null(bnaive)){
    if(is.vector(bnaive)){
      if(length(bnaive) != m*h){
        stop("bnaive must be a vector (mh x 1)", call. = FALSE)
      }
    }else{
      stop("bnaive must be a vector (mh x 1)", call. = FALSE)
    }

    Db <- Dmat(h = h, m = m, n = 1)
    bm <- matrix(Db%*%bnaive, nrow = h, byrow = TRUE)
  }else{
    bm <- b
  }

  if(is.vector(weights)){
    weights <- Diagonal(x = weights)
  }

  nn_type <- match.arg(nn_type, c("osqp", "fbpp", "KAnn", "sntz"))
  if(nn_type == "sntz"){
    nn <- FALSE
  }

  if(const=="exogenous"){
    if(NCOL(weights) != m | NROW(weights) != m){
      stop("weights must be a m x m matrix", call. = FALSE)
    }else if(!is(weights, "Matrix")){
      weights <- Matrix(weights, sparse = TRUE)
    }
    out <- Map(function(al, cl){
      recoLEV(al = al, Wb = weights, bl = bm, cl = cl, nn = nn,
              settings = settings)
    }, al = lal, cl = Kl)

    bl <- lapply(out, function(x){
      extract_data(x = x, name = "recf")
    })

    info <- lapply(out, function(x){
      extract_data(x = x, name = "info")
    })
    names(info) <- paste0("id_k_", kset[-length(kset)])
    info <- info[!is.na(info)]

    uts0 <- do.call(c, lapply(out, function(x){
      extract_data(x = x, name = "uts0")
    }))
    if(any(uts0)){
      message("Some aggregation orders (greater than m) have base forecasts less then 0,",
              " then they have been set equal to zero")
    }
  }else{
    if(NCOL(weights) != kt | NROW(weights) != kt){
      stop("weights must be a (k* + m) x (k* + m) matrix", call. = FALSE)
    }else if(!is(weights, "Matrix")){
      weights <- Matrix(weights, sparse = TRUE)
    }

    id_nl <- rep(1:length(kset[-1]), rev(kset[-1]))
    Ol_id <- lapply(1:length(kset[-1]), function(x) c(which(id_nl == x), (ks+1):(m+ks)))
    out <- Map(function(al, cl, w){
      suppressMessages(htsrec(basef = cbind(al, bm),
                              C = cl,
                              comb = "w",
                              W = weights[w,w,drop = FALSE],
                              nn = nn,
                              settings = settings,
                              nn_type = nn_type, ...))
    }, al = lal, cl = Kl, w = Ol_id)

    bl <- lapply(out, function(x){
      out <- extract_data(x = x, name = "recf")
      out <- out[,-(1:(NCOL(out)-m)), drop = FALSE]
    })

    info <- lapply(out, function(x){
      extract_data(x = x, name = "info")
    })
    names(info) <- paste0("id_k_", kset[-length(kset)])
    info <- info[!is.na(info)]
  }

  if(nn){
    b[b<0] <- 0
  }
  bl[[length(bl)+1]] <- b

  yl <- lapply(bl, function(x){
    out <- as.matrix(x%*%t(R))
    return(out)
  })

  # return condition
  names(yl) <- paste0("id_k_", kset)
  if(CCC){
    out <- list()
    if(nn_type == "sntz"){
      b_sntz <- Reduce("+", bl)/length(bl)
      b_sntz <- b_sntz*(b_sntz>0)
      ccc_recf <- as.matrix(b_sntz%*%t(R))
    }else{
      ccc_recf <- Reduce("+", yl)/length(yl)
    }

    out$recf <- stats::setNames(as.vector(t(Dh) %*% as.vector(t(ccc_recf))), paste0("k", rep(kset, h * (m/kset)), "h",
                                                                                    do.call("c", as.list(sapply(
                                                                                      (m/kset) * h,
                                                                                      function(x) seq(1:x)
                                                                                    )))))
    out$levrecf <- lapply(yl, function(x){
      y <- as.vector(t(Dh) %*% as.vector(t(x)))
      y <- stats::setNames(y, paste0("k", rep(kset, h * (m/kset)), "h",
                                     do.call("c", as.list(sapply(
                                       (m/kset) * h,
                                       function(x) seq(1:x)
                                     )))))
      return(y)
    })
    if(length(info)>0){
      out$info <- info
    }
    return(out)
  }else{
    if(length(info)>0){
      lev <- lapply(yl, function(x){
        y <- as.vector(t(Dh) %*% as.vector(t(x)))
        y <- stats::setNames(y, paste0("k", rep(kset, h * (m/kset)), "h",
                                       do.call("c", as.list(sapply(
                                         (m/kset) * h,
                                         function(x) seq(1:x)
                                       )))))
        return(y)
      })
      out <- list()
      #out$recf <- lev[[1]]
      out$recf <- Reduce("+", lev[-length(lev)])/(length(lev)-1)
      out$levrecf <- lev
      out$info <- info
    }else if(nn_type == "sntz"){
      yl0 <- lapply(bl, function(x){
        out <- as.matrix((x*(x>0))%*%t(R))
        return(out)
      })
      lev <- lapply(yl0, function(x){
        y <- as.vector(t(Dh) %*% as.vector(t(x)))
        y <- stats::setNames(y, paste0("k", rep(kset, h * (m/kset)), "h",
                                       do.call("c", as.list(sapply(
                                         (m/kset) * h,
                                         function(x) seq(1:x)
                                       )))))
        return(y)
      })
      out <- list()
      #out$recf <- lev[[1]]
      out$recf <- Reduce("+", lev[-length(lev)])/(length(lev)-1)
      out$levrecf <- lev
    }else{
      lev <- lapply(yl, function(x){
        y <- as.vector(t(Dh) %*% as.vector(t(x)))
        y <- stats::setNames(y, paste0("k", rep(kset, h * (m/kset)), "h",
                                       do.call("c", as.list(sapply(
                                         (m/kset) * h,
                                         function(x) seq(1:x)
                                       )))))
        return(y)
      })
      out <- list()
      #out$recf <- lev[[1]]
      out$recf <- Reduce("+", lev[-length(lev)])/(length(lev)-1)
      out$levrecf <- lev
    }
    return(out)
  }
}

lccrec_hts <- function(basef, C, nl, weights, bnaive = NULL,
                       const = "exogenous", CCC = TRUE, nn = FALSE,
                       nn_type = "osqp", settings = list(), ...){

  if(is.vector(basef)){
    basef <- t(basef)
  }

  h <- NROW(basef)
  n <- NCOL(basef)

  if(is.null(C)){
    utd <- TRUE
    C <- list(Matrix(1, ncol = n-1))
    nl <- 1
  }else{
    utd <- FALSE
  }

  if(is.list(C)){
    nb <- NCOL(C[[1]])
    nl <- unlist(Map("NROW", C))
    na <- sum(nl)

    if(n != nb + na){
      stop("columns numb of basef != nb + na!", call. = FALSE)
    }
    Cl <- lapply(C, function(x) Matrix(x, sparse = TRUE))
    C <- do.call("rbind", Cl)
    S <- rbind(C, Diagonal(nb))
  }else{
    nb <- NCOL(C)
    na <- NROW(C)

    if(n != nb + na){
      stop("columns numb of basef != nb + na!", call. = FALSE)
    }

    if(missing(nl)){
      warning("nl? There is only one level of upper time series?", call. = FALSE)
      nl <- na
    }else if(sum(nl) != na){
      stop("Please, provide a valid nl vector s.t. sum(nl) == na", call. = FALSE)
    }else if(nl[1] != 1){
      stop("Please, provide a valid nl vector s.t. nl[1] == 1", call. = FALSE)
    }


    idr <- rep(1:length(nl), nl)
    if(!is(C, "Matrix")){
      C <- Matrix(C, sparse = TRUE)
    }

    Cl <- lapply(unique(idr), function(i) C[idr==i, , drop = FALSE])
    S <- rbind(C, Diagonal(nb))
  }

  check_bil <- any(sapply(Cl, function(x) any(colSums(x)!=1)))
  if(check_bil){
    warning("The hierarchy is not balanced. \nThis could produce results that are inconsistent with the reconciliation procedure.", call. = FALSE)
  }

  idc <- rep(1:(length(nl)+1), c(nl, nb))
  lal <- lapply(unique(idc), function(i) basef[,idc==i, drop = FALSE])
  b <- lal[[max(idc)]]
  lal <- lal[-max(idc)]

  if(!is.null(bnaive)){
    if(is.vector(bnaive)){
      bnaive <- rbind(bnaive)
    }
    if(is.matrix(bnaive)){
      if(any(dim(b) != dim(bnaive))){
        bnaive <- t(bnaive)
      }

      if(any(dim(b) != dim(bnaive))){
        stop("bnaive must be a matrix (h x nb)")
      }
    }else{
      stop("bnaive must be a matrix (h x nb)")
    }
    bm <- bnaive
  }else{
    bm <- b
  }

  if(is.vector(weights)){
    weights <- .sparseDiagonal(x = weights)
  }

  nn_type <- match.arg(nn_type, c("osqp", "fbpp", "KAnn", "sntz"))
  if(nn_type == "sntz"){
    nn <- FALSE
  }

  if(const=="exogenous"){
    if(NCOL(weights) != nb | NROW(weights) != nb){
      stop("weights must be a nb x nb matrix", call. = FALSE)
    }else if(!is(weights, "Matrix")){
      weights <- Matrix(weights, sparse = TRUE)
    }
    out <- Map(function(al, cl){
      recoLEV(al = al, Wb = weights, bl = bm, cl = cl, nn = nn,
              settings = settings)
    }, al = lal, cl = Cl)

    bl <- lapply(out, function(x){
      extract_data(x = x, name = "recf")
    })

    info <- lapply(out, function(x){
      extract_data(x = x, name = "info")
    })
    names(info) <- paste0("id_l_", 1:length(bl))
    info <- info[!is.na(info)]

    uts0 <- do.call(c, lapply(out, function(x){
      extract_data(x = x, name = "uts0")
    }))
    if(any(uts0)){
      message("Some upper times series have base forecasts less then 0,",
              " then they have been set equal to zero")
    }

  }else{
    if(NCOL(weights) != n | NROW(weights) != n){
      stop("weights must be a n x n matrix", call. = FALSE)
    }else if(!is(weights, "Matrix")){
      weights <- Matrix(weights, sparse = TRUE)
    }

    id_nl <- rep(1:length(nl), nl)
    Wl_id <- lapply(1:length(nl), function(x) c(which(id_nl == x), (na+1):(nb+na)))
    out <- Map(function(al, cl, w){
      suppressMessages(htsrec(basef = cbind(al, bm),
                              C = cl,
                              comb = "w",
                              W = weights[w,w,drop = FALSE],
                              nn = nn,
                              settings = settings,
                              nn_type = nn_type, ...))
    }, al = lal, cl = Cl, w = Wl_id)

    bl <- lapply(out, function(x){
      out <- extract_data(x = x, name = "recf")
      out <- out[,-(1:(NCOL(out)-nb)), drop = FALSE]
    })

    info <- lapply(out, function(x){
      extract_data(x = x, name = "info")
    })
    names(info) <- paste0("id_l_", 1:length(bl))
    info <- info[!is.na(info)]
  }

  if(!utd){
    if(nn){
      b[b<0] <- 0
    }
    bl[[length(bl)+1]] <- b
  }

  yl <- lapply(bl, function(x){
    out <- as.matrix(x%*%t(S))
    rownames(out) <- paste0("h", 1:h)
    colnames(out) <- if(is.null(colnames(basef))) paste("serie", 1:n, sep = "") else colnames(basef)
    return(out)
  })

  # return condition
  if(length(yl)==1){
    if(length(info)>0){
      return(list(recf = yl[[1]],
                  info = info))
    }else if(nn_type == "sntz"){
      b_sntz <- bl[[1]]*(bl[[1]]>0)
      y_sntz <- as.matrix(b_sntz%*%t(S))
      rownames(y_sntz) <- paste0("h", 1:h)
      colnames(y_sntz) <- if(is.null(colnames(basef))) paste("serie", 1:n, sep = "") else colnames(basef)
      return(y_sntz)
    }else{
      return(yl[[1]])
    }
  }else{
    names(yl) <- paste0("id_l_", 1:length(yl))
    if(CCC){
      out <- list()
      if(nn_type == "sntz"){
        b_sntz <- Reduce("+", bl)/length(bl)
        b_sntz <- b_sntz*(b_sntz>0)
        out$recf <- as.matrix(b_sntz%*%t(S))
        rownames(out$recf) <- paste0("h", 1:h)
        colnames(out$recf) <- if(is.null(colnames(basef))) paste("serie", 1:n, sep = "") else colnames(basef)
      }else{
        out$recf <- Reduce("+", yl)/length(yl)
      }

      out$levrecf <- yl
      if(length(info)>0){
        out$info <- info
      }
      return(out)
    }else{
      out <- list()
      if(length(info)>0){
        #out$recf <- yl[[1]]
        out$recf <- Reduce("+", yl[-length(yl)])/(length(yl)-1)
        out$levrecf <- yl
        out$info <- info
      }else if(nn_type == "sntz"){
        yl0 <- lapply(bl, function(x){
          out <- as.matrix((x*(x>0))%*%t(S))
          rownames(out) <- paste0("h", 1:h)
          colnames(out) <- if(is.null(colnames(basef))) paste("serie", 1:n, sep = "") else colnames(basef)
          return(out)
        })
        #out$recf <- yl0[[1]]
        out$recf <- Reduce("+", yl0[-length(yl0)])/(length(yl0)-1)
        out$levrecf <- yl0
      }else{
        #out$recf <- yl[[1]]
        out$recf <- Reduce("+", yl[-length(yl)])/(length(yl)-1)
        out$levrecf <- yl
      }
      return(out)
    }
  }
}

lccrec_ctf <- function(basef, C, nl, m, weights, bnaive = NULL,
                       const = "exogenous", CCC = TRUE, nn = FALSE,
                       nn_type = "osqp", settings = list(), ...){

  ctf <- suppressMessages(ctf_tools(C = C, m = m,  h = 1, sparse = TRUE))
  S <- ctf$ctf$Fmat

  n <- ctf$hts$n
  na <- ctf$hts$na
  nb <- ctf$hts$nb

  kset <- ctf$thf$kset
  m <- ctf$thf$m
  p <- ctf$thf$p
  kt <- ctf$thf$kt
  ks <- ctf$thf$ks

  if (NROW(basef) != n) {
    stop("Incorrect dimension of Ut or basef (they don't have same columns).", call. = FALSE)
  }

  # Base forecast
  if (NCOL(basef) %% kt != 0) {
    stop("basef has a number of row not in line with the frequency of the series.", call. = FALSE)
  }

  h <- NCOL(basef) / kt
  Dh <- Dmat(h = h, m = kset, n = n)
  Ybase <- matrix(Dh %*% as.vector(t(basef)), nrow = h, byrow = TRUE)

  if(missing(nl)){
    warning("nl? There is only one level of upper time series?", call. = FALSE)
    nl <- na
  }else if(sum(nl) != na){
    stop("Please, provide a valid nl vector s.t. sum(nl) == na", call. = FALSE)
  }else if(nl[1] != 1){
    stop("Please, provide a valid nl vector s.t. nl[1] == 1", call. = FALSE)
  }

  nk_nb <- c(nl, nb)
  nLk <- kronecker(nk_nb, m/kset)
  MLk <- matrix(1:length(nLk), ncol = length(kset),  byrow = TRUE)
  id <- as.vector(sapply(rep(split(MLk, 1:NROW(MLk)), nk_nb),
                         function(x) rep(x, m/kset)))
  lal <- lapply(unique(id), function(i) Ybase[,id==i, drop = FALSE])
  b <- lal[[max(id)]]
  lal <- lal[-max(id)]

  Cl <- lapply(1:(max(id)-1), function(x) S[id==x,, drop = FALSE])
  check_bil <- any(sapply(Cl, function(x) any(colSums(x)!=1)))
  if(check_bil){
    warning("The hierarchy is not balanced. \nThis could produce results that are inconsistent with the reconciliation procedure.", call. = FALSE)
  }

  names_lev <- t(sapply(unique(id), function(x) which(MLk == x, arr.ind = T)))
  names_lev[,2] <- kset[names_lev[,2]]

  if(!is.null(bnaive)){
    if(is.vector(bnaive)){
      bnaive <- rbind(bnaive)
    }
    if(is.matrix(bnaive)){
      if(any(c(nb, m*h) != dim(bnaive))){
        bnaive <- t(bnaive)
      }

      if(any(c(nb, m*h) != dim(bnaive))){
        stop("bnaive must be a matrix (nb x mh)", call. = FALSE)
      }
    }else{
      stop("bnaive must be a matrix (nb x mh)", call. = FALSE)
    }
    Db <- Dmat(h = h, m = m, n = nb)
    bm <- matrix(Db%*%as.vector(t(bnaive)), nrow = h, byrow = TRUE)
  }else{
    bm <- b
  }

  if(is.vector(weights)){
    weights <- Diagonal(x = weights)
  }

  nn_type <- match.arg(nn_type, c("osqp", "fbpp", "KAnn", "sntz"))
  if(nn_type == "sntz"){
    nn <- FALSE
  }

  if(const=="exogenous"){
    if(NCOL(weights) != nb*m | NROW(weights) != nb*m){
      stop("weights must be a nb*m x nb*m matrix", call. = FALSE)
    }else if(!is(weights, "Matrix")){
      weights <- Matrix(weights, sparse = TRUE)
    }
    out <- Map(function(al, cl){
      recoLEV(al = al, Wb = weights, bl = bm, cl = cl, nn = nn,
              settings = settings)
    }, al = lal, cl = Cl)

    bl <- lapply(out, function(x){
      extract_data(x = x, name = "recf")
    })

    info <- lapply(out, function(x){
      extract_data(x = x, name = "info")
    })
    names(info) <- paste0("id_l_", names_lev[-NROW(names_lev),1], "_k_", names_lev[-NROW(names_lev),2])
    info <- info[!is.na(info)]

    uts0 <- do.call(c, lapply(out, function(x){
      extract_data(x = x, name = "uts0")
    }))
    if(any(uts0)){
      message("Some times series (except for the high frequency bottom time ",
              "series) have base forecasts less then 0,",
              " then they have been set equal to zero")
    }
  }else{
    if(NCOL(weights) != kt*n | NROW(weights) != kt*n){
      stop("weights must be a kt*n x kt*n matrix", call. = FALSE)
    }else if(!is(weights, "Matrix")){
      weights <- Matrix(weights, sparse = TRUE)
    }

    Wl_id <- lapply(1:(max(id)-1), function(x)
      c(which(id == x), which(id == max(id))))
    out <- Map(function(al, cl, w){
      suppressMessages(htsrec(basef = cbind(al, bm),
                              C = cl,
                              comb = "w",
                              W = weights[w,w,drop = FALSE],
                              nn = nn,
                              settings = settings,
                              nn_type = nn_type, ...))
    }, al = lal, cl = Cl, w = Wl_id)

    bl <- lapply(out, function(x){
      out <- extract_data(x = x, name = "recf")
      out <- out[,-(1:(NCOL(out)-NCOL(b))), drop = FALSE]
    })

    info <- lapply(out, function(x){
      extract_data(x = x, name = "info")
    })
    names(info) <- paste0("id_l_", names_lev[-NROW(names_lev),1], "_k_", names_lev[-NROW(names_lev),2])
    info <- info[!is.na(info)]
  }

  if(nn){
    b[b<0] <- 0
  }
  bl[[length(bl)+1]] <- b

  yl <- lapply(bl, function(x){
    out <- as.matrix(x%*%t(S))
    return(out)
  })

  names(yl) <- paste0("id_l_", names_lev[,1], "_k_", names_lev[,2])

  if(CCC){
    out <- list()
    if(nn_type == "sntz"){
      b_sntz <- Reduce("+", bl)/length(bl)
      b_sntz <- b_sntz*(b_sntz>0)
      out$recf <- as.matrix(b_sntz%*%t(S))
    }else{
      out$recf <- Reduce("+", yl)/length(yl)
    }
    out$recf <- v2m_oct(Dh = Dh, y = out$recf, n = n,
                        nam = rownames(basef), m = m,
                        kset = kset, h = h)
    out$levrecf <- lapply(yl, v2m_oct, Dh = Dh, n = n,
                      nam = rownames(basef), m = m,
                      kset = kset, h = h)
    if(length(info)>0){
      out$info <- info
    }
    return(out)
  }else{
    out <- list()
    if(length(info)>0){
      yl <- lapply(yl, v2m_oct, Dh = Dh, n = n,
                   nam = rownames(basef), m = m,
                   kset = kset, h = h)
      #out$recf <- yl[[1]]
      out$recf <- Reduce("+", yl[-length(yl)])/(length(yl)-1)
      out$levrecf <- yl
      out$info <- info
    }else if(nn_type == "sntz"){
      yl <- lapply(bl, function(x){
        out <- as.matrix((x*(x>0))%*%t(S))
        return(out)
      })
      yl <- lapply(yl, v2m_oct, Dh = Dh, n = n,
                   nam = rownames(basef), m = m,
                   kset = kset, h = h)
      #out$recf <- yl[[1]]
      out$recf <- Reduce("+", yl[-length(yl)])/(length(yl)-1)
      out$levrecf <- yl
    }else{
      yl <- lapply(yl, v2m_oct, Dh = Dh, n = n,
                   nam = rownames(basef), m = m,
                   kset = kset, h = h)
      #out$recf <- yl[[1]]
      out$recf <- Reduce("+", yl[-length(yl)])/(length(yl)-1)
      out$levrecf <- yl
    }
    return(out)
  }
}

v2m_oct <- function(y, Dh, n, nam, m, kset, h){
  out <- matrix(t(Dh) %*% as.vector(t(y)), nrow = n, byrow = TRUE)
  rownames(out) <- if (is.null(nam)) paste("serie", 1:n, sep = "") else nam
  colnames(out) <- paste("k", rep(kset, h * (m/kset)), "h",
                         do.call("c", as.list(sapply(
                           (m/kset) * h,
                           function(x) seq(1:x)
                         ))),
                         sep = ""
  )
  return(out)
}

# Generic lccrec
recoLEV <- function(al, Wb, bl, cl, nn = FALSE, settings) {
  out <- list()
  rh <- cl%*%Wb%*%t(cl)
  lh <- t(al) - cl%*%t(bl)
  out$recf <- t(bl) + Wb%*%t(cl)%*%solveLin(rh, lh, verb = FALSE)
  out$recf <- t(out$recf)
  out$uts0 <- FALSE
  if(nn & any(out$recf < (-sqrt(.Machine$double.eps)))){

    if(isDiagonal(Wb)){
      P <- .sparseDiagonal(x = diag(Wb)^(-1))
    }else{
      R <- chol(Wb)
      P <- chol2inv(R)
    }
    id <- which(rowSums(out$recf<(-sqrt(.Machine$double.eps)))!=0)

    if(any(al[id,]<0)){
      out$uts0 <- TRUE
      al[id,] <- al[id,]*(al[id,]>0)
    }

    rec <- Map(function(ah, bh){
      lev_osqp(a = ah, b = bh, P = P, C = cl, settings = settings)
    }, ah = split(al[id,,drop = FALSE], id), bh = split(bl[id,,drop = FALSE], id))

    rec <- do.call("rbind",rec)
    out$recf[id,] <- do.call("rbind", rec[,"recf"])
    out$info <- do.call("rbind", rec[,"info"])
    colnames(out$info) <- c("obj_val", "run_time", "iter", "pri_res",
                            "status", "status_polish")
    rownames(out$info) <- id

    out$recf[-id,,drop=FALSE] <- out$recf[-id,,drop=FALSE] * (out$recf[-id,,drop=FALSE] > 0)
  }else if(nn){
    out$recf <- out$recf * (out$recf > 0)
  }

  return(out)
}

# lccrec osqp
lev_osqp <- function(a, b, P, C, settings) {
  nb <- NCOL(C)

  l <- c(a, rep(0, nb))
  u <- c(a, rep(Inf, nb))
  A <- rbind(C, .sparseDiagonal(nb))

  q <- (-1) * t(P) %*% as.vector(b)

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

  out$recf[which(out$recf < 0)] <- 0

  out$info <- c(rec$info$obj_val, rec$info$run_time, rec$info$iter, rec$info$pri_res,
                rec$info$status_val, rec$info$status_polish)

  return(out)
}
