##########################################################################
#### 7 - VCdeco -  weighted kdecop
##########################################################################

### The following function is version of Nagler's kdecop() function supporting
### weights. Currently, only the transformation estimator with log-quadratic 
### local likelihood estimation and nearest-neighbor bandwidths (TLL*nn) is 
### implemented (Nagler's default). The method is described in Geenenens et al., 2014. 
### It estimates the density of bivariate quasi-copula data that is twice transformed.
### First, the copula ranks are transformed to standard normal variable. Second, 
### a PCA is performed on the two standard normal variables. The LL*nn estimation
### is performend ont pca score variables. 

require("locfit")

wtd.kdecop  <-  function (udata, bw = NA, mult = 1, method = "TLL2nn", knots = 30, 
                          renorm.iter = 3L, info = TRUE, weights=NULL) 
{
  udata <- as.matrix(udata)
  n <- NROW(udata)
  d <- NCOL(udata)
  if(is.null(weights)){weights <- rep(1,n)}
  
  ### Method implemented?
  if (!(method %in% c("TLL1nn", "TLL2nn"))) 
    stop("method not implemented")
  
  ## This function checks if input data is correct
  kdecopula:::check_kdecop_input(as.list(environment()))
  if (any(is.na(bw))) 
    # This function selects the bandwidth depending on the method
    # bw <- kdecopula:::bw_select(udata, method)
    # Use directly the weighted bw choice function for TLL*nn
    bw <- bw_tll_nn_weighted(udata, weights, deg= as.numeric(substr(method,4, 4)))
  
  # This function multiplies the bandwidth if required
  bw <- kdecopula:::multiply_bw(bw, mult, method, d)
  
  # Checks if bw is correctly defined
  kdecopula:::check_bw(bw, method)
  
  # If linear local regression
  if (method %in% c("TLL1", "TLL2", "TLL1nn", "TLL2nn")) {
    #lfit <- kdecopula:::my_locfit(udata, bw, deg = as.numeric(substr(method, 
    #                                                     4, 4)))
    #use here instead
    lfit <- my_locfitnn_weighted(udata, weights,  bw$B, bw$alpha, bw$kappa, deg=as.numeric(substr(method,4,4)))
  }
  else {
    lfit <- NULL
  }
  grid <- kdecopula:::make_grid(knots, d)
  # evall_tll transforms the score data back to normal distributed data
  evalf <- kdecopula:::eval_func(method)
  object <- list(udata = udata, bw = bw, lfit = lfit, method = method)
  # Finally here, the results are returned
  vals <- array(evalf(grid$expanded, obj = object), dim = rep(knots, 
                                                              d))
  rm("object")
  vals <- kdecopula:::renorm2unif(vals, grid, renorm.iter)
  res <- list(udata = udata, grid = grid$pnts, estimate = vals, 
              bw = bw, mult = mult, method = method, weights=weights)
  class(res) <- "kdecopula"
  with_fit_info_weighted(res, info, lfit, weights)
}


#################################################################
####  bw_tll_nn_weighted:  weighted bandwith selection function 

# Idea of the procedure is double transformation:
# 1) density is estimate not on quasi-copula data but 
#    on normally transformed variable (i.e. qnorm(u) data)
# 2) optimal bandwith is selected on score values of principal
#    components of normally transformed variable. Idea:
#    score values are uncorrelated and we do not have to 
#    estimated a optimal bandwithd for the "correlation" entry
#    but only univariate ones in each direction of the scores.

bw_tll_nn_weighted <- function (udata,weights,deg) 
{
  zdata <- qnorm(udata)
  n <- nrow(zdata)
  d <- ncol(zdata)
  pca <- princomp(zdata,covmat=cov.wt(zdata,wt=weights))
  B <- unclass(pca$loadings)
  B <- B * sign(diag(B))
  qrs <- unclass(pca$scores)
  alphsq <- seq(nrow(zdata)^(-1/5), 1, l = 50)
  
  # Loop to find optimal alpha (smoothing parameter),
  # lscvplot function needs weights as well. 
  # LSCV: Least squares cross validation statistic
  # The statistic is conceptionalized in order to find
  # optimal smoothing parameter in 1 dimension
  opt <- function(i) {
    val <- locfit::lscvplot(~qrs[, i], weights=weights, alpha = alphsq, deg = deg, 
                    kern = "gauss", maxk = 512)$value
    mean(alphsq[which.min(val)])
  }
  alpha.vec <- sapply(1:d, opt)
  kappa <- alpha.vec[1]/alpha.vec
  dimnames(B) <- NULL
  if (deg == 1) {
    alpha <- n^(1/5 - d/(4 + d)) * alpha.vec[1]
  } else {
    alpha <- n^(1/9 - d/(8 + d)) * alpha.vec[1]
  }
  list(B = B, alpha = alpha, kappa = kappa)
}



#################################################################
#### my_locfitnn function

my_locfitnn_weighted <- function (udata, weights, B, alpha, kappa, deg) 
{
  zdata <- qnorm(udata)
  qrs <- t(solve(B) %*% t(zdata))
  gr <- do.call(expand.grid, lapply(1:ncol(qrs), function(i) c(-4, 
                                                               4)))
  qgr <- t(solve(B) %*% t(as.matrix(gr)))
  lims <- apply(qgr, 2L, range)
  cl.lst <- split(as.vector(qrs), rep(1:ncol(qrs), each = nrow(qrs)))
  cl.lst$nn <- alpha
  cl.lst$deg <- deg
  cl.lst$scale <- kappa
  
  # The sample weights must be introduced in in the lf.lst object after the
  lf.lst <- list(~do.call(locfit::lp, cl.lst), weights=weights, maxk = 1024, kern = "gauss", 
                 ev = lfgrid(mg = 50, ll = lims[1L, ], ur = lims[2L, 
                                                                 ]))
  suppressWarnings(do.call(locfit::locfit, lf.lst))
}

#################################################################
#### bw_tll_nn bandwith selection function 

with_fit_info_weighted <- function (res, info, lfit, weights) 
{
  if (info) {
    likvalues <- dkdecop(res$udata, res)
    loglik <- sum(log(likvalues)*weights)
    if (!(res$method %in% c("TTPI", "TTCV"))) {
      effp <- kdecopula:::eff_num_par(res$udata, likvalues, res$bw, 
                          res$method, lfit)
      AIC <- -2 * loglik + 2 * effp
      n <- nrow(res$udata)
      cAIC <- AIC + (2 * effp * (effp + 1))/(n - effp - 
                                               1)
      BIC <- -2 * loglik + log(n) * effp
    }
    else {
      effp <- AIC <- cAIC <- BIC <- NA
    }
    res$info <- list(likvalues = likvalues, loglik = loglik, 
                     effp = effp, AIC = AIC, cAIC = cAIC, BIC = BIC)
  }
  res
}


