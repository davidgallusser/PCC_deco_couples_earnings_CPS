##########################################################################
#### 10 - VCdeco -  weighted kdvinecop
##########################################################################



data(wdbc, package = "kdecopula")
# rank-transform to copula data (margins are uniform)
u <- VineCopula::pobs(wdbc[, 5:7], ties = "average")

fit <- kdevinecop(u, method = "TLL2nn")   
dkdevinecop(c(0.1, 0.1, 0.1), fit)
contour(fit)
fit$matrix
summary(fit)

fit2 <- wtd.kdevinecop(u)
dkdevinecop(c(0.1, 0.1, 0.1), fit2)

fit3 <- wtd.kdevinecop(u, weights=runif(nrow(u),0.5,1.5),
                       treecrit="rho", info=TRUE)
dkdevinecop(c(0.1, 0.1, 0.1), fit3)
rkdevinecop(10,fit3, quasi=TRUE)

setwd(paste0(sav_d,datafiles_d))
load(paste0("CPScouples",2016,".rda"))
jit <- c(2,4)
data1 <- cdata[sample(1:nrow(cdata),1000),c("EARNINGS","WAGE_m","ANNUALHOURS_m", "WAGE_f","ANNUALHOURS_f","ASECWTH_m")]
w <- data1[,"ASECWTH_m"]
names(data1)[6] <- ("(weights)")
data1[,jit+1] <- apply(data1[,jit+1], 2, function(x) add_jitter(x))
u1  <- apply(data1[,2:5], 2, function(x) rank.wt(x,w))

VCfit1 <- wtd.kdevinecop(u1, weights=w,
               treecrit="rho", info=TRUE)
m <- VCfit1$matrix

VCfit1 <- wtd.kdevinecop(u1, matrix=m, weights=w,
                         treecrit="rho", info=TRUE)

load(paste0("CPScouples",1976,".rda"))
jit <- c(2,4)
data0 <- cdata[sample(1:nrow(cdata),1000),c("EARNINGS","WAGE_m","ANNUALHOURS_m", "WAGE_f","ANNUALHOURS_f","ASECWTH_m")]
w <- data0[,"ASECWTH_m"]
names(data0)[6] <- ("(weights)")
data0[,jit+1] <- apply(data0[,jit+1], 2, function(x) add_jitter(x))
u0  <- apply(data0[,2:5], 2, function(x) rank.wt(x,w))

VCfit0 <- wtd.kdevinecop(u0, matrix = m, weights=w,
                         treecrit="rho", info=TRUE)
VCfit0$matrix
VCfit0$info$cAIC


require("kdevine")
getAnywhere(vinecop)
library(rvinecopulib)
fit <- vinecop(u0,family_set = "tll", nonpar_method="quadratic", weights=w, tree_crit="rho", show_trace=TRUE)
plot(fit$pair_copulas[[1]][[3]])
summary(fit)

fit2 <- vinecop(u0,family_set = "tll", nonpar_method="quadratic", tree_crit="rho", show_trace=TRUE)
summary(fit2)

fit$controls
plot(fit)
rvinecopulib:::vinecop_select_cpp

##########################################################################
## weighted kdvinecop
wtd.kdevinecop <- function (data, matrix = NA, method = "TLL2nn", renorm.iter = 3L, 
                            mult = 1, test.level = NA, trunc.level = NA, treecrit = "tau", 
                            cores = 1, info = FALSE, weights=NULL) 
{
  if (is.null(info)) 
    info <- FALSE
  if (is.null(matrix)) 
    matrix <- NA
  if (is.null(method)) 
    method <- "TLL2"
  if (is.null(mult)) 
    mult <- 1
  if (is.null(test.level)) 
    test.level <- 1
  if (is.na(test.level)) 
    test.level <- 1
  if (is.null(trunc.level)) 
    trunc.level <- ncol(data)
  if (is.na(trunc.level)) 
    trunc.level <- ncol(data)
  if (is.null(treecrit)) 
    treecrit <- "tau"
  if (is.na(treecrit)) 
    treecrit <- "tau"
  if (is.null(cores)) 
    cores <- 1
  if (is.null(weights)) # Add weights
    weights <- rep(1,nrow(data))
  data <- as.matrix(data)
  data <- pobs(data, ties.method = "first")
  matrix <- as.matrix(matrix)
  d <- ncol(data)
  n <- nrow(data)
  if (any(data >= 1) || any(data < 0)) 
    stop("Data has be in the interval [0,1].")
  if (n < 2) 
    stop("Number of observations has to be at least 2.")
  if (d < 2) 
    stop("Dimension has to be at least 2.")
  if (!(treecrit %in% c("tau", "rho", "AIC", "cAIC", 
                        "hoeffd"))) 
    stop("'treecrit' not available; please choose either 'tau', 'AIC', 'cAIC', or 'hoeffd'")
  if (any(is.na(matrix)) & d > 2) {
    ################################################3
    # Needs weights:
    return(wtd.structure_select2(data, type = 0, method = method, 
                                       mult = mult, info = info, struct.crit = treecrit, 
                                       test.level = test.level, trunc.level = trunc.level, 
                                       renorm.iter = renorm.iter, cores = cores, weights=weights))
  }
  else if (any(is.na(matrix)) & d == 2) {
    matrix <- matrix(c(2, 1, 0, 1), 2, 2)
  }
  else if (nrow(matrix) == 1 && matrix == 1) {
    ################################################
    # Needs weights:
    return(wtd.structure_select2(data, type = 1, method = method, 
                                    mult = mult, info = info, struct.crit = treecrit, 
                                    test.level = test.level, trunc.level = trunc.level, 
                                    renorm.iter = renorm.iter, cores = cores, weights=weights))
  }
  if (nrow(matrix) != ncol(matrix)) 
    stop("Structure matrix has to be quadratic.")
  if (ncol(data) != ncol(matrix)) 
    stop("Dimensions of data and matrix don't match.")
  if (max(matrix) > nrow(matrix)) 
    stop("Error in the structure matrix.")
  if (cores != 1 | is.na(cores)) {
    if (is.na(cores)) 
      cores <- max(1, detectCores() - 1)
    if (cores > 1) {
      cl <- makeCluster(cores)
      registerDoParallel(cl)
      on.exit(try(stopCluster(), silent = TRUE))
      on.exit(try(closeAllConnections(), silent = TRUE), 
              add = TRUE)
    }
  }
  M <- kdevine:::ToLowerTri(matrix)
  Mold <- M
  o <- diag(M)
  M <- kdevine:::reorderRVineMatrix(M)
  data <- data[, o[length(o):1]]
  MaxMat <- kdevine:::createMaxMat(M)
  CondDistr <- kdevine:::neededCondDistr(M)
  res <- as.list(numeric(d - 1))
  for (i in 1:(d - 1)) res[[i]] <- as.list(numeric(d - i))
  llikv <- array(0, dim = c(d, d, n))
  llik <- matrix(0, d, d)
  effp <- matrix(0, d, d)
  AIC <- matrix(0, d, d)
  cAIC <- matrix(0, d, d)
  BIC <- matrix(0, d, d)
  nms <- matrix("", d - 1, d - 1)
  V <- list()
  V$direct <- array(NA, dim = c(d, d, n))
  V$indirect <- array(NA, dim = c(d, d, n))
  V$direct[d, , ] <- t(data[, d:1])
  indepinfo <- list(effp = 0, likvalues = rep(1, n), loglik = 0, 
                    effp = 0, AIC = 0, cAIC = 0, BIC = 0)
  for (k in d:2) {
    doEst <- function(i) {
      if (k > i) {
        m <- MaxMat[k, i]
        zr1 <- V$direct[k, i, ]
        zr2 <- if (m == M[k, i]) {
          V$direct[k, (d - m + 1), ]
        }
        else {
          V$indirect[k, (d - m + 1), ]
        }
        samples <- cbind(zr2, zr1)
        indep <- ifelse(test.level < 1, BiCopIndTest(zr2, 
                                                     zr1)$p.value >= test.level, FALSE)
        if (trunc.level <= (d - k)) 
          indep <- TRUE
        if (indep) {
          cfit <- list()
          if (info) 
            cfit$info <- indepinfo
          class(cfit) <- c("kdecopula", "indep.copula")
        }
        else {
          cfit <- wtd.kdecop(samples, mult = mult, method = method, 
                            renorm.iter = renorm.iter, info = info,
                            weights=weights)
        }
        hfit <- list()
        direct <- indirect <- NULL
        if (CondDistr$direct[k - 1, i]) {
          direct <- hkdecop(samples, obj = cfit, cond.var = 1L)
        }
        if (CondDistr$indirect[k - 1, i]) {
          indirect <- hkdecop(samples, obj = cfit, cond.var = 2L)
        }
        names <- kdevine:::naming(Mold[c(i, k:d), i])
        res.ki <- list(c = cfit, name = names)
        return(list(direct = direct, indirect = indirect, 
                    res.ki = res.ki))
      }
      else {
        return(NULL)
      }
    }
    res.k <- if (cores > 1) {
      foreach(i = 1:(k - 1), .export = c("naming"), 
              .packages = c("kdevine", "kdecopula")) %dopar% 
        doEst(i)
    }
    else lapply(1:(k - 1), doEst)
    for (i in 1:(d - 1)) {
      nums <- Mold[c(i, k:d), i]
      name <- kdevine:::naming(nums)
      if (any(kdevine:::extract_nums(res.k) == name)) {
        ind <- which(kdevine:::extract_nums(res.k) == name)
        res.ki <- res.k[[ind]]
        res[[d + 1 - k]][[i]] <- res.ki$res.ki
        if (!is.null(res.ki$direct)) {
          V$direct[k - 1, i, ] <- res.ki$direct
        }
        if (!is.null(res.ki$indirect)) {
          V$indirect[k - 1, i, ] <- res.ki$indirect
        }
        if (info) {
          cfit <- res.ki$res.ki$c
          llikv[k, i, ] <- log(cfit$info$likvalues)
          llik[k, i] <- cfit$info$loglik
          effp[k, i] <- cfit$info$effp
          AIC[k, i] <- -2 * cfit$info$loglik + 2 * effp[k, 
                                                        i]
          cAIC[k, i] <- AIC[k, i] + (2 * effp[k, i] * 
                                       (effp[k, i] + 1))/(n - effp[k, i] - 1)
          BIC[k, i] <- -2 * cfit$info$loglik + log(n) * 
            effp[k, i]
        }
      }
    }
  }
  res[[d]] <- data
  res[[d + 1]] <- Mold
  res[[d + 2]] <- list(llikv = apply(llikv, 3, sum), loglik = sum(llik), 
                       pair.loglik = llik, effp = sum(effp), pair.effp = effp, 
                       AIC = sum(AIC), pair.AIC = AIC, cAIC = sum(AIC) + (2 * 
                                                                            sum(effp) * (sum(effp) + 1))/(n - sum(effp) - 1), 
                       pair.cAIC = cAIC, BIC = sum(BIC), pair.BIC = BIC)
  names(res) <- vapply(1:(d - 1), function(x) paste("T", 
                                                    x, sep = ""), "")
  names(res)[d:(d + 2)] <- c("data", "matrix", 
                             "info")
  if (!info) 
    res[[d + 2]] <- NULL
  class(res) <- "kdevinecop"
  res
}


##########################################################################
## weighted structure_select2 (kdevine:::structure_select2)
# Applies Dissmann algorithm
wtd.structure_select2 <- function(data, type, method, mult, struct.crit, test.level, 
                                  trunc.level, renorm.iter, cores, info, progress = FALSE, 
                                  weights) {
  if (type == 0) {
    type <- "RVine"
  }
  else if (type == 1) {
    type <- "CVine"
  }
  if (type != "RVine" & type != "CVine") 
    stop("Vine model not implemented.")
  n <- nrow(data)
  d <- ncol(data)
  if (n < 2) 
    stop("Number of observations has to be at least 2.")
  if (d < 2) 
    stop("Dimension has to be at least 2.")
  if (any(data > 1) || any(data < 0)) 
    stop("Data has be in the interval [0,1].")
  if (is.na(test.level)) 
    test.level <- 1
  if (!(struct.crit %in% c("tau", "rho", "AIC", "cAIC", 
                           "hoeffd"))) 
    stop("'struct.crit' has to be one of 'tau', 'rho', 'AIC', 'cAIC', or 'hoeffd'")
  if (is.null(colnames(data))) 
    colnames(data) = paste("V", 1:d, sep = "")
  RVine <- list(Tree = NULL, Graph = NULL)
  res <- as.list(numeric(d - 1))
  for (i in 1:(d - 1)) res[[i]] <- as.list(numeric(d - i))
  llikv <- array(0, dim = c(d, d, n))
  llik <- matrix(0, d, d)
  effp <- matrix(0, d, d)
  AIC <- matrix(0, d, d)
  cAIC <- matrix(0, d, d)
  BIC <- matrix(0, d, d)
  if (cores != 1 | is.na(cores)) {
    if (is.na(cores)) 
      cores <- max(1, detectCores() - 1)
    if (cores > 1) {
      cl <- makeCluster(cores)
      registerDoParallel(cl)
      on.exit(try(stopCluster(), silent = TRUE))
      on.exit(try(closeAllConnections(), silent = TRUE), 
              add = TRUE)
    }
  }
  g <- wtd.initializeFirstGraph2(data, struct.crit = struct.crit, 
                             weights = weights)
  mst <- kdevine:::findMaxTree2(g, mode = type)
  est <- wtd.est.FirstTreeCopulas2(mst, data, method = method, 
                                   mult = mult, info = info, test.level = test.level, renorm.iter = renorm.iter, 
                                   parallel = cores > 1, weights = weights)
  res[[1]] <- est$est
  VineTree <- est$mst
  RVine$Tree[[1]] <- VineTree
  RVine$Graph[[1]] <- g
  for (k in 2:(d - 1)) {
    g <- wtd.buildNextGraph2(VineTree, weights = weights, struct.crit = struct.crit, 
                         parallel = cores > 1)
    mst <- kdevine:::findMaxTree2(g, mode = type)
    est <- wtd.est.TreeCopulas2(mst, k = k, d2 = d, data = data, 
                            oldVineGraph = VineTree, method = method, mult = mult, 
                            info = info, test.level = test.level, renorm.iter = renorm.iter, 
                            parallel = cores > 1, truncate = trunc.level <= k, weights = weights)
    res[[k]] <- est$est
    VineTree <- est$tree
    RVine$Tree[[k]] <- VineTree
    RVine$Graph[[k]] <- g
  }
  res[[d]] <- data
  res[[d + 1]] <- M <- kdevine:::as.RVMKernel2(RVine)$Matrix
  for (i in 1:(d - 1)) {
    for (k in (i + 1):d) {
      hfit <- res[[i]][[k - i]]$h
      if (info) {
        pcfit <- res[[i]][[k - i]]$c
        llikv[k, i, ] <- log(pcfit$info$likvalues)
        llik[k, i] <- pcfit$info$loglik
        effp[k, i] <- pcfit$info$effp
        AIC[k, i] <- -2 * pcfit$info$loglik + 2 * effp[k, 
                                                       i]
        cAIC[k, i] <- AIC[k, i] + (2 * effp[k, i] * (effp[k, 
                                                          i] + 1))/(n - effp[k, i] - 1)
        BIC[k, i] <- -2 * pcfit$info$loglik + log(n) * 
          effp[k, i]
      }
    }
  }
  res[[d + 2]] <- if (info) {
    list(llikv = apply(llikv, 3, sum), loglik = sum(llik), 
         pair.loglik = llik, effp = sum(effp), pair.effp = effp, 
         AIC = sum(AIC), pair.AIC = AIC, cAIC = sum(AIC) + 
           (2 * sum(effp) * (sum(effp) + 1))/(n - sum(effp) - 
                                                1), pair.cAIC = cAIC, BIC = sum(BIC), pair.BIC = BIC)
  } else list(NULL)
  names(res) <- vapply(1:(d - 1), function(x) paste("T", 
                                                    x, sep = ""), "")
  names(res)[d:(d + 2)] <- c("data", "matrix", 
                             "info")
  class(res) <- "kdevinecop"
  res
}

##########################################################################
## weighted est.FirstTreeCopulas2 (kdevine:::est.FirstTreeCopulas2)
# Estimates copula of first tree in Dissman algorithm
wtd.est.FirstTreeCopulas2 <- function (mst, data.univ, method, mult, test.level, renorm.iter, 
                                       info, parallel, weights) 
{
  d <- nrow(mst$E$nums)
  pkgs <- c("kdevine", "kdecopula")
  indepinfo <- list(effp = 0, likvalues = rep(1, nrow(data.univ)), 
                    loglik = 0, effp = 0, AIC = 0, cAIC = 0, BIC = 0)
  doEst <- function(i) {
    a <- rev(mst$E$nums[i, ])
    Copula.Data.1 = list(data.univ[, a[1]])
    Copula.Data.2 = list(data.univ[, a[2]])
    if (is.null(mst$V$names[a[1]])) {
      Copula.CondName.1 <- a[1]
    }
    else {
      Copula.CondName.1 <- mst$V$names[a[1]]
    }
    if (is.null(mst$V$names[a[2]])) {
      Copula.CondName.2 <- a[2]
    }
    else {
      Copula.CondName.2 <- mst$V$names[a[2]]
    }
    if (is.null(mst$V$names[a[1]]) || is.null(mst$V$names[a[2]])) {
      Copula.Name <- paste(a[1], a[2], sep = " , ")
    }
    else {
      Copula.Name <- paste(mst$V$names[a[1]], mst$V$names[a[2]], 
                           sep = " , ")
    }
    nums <- paste(a[2], a[1], sep = ",")
    s <- cbind(data.univ[, a[2]], data.univ[, a[1]])
    indep <- ifelse(test.level < 1, BiCopIndTest(s[, 1], 
                                                 s[, 2])$p.value >= test.level, FALSE)
    if (is.na(indep)) 
      browser()
    if (indep) {
      pcfit <- list()
      if (info) 
        pcfit$info <- indepinfo
      class(pcfit) <- c("kdecopula", "indep.copula")
    }
    else {
      pcfit <- wtd.kdecop(s, method = method, mult = mult, 
                      renorm.iter = renorm.iter, info = info, weights=weights)
    }
    if (indep == TRUE) {
      Copula.CondData.1 <- s[, 1]
      Copula.CondData.2 <- s[, 2]
    }
    else {
      Copula.CondData.1 <- list(hkdecop(s, obj = pcfit, 
                                        cond.var = 2L))
      Copula.CondData.2 <- list(hkdecop(s, obj = pcfit, 
                                        cond.var = 1L))
    }
    resi <- list(c = pcfit, name = nums)
    list(Copula.Data.1 = Copula.Data.1, Copula.Data.2 = Copula.Data.2, 
         Copula.CondName.1 = Copula.CondName.1, Copula.CondName.2 = Copula.CondName.2, 
         Copula.Name = Copula.Name, Copula.CondData.1 = Copula.CondData.1, 
         Copula.CondData.2 = Copula.CondData.2, resi = resi)
  }
  res <- if (parallel) {
    foreach(i = 1:d, .export = c("d"), .packages = pkgs) %dopar% 
      doEst(i)
  } else {
    lapply(1:d, doEst)
  }
  for (i in 1:d) {
    mst$E$Copula.Data.1[[i]] = res[[i]]$Copula.Data.1
    mst$E$Copula.Data.2[[i]] = res[[i]]$Copula.Data.2
    mst$E$Copula.CondName.1[[i]] = res[[i]]$Copula.CondName.1
    mst$E$Copula.CondName.2[[i]] = res[[i]]$Copula.CondName.2
    mst$E$Copula.Name[[i]] = res[[i]]$Copula.Name
    mst$E$Copula.CondData.1[[i]] = res[[i]]$Copula.CondData.1
    mst$E$Copula.CondData.2[[i]] = res[[i]]$Copula.CondData.2
    mst$E$Copula.Nums.1 = res[[i]]$resi$name
    res[[i]] <- res[[i]]$resi
  }
  list(mst = mst, est = res)
}

##########################################################################
## weighted est.TreeCopulas2 (kdevine:::est.TreeCopulas2)
# Estimates following trees in Dissman algorithm
wtd.est.TreeCopulas2 <- function (mst, k, d2, data, oldVineGraph, method, mult, info, 
                                  test.level, renorm.iter, weights, parallel, truncate){
  d <- nrow(mst$E$nums)
  indepinfo <- list(effp = 0, likvalues = rep(1, nrow(data)), 
                    loglik = 0, effp = 0, AIC = 0, cAIC = 0, BIC = 0)
  exp <- c("split_name", "split_num", "naming")
  pkgs <- c("kdecopula")
  doEst <- function(i, mst) {
    con <- rev(mst$E$nums[i, ])
    temp <- oldVineGraph$E$nums[con, ]
    if ((temp[1, 1] == temp[2, 1]) || (temp[1, 2] == temp[2, 
                                                          1])) {
      same <- temp[2, 1]
    }
    else {
      if ((temp[1, 1] == temp[2, 2]) || (temp[1, 2] == 
                                         temp[2, 2])) {
        same <- temp[2, 2]
      }
    }
    other1 <- temp[1, temp[1, ] != same]
    other2 <- temp[2, temp[2, ] != same]
    if (temp[1, 1] == same) {
      zr1 <- oldVineGraph$E$Copula.CondData.2[[con[1]]]
      n1 <- oldVineGraph$E$Copula.CondName.2[[con[1]]]
    }
    else {
      zr1 <- oldVineGraph$E$Copula.CondData.1[[con[1]]]
      n1 <- oldVineGraph$E$Copula.CondName.1[[con[1]]]
    }
    if (temp[2, 1] == same) {
      zr2 <- oldVineGraph$E$Copula.CondData.2[[con[2]]]
      n2 <- oldVineGraph$E$Copula.CondName.2[[con[2]]]
    }
    else {
      zr2 <- oldVineGraph$E$Copula.CondData.1[[con[2]]]
      n2 <- oldVineGraph$E$Copula.CondName.1[[con[2]]]
    }
    zr1a <- if (is.list(zr1)) 
      as.vector(zr1[[1]])
    else zr1
    zr2a <- if (is.list(zr2)) 
      as.vector(zr2[[1]])
    else zr2
    n1a <- if (is.list(n1)) 
      as.vector(n1[[1]])
    else n1
    n2a <- if (is.list(n2)) 
      as.vector(n2[[1]])
    else n2
    Copula.Data.1 <- zr1a
    Copula.Data.2 <- zr2a
    Copula.CondName.2 <- n1a
    Copula.CondName.1 <- n2a
    samples <- cbind(zr1a, zr2a)
    indep <- ifelse(test.level < 1, BiCopIndTest(samples[, 
                                                         1], samples[, 2])$p.value >= test.level, FALSE)
    if (truncate) 
      indep <- TRUE
    if (indep) {
      pcfit <- list()
      if (info) 
        pcfit$info <- indepinfo
      class(pcfit) <- c("kdecopula", "indep.copula")
    }
    else {
      pcfit <- wtd.kdecop(samples, mult = mult, method = method, 
                      renorm.iter = renorm.iter, info = info, weights=weights)
    }
    tmpname <- mst$E$names[i]
    namesplt <- kdevine:::split_name(tmpname)
    numsplt <- sapply(namesplt, function(x) which(colnames(data) == 
                                                    x))
    nums <- kdevine:::naming(sprintf("%d", numsplt))
    if (indep == TRUE) {
      Copula.CondData.1 <- samples[, 2]
      Copula.CondData.2 <- samples[, 1]
    }
    else {
      Copula.CondData.1 <- hkdecop(samples, obj = pcfit, 
                                   cond.var = 1L)
      Copula.CondData.2 <- hkdecop(samples, obj = pcfit, 
                                   cond.var = 2L)
    }
    resi <- list(c = pcfit, name = nums)
    list(Copula.Data.1 = Copula.Data.1, Copula.Data.2 = Copula.Data.2, 
         Copula.CondName.1 = Copula.CondName.1, Copula.CondName.2 = Copula.CondName.2, 
         Copula.CondData.1 = Copula.CondData.1, Copula.CondData.2 = Copula.CondData.2, 
         resi = resi)
  }
  res.k <- if (parallel) {
    foreach(i = 1:d, .export = exp, .packages = pkgs) %dopar% 
      doEst(i, mst)
  }
  else {
    lapply(1:d, doEst, mst = mst)
  }
  for (i in 1:d) {
    mst$E$Copula.Data.1[[i]] = res.k[[i]]$Copula.Data.1
    mst$E$Copula.Data.2[[i]] = res.k[[i]]$Copula.Data.2
    mst$E$Copula.CondName.1[[i]] = res.k[[i]]$Copula.CondName.1
    mst$E$Copula.CondName.2[[i]] = res.k[[i]]$Copula.CondName.2
    mst$E$Copula.CondData.1[[i]] = res.k[[i]]$Copula.CondData.1
    mst$E$Copula.CondData.2[[i]] = res.k[[i]]$Copula.CondData.2
    res.k[[i]] <- res.k[[i]]$resi
  }
  list(tree = mst, est = res.k)
}

##########################################################################
## getEdgeInfo2() with weights and Spearman rho, kdevine:::getEdgeInfo2() 
# Computes weight of an edge, Spearmans rho as tree crit added
wtd.getEdgeInfo2  <- function (i, g, oldVineGraph, weights, struct.crit = "tau") 
{
  con <- g$E$nums[i, ]
  temp <- oldVineGraph$E$nums[con, ]
  ok <- FALSE
  if ((temp[1, 1] == temp[2, 1]) || (temp[1, 2] == temp[2, 
                                                        1])) {
    ok <- TRUE
    same <- temp[2, 1]
  }
  else {
    if ((temp[1, 1] == temp[2, 2]) || (temp[1, 2] == temp[2, 
                                                          2])) {
      ok <- TRUE
      same <- temp[2, 2]
    }
  }
  w <- nedSet <- ningSet <- name <- NA
  todel <- TRUE
  if (ok) {
    if (temp[1, 1] == same) {
      zr1 <- oldVineGraph$E$Copula.CondData.2[[con[1]]]
    }
    else {
      zr1 <- oldVineGraph$E$Copula.CondData.1[[con[1]]]
    }
    if (temp[2, 1] == same) {
      zr2 <- oldVineGraph$E$Copula.CondData.2[[con[2]]]
    }
    else {
      zr2 <- oldVineGraph$E$Copula.CondData.1[[con[2]]]
    }
    zr1a <- if (is.list(zr1)) 
      as.vector(zr1[[1]])
    else zr1
    zr2a <- if (is.list(zr2)) 
      as.vector(zr2[[1]])
    else zr2
    keine_nas <- !(is.na(zr1a) | is.na(zr2a))
    if (struct.crit == "tau") {
      w <- kdevine:::fasttau(zr1a[keine_nas], zr2a[keine_nas], weights[keine_nas])
    }
    else if (struct.crit == "rho") {
      w <- cov.wt(cbind(zr1a[keine_nas], zr2a[keine_nas]),wt=weights[keine_nas],cor=TRUE)$cor[1, 2]
    }
    else if (struct.crit == "AIC") {
      w <- wtd.kdecop(cbind(zr1a[keine_nas], zr2a[keine_nas]), 
                  info = TRUE, weights=weights[keine_nas])$info$AIC
    }
    else if (struct.crit == "cAIC") {
      w <- wtd.kdecop(cbind(zr1a[keine_nas], zr2a[keine_nas]), 
                  info = TRUE, weights=weights[keine_nas])$info$cAIC
    }
    else if (struct.crit == "hoeffd") {
      w <- abs(hoeffd(cbind(zr1a[keine_nas], zr2a[keine_nas])))
    }
    name.node1 <- strsplit(g$V$names[con[1]], split = " *[,;] *")[[1]]
    name.node2 <- strsplit(g$V$names[con[2]], split = " *[,;] *")[[1]]
    l1 <- c(g$V$conditionedSet[[con[1]]], g$V$conditioningSet[[con[1]]])
    l2 <- c(g$V$conditionedSet[[con[2]]], g$V$conditioningSet[[con[2]]])
    nedSet <- c(setdiff(l1, l2), setdiff(l2, l1))
    ningSet <- intersect(l1, l2)
    nmdiff <- c(setdiff(name.node1, name.node2), setdiff(name.node2, 
                                                         name.node1))
    nmsect <- intersect(name.node1, name.node2)
    name <- paste(paste(nmdiff, collapse = ","), paste(nmsect, 
                                                       collapse = ","), sep = " ; ")
    todel <- FALSE
  }
  list(w = w, nedSet = nedSet, ningSet = ningSet, name = name, 
       todel = todel)
}

##########################################################################
## initializeFirstGraph2 with Spearman rho, kdevine:::initializeFirstGraph2() 
# Searches for first graph based on tree crit
wtd.initializeFirstGraph2 <- function (data.univ, weights, struct.crit = "tau") 
{
  q <- dim(data.univ)[2]
  C <- matrix(rep(1, q * q), ncol = q)
  for (i in 1:(q - 1)) {
    for (j in (i + 1):q) {
      if (struct.crit == "tau") {
        crit <- kdevine:::fasttau(data.univ[, i], data.univ[, 
                                                  j], weights)
      }
      else if (struct.crit == "rho") {
        crit <- cov.wt(data.univ[, c(i, j)],wt=weights,cor=TRUE)$cor[1, 2]
      }
      else if (struct.crit == "AIC") {
        crit <- wtd.kdecop(data.univ[, c(i, j)], info = TRUE, weights=weights)$info$AIC
      }
      else if (struct.crit == "cAIC") {
        crit <- wtd.kdecop(data.univ[, c(i, j)], info = TRUE, weights=weights)$info$cAIC
      }
      else if (struct.crit == "hoeffd") {
        crit <- abs(kdevine:::hoeffd(data.univ[, c(i, j)]))
      }
      C[i, j] <- crit
      C[j, i] <- crit
    }
  }
  rownames(C) <- colnames(C) <- colnames(data.univ)
  rownames(C) <- colnames(C) <- colnames(data.univ)
  kdevine:::graphFromTauMatrix(C)
}


##########################################################################
## buildNextGraph2 with weights and Spearman rho, kdevine:::buildNextGraph2() 
wtd.buildNextGraph2 <- function(oldVineGraph, weights, struct.crit = "tau", parallel) 
{
  d <- nrow(oldVineGraph$E$nums)
  g <- kdevine:::makeFullGraph2(d)
  g$V$names <- oldVineGraph$E$names
  g$V$conditionedSet <- oldVineGraph$E$conditionedSet
  g$V$conditioningSet <- oldVineGraph$E$conditioningSet
  if (parallel) {
    i <- NULL
    out <- foreach(i = 1:nrow(g$E$nums)) %dopar% wtd.getEdgeInfo2(i, 
                                                              g = g, oldVineGraph = oldVineGraph, weights = weights, 
                                                              struct.crit = struct.crit)
  }
  else {
    out <- lapply(1:nrow(g$E$nums), wtd.getEdgeInfo2, g = g, 
                  oldVineGraph = oldVineGraph, weights = weights, 
                  struct.crit = struct.crit)
  }
  g$E$weights <- sapply(out, function(x) x$w)
  g$E$names <- sapply(out, function(x) x$name)
  g$E$conditionedSet <- lapply(out, function(x) x$nedSet)
  g$E$conditioningSet <- lapply(out, function(x) x$ningSet)
  g$E$todel <- sapply(out, function(x) x$todel)
  kdevine:::deleteEdges(g)
}

