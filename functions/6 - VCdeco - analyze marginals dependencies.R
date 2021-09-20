##########################################################################
#### 6 - VCdeco - analyze marginals and dependencies
##########################################################################

#### The following functions plot the copula density estimates and 
#### estimate the statistics of the marginal distribution

############################################################################
#### Compute descriptive statistics of marginal distributions

VC_stats_marginals <- function(formula,data,group,weights,
                               tau=seq(5,95,5)/100,
                               reference=NULL,
                               kendall=TRUE,          # Estimate Kendall's tau?
                               bootstrap=FALSE,       # boostrap 
                               it=100,                # boostrap iterations
                               psu=NULL,              # psu 
                               stratum=NULL,          # stratum   
                               alpha=0.05,            # significance for uniform ci
                               ncore=1,               # number of core for bootstrap
                               trimgini=0,            # trim gini?
                               trimcor=0              # correlation coef?
                               ){

  #Use match.call to create data.frame
  mf = match.call()
  m = match(c("formula", "data", "weights", "na.action","group"), names(mf), 0)
  mf = mf[c(1, m)]
  mf$drop.unused.levels = TRUE
  mf[[1]] = as.name("model.frame")
  
  #Add psu and stratum to formula in order to integrate it to model frame
  if(is.null(psu)==TRUE&is.null(stratum)==FALSE){
    mf$formula <- as.formula(paste0(paste0(as.character(formula)[c(2,1,3)], collapse=" ")," +",stratum)) 
  }else if(is.null(psu)==FALSE&is.null(stratum)==TRUE){
    mf$formula <- as.formula(paste0(paste0(as.character(formula)[c(2,1,3)], collapse=" ")," +",psu))   
  }else if(is.null(psu)==FALSE&is.null(stratum)==FALSE){
    mf$formula <- as.formula(paste0(paste0(as.character(formula)[c(2,1,3)], collapse=" ")," +",psu," +", stratum))   
  }
  
  # Evaluate model.frame
  mf = eval.parent(mf)
  mt = attr(mf, "terms")
 
  #Retrieve weights
  if(is.null(model.weights(mf))){
    mf$`(weights)` <- rep(1,nrow(mf))
  }
  
  # Set reference group 
  if(is.null(reference)==FALSE){
    sel0 <- reference
  } else if(is.factor(mf[,"(group)"])){
    sel0 <- levels(mf[,"(group)"])[1]
  }else{
    sel0 <- sort(unique(mf[,"(group)"]))[1]
  }
  
  #drop observations with zero weight
  cat(paste0(c("\nDropped because of 0 household weight:", length(which(mf[,"(weights)"]==0&mf[,"(group)"]==sel0)),"in group 0 and",
               length(which(mf[,"(weights)"]==0&mf[,"(group)"]!=sel0)) ,"in group 1.\n")))
  mf <- mf[which(mf[,"(weights)"]!=0),]
  
  # Select data sets
  data0 <- mf[which(mf[,"(group)"]==sel0),]
  data1 <- mf[which(mf[,"(group)"]!=sel0),]

  # Estimate statistics
  res <- stats_marginals_est(data0,data1,tau,kendall,trimgini,trimcor)

  # Bootstrap se 
  if(bootstrap){
    cat("\nBootstrapping...\n")
    se <- stats_marginals_bs(res,it=it, ncore=ncore,
                             psu=psu, stratum=stratum,        
                             alpha=alpha)
  }else{
    se <- NULL  
  }
  
  return(list(estimates=res,
              se=se))
}

############################################################################
#### Function of estimating descriptive statistics of marginals

stats_marginals_est <- function(data0, data1, tau, kendall,trimgini,trimcor){
  
  # Get variable names
  vn <- names(data0)[1:5]
  vA <- NULL
  vB <- NULL
  for(i in 1:3){
    vA <- c(vA,rep(i,4-i))
    for(j in (i+1):4){
      vB <- c(vB,j) 
    }
  }
  vcomb <- paste0(vn[vA+1],rep("-",6),vn[vB+1])
  
  # Compute weighted percentile ranks
  u0 <- apply(data0[,vn[-1]],2, function(x) rank.wt(x,data0[,"(weights)"]))
  u1 <- apply(data1[,vn[-1]],2, function(x) rank.wt(x,data1[,"(weights)"]))
  
  # Quantiles
  quants0 <- as.data.frame(sapply(vn, function(x) wtd.quantile(data0[,x],weights=data0[,"(weights)"], probs=tau)))
  quants1 <- as.data.frame(sapply(vn, function(x) wtd.quantile(data1[,x],weights=data1[,"(weights)"], probs=tau)))
  
  quants0 <- cbind(data.frame(group=0,tau=tau),quants0)
  quants1 <- cbind(data.frame(group=1,tau=tau),quants1)
  
  # Summary statistics
  stats0 <- as.data.frame(sapply(vn, function(x) wtd.sum(data0[,x],data0[,"(weights)"],trim=trimgini)))
  stats1 <- as.data.frame(sapply(vn, function(x) wtd.sum(data1[,x],data1[,"(weights)"],trim=trimgini)))
  
  stats0 <- cbind(data.frame(group=0,stat=rownames(stats0)),stats0)
  stats1 <- cbind(data.frame(group=1,stat=rownames(stats1)),stats1)
  
  ## Correlations
  s0 <- unique (unlist (lapply (data0[,vn[-1]], function (x) which (is.na (x)))))
  s0 <- setdiff(1:nrow(data0),s0)
  s1 <- unique (unlist (lapply (data1[,vn[-1]], function (x) which (is.na (x)))))
  s1 <- setdiff(1:nrow(data1),s1)
  
  # Trim data set if requested for Pearson rho
  if(length(trimcor)==1&trimcor[1]==0){
    sel0 <- 1:nrow(data0[s0,])
    sel1 <- 1:nrow(data1[s1,])
  }else{
    sel0 <- trims(data0[s0,vn[-1]], data0[s0,"(weights)"], trimcor)
    sel1 <- trims(data1[s1,vn[-1]], data1[s1,"(weights)"], trimcor)
  }
  
  rho_pearson0 <- cov.wt(data0[s0,vn[-1]][sel0,], wt=data0[s0,"(weights)"][sel0], cor=TRUE)$cor
  rho_pearson1 <- cov.wt(data1[s1,vn[-1]][sel1,], wt=data1[s1,"(weights)"][sel1], cor=TRUE)$cor
  
  rho_spearman0 <- cov.wt(u0[s0,], wt=data0[s0,"(weights)"], cor=TRUE)$cor
  rho_spearman1 <- cov.wt(u1[s1,], wt=data1[s1,"(weights)"], cor=TRUE)$cor

  if(kendall){
   tau_kendalls0 <- VineCopula::TauMatrix(u0[s0,], weights=data0[s0,"(weights)"])
   tau_kendalls1 <- VineCopula::TauMatrix(u1[s1,], weights=data1[s1,"(weights)"])
  }else{
   tau_kendalls0 <- matrix(NA,4,4)
   tau_kendalls1 <- matrix(NA,4,4)
  }
  
  corrs0 <-  data.frame(group=0,
                        vars=rep(vcomb,3),
                        rho=c(rho_pearson0[lower.tri(rho_pearson0)],
                              rho_spearman0[lower.tri(rho_spearman0)],
                              tau_kendalls0[lower.tri(tau_kendalls0)]),
                        type=rep(c("cor.","rank cor.","tau"), each=length(vcomb)))
  
  #Order dependencies by Kendall's tau in decreasing order  
  tau_kendalls0 <- tau_kendalls0[lower.tri(tau_kendalls0)]
  #cat(" \n")
  #cat("Dependencies group 0:\n")
  #cat(paste0(paste0(1:6," ",vcomb[order(abs(tau_kendalls0),decreasing = TRUE)]," (",round(tau_kendalls0[order(abs(tau_kendalls0), decreasing = TRUE)],2),")"),collapse="\n"))
  
  corrs1 <-  data.frame(group=1,
                        vars=rep(vcomb,3),
                        rho=c(rho_pearson1[lower.tri(rho_pearson1)],
                              rho_spearman1[lower.tri(rho_spearman1)],
                              tau_kendalls1[lower.tri(tau_kendalls1)]),
                        type=rep(c("cor.","rank cor.","tau"), each=length(vcomb)))
  
  # Order dependencies by Kendall's tau in dreasing order  
  tau_kendalls1 <- tau_kendalls1[lower.tri(tau_kendalls1)]
  
  #cat(" \n")
  #cat("Dependencies group 1:\n")
  #cat(paste0(paste0(1:6," ",vcomb[order(abs(tau_kendalls1),decreasing = TRUE)]," (",round(tau_kendalls1[order(abs(tau_kendalls1), decreasing = TRUE)],2),")"),collapse="\n"))
  
  # Prepare results to export
  quants <- rbind(quants0,quants1)
  stats <- rbind(stats0,stats1)
  corrs <- rbind(corrs0,corrs1)
  
  return(list(quants=quants,
              stats=stats,
              corrs=corrs,
              data0=data0,
              data1=data1,
              varnames=vn,
              tau=tau))
  
}


############################################################################
#### Get names of variables in right order

VC_arrange_var_names <- function(VCopObj, varnames=NULL){
  
  if(is.null(varnames)){
  varnames <- VCopObj$varnames
  }
  
  order <- VCopObj$order
  type <- VCopObj$type
  ovn <- varnames[order]
  
  #D-Vine
  if(type=="DVine"){
  vx <- ovn[c(1,2,3,1,2,1)]
  vy <- ovn[c(2,3,4,3,4,4)]
  vcond <- c(rep("",3),ovn[2],ovn[3],paste0(ovn[2],",",ovn[3]))
  }else{
  vx <- ovn[c(1,1,1,2,2,3)]
  vy <- ovn[c(2,3,4,3,4,4)]
  vcond <- c(rep("",3),ovn[1],ovn[1],paste0(ovn[1],",",ovn[2]))  
  }
  
  # conditioning variables and separators
  seps <- c(rep("",3),rep("|",3))
  
  # Add together
  vcopula <- paste0(vx,rep(",",6),vy,seps,vcond) # Names of edges/copulas
  vx <- paste0(vx,seps,vcond)
  vy <- paste0(vy,seps,vcond)
  
  return(list(vx=vx,vy=vy,vcopula=vcopula))
  
}

############################################################################
#### Collect copula density plots

VC_deco_plot_density <- function(VCopObj, varnames=NULL){
  
  VCopObj <- VCopObj$VCfit
  
  #Arrange variable names
  vn <- VC_arrange_var_names(VCopObj, varnames=varnames)
  vx <- vn$vx
  vy <- vn$vy
  
  # Change padding of plots
  theme.novpadding <- list(
    layout.heights = list(
      top.padding = 0,
      main.key.padding = 0,
      key.axis.padding = 0,
      axis.xlab.padding = 0,
      xlab.key.padding = 0,
      key.sub.padding = 0,
      bottom.padding = 0
    ),
    layout.widths = list(
      left.padding = 0,
      key.ylab.padding = 0,
      ylab.axis.padding = 0,
      axis.key.padding = 0,
      right.padding = 0
    )
  )
  
  # Plotting copulas
  plotlist <- NULL
  k=0
  for(j in c("copula0","copula1")){
    copulas <- VCopObj[[j]][[1]]
    for(i in 1:6){
      copulaplot <- plot(copulas[[i]],
                         #size=25,
                         par.settings = theme.novpadding,
                         zoom=0.75) 
      copulaplot$panel.args.common$xlab[1] <- vx[i]
      copulaplot$panel.args.common$ylab[1] <- vy[i]
      
      #Save to plot list
      plotlist[[i+k]] <- copulaplot 
      names(plotlist)[i+k] <- paste0(j,"-",vx[i],"-",vy[i])
    }
    k=k+6
  }
  for(i in 1:6){
    copulaplot <- VC_deco_plot_density_diff(VCopObj[["copula0"]][[1]][[i]],VCopObj[["copula1"]][[1]][[i]]) 
    copulaplot$par.settings <- theme.novpadding
    copulaplot$zoom <- 0.75
    copulaplot$panel.args.common$xlab[1] <- vx[i]
    copulaplot$panel.args.common$ylab[1] <- vy[i]
    
    #Save to plot list
    plotlist[[i+k]] <- copulaplot 
    names(plotlist)[i+k] <- paste0("diff-",vx[i],"-",vy[i])
  }
  return(plotlist)
}

############################################################################
#### Collect copula density plots

# Comparing copula density plots of two groups
VC_deco_plot_compare <- function(plotlist, edge=1, diff=FALSE,
                                 groups=c(0,1), conditional=c(0,1)){
  
  t0 <- textGrob(as.character(groups[1]), gp=gpar(fontsize=20, col="black", fontface="bold"))
  t1 <- textGrob(as.character(groups[2]), gp=gpar(fontsize=20, col="black", fontface="bold"))
  tdiff <- textGrob("Diff", gp=gpar(fontsize=20, col="black", fontface="bold"))
  
  if(diff){  
  cdplots <- grid.arrange(t0,t1, tdiff,
                         plotlist[[edge]], plotlist[[edge+6]],plotlist[[edge+12]],
                         ncol=3, nrow=2,
                         padding=0,
                         heights=c(0.1,1),#unit(c(0.5,5,5,5),rep("cm",4)),
                         widths=c(1,1,1))
  }else{  
  sel <- 1:3+(conditional[1]*3)
  
  plotlist0 <- plotlist[1:6]
  plotlist1 <- plotlist[7:12]
  
  cdplots <- grid.arrange(t0,t1,
                          plotlist0[[sel[1]]], plotlist1[[sel[1]]],
                          plotlist0[[sel[2]]], plotlist1[[sel[2]]],
                          plotlist0[[sel[3]]], plotlist1[[sel[3]]],
                          ncol=2, nrow=4,
                          padding=0,
                          heights=c(0.1,1,1,1),#unit(c(0.5,5,5,5),rep("cm",4)),
                          widths=c(1,1)#unit(c(5,5),rep("cm",2))
  )
  }
  
  return(cdplots)
}

############################################################################
#### Function to plot differences between densities

# Based on code by Naglers plot.kdecopula() function
VC_deco_plot_density_diff <- function(copula0, copula1, size=20){
  
  xylim <- c(0.01, 1 - 0.01)
  points <- 1:size/(size + 1)
  g <- as.matrix(expand.grid(points, points))
  points <- g[1L:size, 1L]
  adj <- 1
  gu <- g[, 1L]
  gv <- g[, 2L]
  levels <- c(0.2, 0.6, 1, 1.5, 2, 3, 5, 10, 20)
  xlim <- ylim <- xylim
  at <- c(seq(0, 3, length.out = 50), seq(5, 100, length.out = 50))
  at <- c(-at[100:1],at)
  
  # Density difference
  vals <- dkdecop(g, copula1, stable = FALSE) - dkdecop(g, copula0, stable = FALSE)
  cop <- matrix(vals, size, size)
  nms <- colnames(copula0$udata)
  xlab <- nms[1]
  ylab <- nms[2]
  # Plot details
  lst <- list(u = gu, v = gv, c = as.vector(cop) * as.vector(adj))
  TUMblue <- rgb(0, 103/255, 198/255)
  wireframe(x = c ~ u * v, data = lst, scales = list(arrows = FALSE),
            drape = TRUE, colorkey = FALSE, screen = list(z = 25, x = -55),
            shade = FALSE, aspect = c(1, 1), light.source = c(10, 0, 10),
            zoom = 0.85, par.settings = list(axis.line = list(col = "transparent")),
           # at = at,
            col.regions = c(colorRampPalette(c(kdecopula:::tint(TUMblue, 0.5), "white"))(50),
                            rep("white", 50)), 
            xlab = list(xlab, rot = 20), ylab = list(ylab,  rot = 30 - 90), 
            zlab = "", zlim = c(1.1 * min(lst$c), 1.1 * max(lst$c)))
}

############################################################################
#### Function to plot observed and fitted quantiles 

VC_deco_plot_obs_vs_fit <- function(VCopObj, decomposition="aggregate"){
  
  plotdf <- VCopObj$VCdeco[[decomposition]]$quants_sim_vs_obs
  plotdf <- melt(plotdf, id.vars=names(plotdf)[1], 
                 measure.vars=names(plotdf)[-1],
                 variable.name = "group",
                 value.name = "delta")
  plotdf$tau <- as.numeric(as.character(plotdf$tau))

  if(is.null(VCopObj$VCdeco_se$aggregate)==FALSE){
    #Standard error
    plotdfse <- VCopObj$VCdeco_se[[decomposition]]$quants_sim_vs_obs_se
    plotdfse <- melt(plotdfse, id.vars=names(plotdfse)[1], 
                     measure.vars=names(plotdfse)[-1],
                     variable.name = "group",
                     value.name = "delta")
    plotdfse$tau <- as.numeric(as.character(plotdfse$tau))
    plotdf$se <- plotdfse$delta
    
    #c-value for uniform ci
    plotdf[which(plotdf$group=="group1"),"t"] <- VCopObj$VCdeco_se[[decomposition]]$quants_sim_vs_obs_critval[1,2]
    plotdf[which(plotdf$group=="group0"),"t"] <- VCopObj$VCdeco_se[[decomposition]]$quants_sim_vs_obs_critval[2,2]
    plotdf[which(plotdf$group=="DiffInDiff"),"t"] <- VCopObj$VCdeco_se[[decomposition]]$quants_sim_vs_obs_critval[3,2]
     
  }else{
    plotdf$se <- 0 
    plotdf$t <- 0 
  }
  
  plot <- ggplot(plotdf, aes(tau, delta, color=group, fill=group)) +  
            geom_line() + geom_point() + 
            geom_hline(yintercept=0, col="lightgrey") + 
            geom_ribbon(aes(ymin=delta - 1.96*se,ymax=delta + 1.96*se), alpha=0.3) +
            geom_ribbon(aes(ymin=delta - t*se,ymax=delta + t*se), alpha=0.1) +
            facet_wrap(~group, ncol=3)
  return(plot)
}



