##########################################################################
#### 4 - VCdeco - summary functions for VC_deco
##########################################################################

#### The following functions plot the copula density estimates and 
#### estimate the statistics of the marginal distribution

############################################################################
#### Get density plots
get_density_plots <- function(VCopObj, varnames=NULL, colorpalette=NULL){
  
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
    ),
    axis.line=list(col="transparent")
  )

  fit0 <- VCopObj$VCfit$VCfit0
  fit1 <- VCopObj$VCfit$VCfit1
  if(is.null(varnames)){varnames <- fit0$names}
  edges0 <- get_edge_names(fit0)
  edges1 <- get_edge_names(fit1)
  plotlist <- list()
  k=1
  
  # Change colors if required
  if(is.null(colorpalette)){
    jet.colors <- colorRampPalette(c("#00007F", "blue", 
                                                      "#007FFF", "cyan", "#7FFF7F", "yellow", 
                                                      "#FF7F00", "red", "#7F0000"), bias = 2)
    colorpalette  <- jet.colors(100)
  }else{
    colorpalette <- colorRampPalette(colorpalette, bias = 2)
    colorpalette <- colorpalette(100)
  }
  
  
  for(i in 1:length(edges1[[2]])){
   # Create first plot
   bic1 <- get_pair_copula(fit1, tree=edges1[[3]][i], edge=edges1[[4]][i])
   vx <- paste0(varnames[edges1[[6]][[i]][1]], collapse=",")
   vx <- ifelse(is.na(edges1[[7]][[i]][1]),vx,paste0(vx,"|",paste0(varnames[edges1[[7]][[i]]], collapse=",")))
   vy <- paste0(varnames[edges1[[6]][[i]][2]], collapse=",")
   vy <- ifelse(is.na(edges1[[7]][[i]][1]),vy,paste0(vy,"|",paste0(varnames[edges1[[7]][[i]]], collapse=",")))
   pbic1 <- plot(bic1)
   pbic1$panel.args.common$xlab <- list(vx,rot=20)
   pbic1$panel.args.common$ylab <- list(vy,rot=-61)
   pbic1$par.settings <- theme.novpadding
   pbic1$panel.args.common$col.regions <- colorpalette
   
   # Create second plot if available
   sel <- which(edges0[[2]]==edges1[[2]][i])
    if(length(sel)>0){
   bic0 <- get_pair_copula(fit0, tree=edges0[[3]][sel], edge=edges0[[4]][sel])
   if(edges0[[1]][sel]!=edges1[[1]][i]){
      bic0$parameters <- t(bic0$parameters) 
   }
   pbic0 <- plot(bic0)
   pbic0$panel.args.common$xlab <- list(vx,rot=20)
   pbic0$panel.args.common$ylab <- list(vy,rot=-61)
   pbic0$par.settings <- theme.novpadding
   pbic0$panel.args.common$col.regions <- colorpalette
   }else{
   pbic0 <- NA
   }
   plotlist[[k]] <- list(pbic0,pbic1)
   k=k+1
  }
  return(plotlist)
}

############################################################################
# Comparing copula density plots of two groups
compare_plots <- function(plotlist, edge=1, diff=FALSE,
                          groups=c(0,1), colorpalette=NULL,
                          bias=2, zlim=NULL, title=TRUE){
  
  #t0 <- textGrob(as.character(groups[1]), gp=gpar(fontsize=14, col="black", fontface="bold"))
  #t1 <- textGrob(as.character(groups[2]), gp=gpar(fontsize=14, col="black", fontface="bold"))
  #tdiff <- textGrob("Diff", gp=gpar(fontsize=14, col="black", fontface="bold"))
  
  # Change colors if required
  if(is.null(colorpalette)){
    jet.colors <- colorRampPalette(c("#00007F", "blue", 
                                     "#007FFF", "cyan", "#7FFF7F", "yellow", 
                                     "#FF7F00", "red", "#7F0000"), bias = 2)
    colorpalette  <- jet.colors(100)
  }else{
    colorpalette <- colorRampPalette(colorpalette, bias = bias)
    colorpalette <- colorpalette(100)
  }
  
  
  if(diff){  
    cdplots <- grid.arrange(t0,t1, tdiff,
                            plotlist[[edge]], plotlist[[edge+6]],plotlist[[edge+12]],
                            ncol=3, nrow=2,
                            padding=0,
                            heights=c(0.1,1),#unit(c(0.5,5,5,5),rep("cm",4)),
                            widths=c(1,1,1))
  }else{  
    # Harominze col.regions and zlim
    p1 <- plotlist[[edge]][[1]]
    p2 <- plotlist[[edge]][[2]]
    zlim1 <- p1$panel.args.common$zlim
    zlim2 <- p2$panel.args.common$zlim
    if(is.null(zlim)){
    zlim <- max(c(zlim1,zlim2))
    }
    p1$panel.args.common$zlim[2] <- zlim
    p2$panel.args.common$zlim[2] <- zlim
    at1 <- p1$panel.args.common$at
    at2 <- p2$panel.args.common$at
    r <- range(c(at1,at2))
    at <- seq(from=r[1],to=r[2],by=(r[2]-r[1])/(length(at1)-1))
    p1$panel.args.common$at <- at
    p2$panel.args.common$at <- at
    p1$panel.args.common$col.regions <- colorpalette
    p2$panel.args.common$col.regions <- colorpalette
    if(title){
    p1$main <- as.character(groups[1])
    p2$main <- as.character(groups[2])
    }
    cdplots <- grid.arrange(#t0,t1,
                            p1,p2, 
                            ncol=2, nrow=1,
                            padding=0,
                            heights=c(1.5),#unit(c(0.5,5,5,5),rep("cm",4)),
                            widths=c(1,1)#unit(c(5,5),rep("cm",2))
    )
  }
  
  return(cdplots)
}

############################################################################
#### Function to plot observed and fitted quantiles 

VC_deco_plot_obs_vs_fit <- function(VCopObj){
  
  plotdf <- VCopObj$VCdeco$quants_sim_vs_obs
  plotdf <- melt(plotdf, id.vars=names(plotdf)[1], 
                 measure.vars=names(plotdf)[-1],
                 variable.name = "group",
                 value.name = "delta")
  plotdf$tau <- as.numeric(as.character(plotdf$tau))
  
  if(is.null(VCopObj$VCdeco_se)==FALSE){
    #Standard error
    plotdfse <- VCopObj$VCdeco_se$quants_sim_vs_obs_se
    plotdfse <- melt(plotdfse, id.vars=names(plotdfse)[1], 
                     measure.vars=names(plotdfse)[-1],
                     variable.name = "group",
                     value.name = "delta")
    plotdfse$tau <- as.numeric(as.character(plotdfse$tau))
    plotdf$se <- plotdfse$delta
    
    #c-value for uniform ci
    plotdf[which(plotdf$group=="group1"),"t"] <- VCopObj$VCdeco_se$quants_sim_vs_obs_critval[1,2]
    plotdf[which(plotdf$group=="group0"),"t"] <- VCopObj$VCdeco_se$quants_sim_vs_obs_critval[2,2]
    plotdf[which(plotdf$group=="DiffInDiff"),"t"] <- VCopObj$VCdeco_se$quants_sim_vs_obs_critval[3,2]
    
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

############################################################################
#### deco_qqplot
VC_deco_qqplot <- function(VCopObj, formula=NULL, group=c("0","1"),
                           logscale=TRUE, type=1:5,
                           tau=1:99/100){
  if(is.null(formula)){
    y0 <- VCopObj$VCfit$data0[,1]
    y1 <- VCopObj$VCfit$data1[,1]
    ysim0 <- VCopObj$VCdeco$aggregate$ysim[,2]
    ysim1 <- VCopObj$VCdeco$aggregate$ysim[,1]
    w0 <- VCopObj$VCfit$data0[,"(weights)"]
    w1 <- VCopObj$VCfit$data1[,"(weights)"]
  }else{
    y0 <- as.data.frame(do.call("cbind",as.list(VCopObj$VCfit$data0[,1:ncol(VCopObj$VCfit$data0)])))
    y1 <- as.data.frame(do.call("cbind",as.list(VCopObj$VCfit$data1[,1:ncol(VCopObj$VCfit$data1)])))
    w0 <- model.frame(formula,y0, weights=`(weights)`)[,"(weights)"]
    w1 <- model.frame(formula,y1, weights=`(weights)`)[,"(weights)"]
    y0 <- rowSums(model.matrix(formula,y0))-1
    y1 <- rowSums(model.matrix(formula,y1))-1
    ysim0 <- rowSums(model.matrix(formula,VCopObj$VCdeco$aggregate$xsim[[2]]))-1
    ysim1 <- rowSums(model.matrix(formula,VCopObj$VCdeco$aggregate$xsim[[1]]))-1
  }
  
  #ysim0 <- VCopObj$VCdeco$xsim[[2]]
  #sum(y0[which(y0$ANNUALHOURS_m>0&y0$ANNUALHOURS_f>0),"(weights)"])/sum(y0[,"(weights)"])
  #sum(ysim0$ANNUALHOURS_m>0&ysim0$ANNUALHOURS_f>0)/nrow(ysim0)
  
  #ysim1 <- VCopObj$VCdeco$xsim[[1]]
  #sum(y1[which(y1$ANNUALHOURS_m>0&y1$ANNUALHOURS_f>0),"(weights)"])/sum(y1[,"(weights)"])
  #sum(ysim1$ANNUALHOURS_m>0&ysim1$ANNUALHOURS_f>0)/nrow(ysim1)
  
  #Observed
  q0 <- Hmisc::wtd.quantile(y0,w0, probs=tau) 
  q1 <- Hmisc::wtd.quantile(y1,w1, probs=tau) 
  qsim0 <- quantile(ysim0,probs=tau)
  qsim1 <- quantile(ysim1,probs=tau)
  m0 <- Hmisc::wtd.mean(y0,w0) 
  m1 <- Hmisc::wtd.mean(y1,w1) 
  s0 <- sqrt(Hmisc::wtd.var(y0,w0)) 
  s1 <-  sqrt(Hmisc::wtd.var(y1,w1)) 
  msim0 <- mean(ysim0) 
  msim1 <- mean(ysim1) 
  ssim0 <- sqrt(var(ysim0)) 
  ssim1 <-  sqrt(var(ysim1))
  
  
  dfplot <- data.frame(group=rep(group,each=length(tau)),
                       tau=rep(tau,2),
                       observed=c(q0,q1),
                       simulated=c(qsim0,qsim1))
  dfplot$observed_std <- (dfplot$observed-rep(c(m0,m1), each=length(tau)))/rep(c(s0,s1), each=length(tau))
  dfplot$simulated_std <- (dfplot$simulated-rep(c(msim0,msim1), each=length(tau)))/rep(c(ssim0,ssim1), each=length(tau))
  
  
  if(type==1){
  dfplot$observed <- ifelse(dfplot$observed==0&logscale, 1, dfplot$observed)  
  dfplot$simulated <- ifelse(dfplot$simulated==0&logscale, 1, dfplot$simulated)  
  plot <- ggplot(dfplot, aes(observed,simulated, col=group)) + 
    geom_point(shape=1) + ggtitle("Q-Q plot") + 
    geom_abline(intercept=0,slope=1, linetype = "dashed", color="darkgrey") + 
    facet_wrap(~group) + 
    scale_x_continuous(trans = ifelse(logscale,'log10','identity'), labels = scales::number_format(big.mark=",")) +
    scale_y_continuous(trans = ifelse(logscale,'log10','identity'), labels = scales::number_format(big.mark=",")) + 
    theme(legend.position = "none") 
  }else if(type==2){
    plot <- ggplot(dfplot, aes(observed_std,simulated_std, col=group)) + 
      geom_point(shape=1) + ggtitle("Q-Q plot") + 
      geom_abline(intercept=0,slope=1, linetype = "dashed", color="darkgrey") + 
      facet_wrap(~group) + 
      theme(legend.position = "none") 
 }else if(type==3){
    dfplot <- melt(dfplot, id.vars=names(dfplot)[1:2],
                variable.name = "series",
                value.name = "quantile")
    plot <- ggplot(dfplot,aes(tau,quantile, color=series)) + 
            geom_point() + geom_line() + 
            ggtitle("Quantiles") + 
            facet_wrap(~group)
    }else{
     if(logscale){
       dfplot$diff <-  log(dfplot$simulated) - log(dfplot$observed)
     }else{
       dfplot$diff <-  dfplot$simulated - dfplot$observed
     }
     plot <- ggplot(dfplot,aes(tau,diff, color=group)) + 
       geom_point() + geom_smooth(formula= y~x, method=loess, span=0.5) + 
       ggtitle("Simulated minus observed quantiles") + 
       facet_wrap(~group) + geom_hline(yintercept=0, color="darkgrey")
  }
  plot
}

############################################################################
#### create stat table
create_stats_table <- function(VCopObj,
                               stats=c("inequality","distribution"),
                               type=c("aggregate","marginal","copula"),
                               varnames=NULL,
                               tex=FALSE,
                               caption=NULL,
                               col.width=1){
  
  d <- VCopObj$VCfit$VCfit0$structure$d
  deco <- VCopObj$VCdeco$aggregate$stats_deco
  
  deco$names2 <- dplyr::recode(deco$names2,
                               "1 Total difference"="1 Observed difference",
                               "2 Marginals"="2 Aggregate marginal effect",
                               "3 Copula"="3 Aggregate dependence effect")

  # Rename if required
  if(length(varnames)==d){
    original_names <- as.character(deco[which(is.element(deco$names1,paste0("M ",1:d))),"names2"])
    original_names <- substr(original_names,3,nchar(original_names))
    for(i in 1:d){
      deco$names2 <- gsub(original_names[i],varnames[i],deco$names2)
    }
  }
  if(tex){
    deco$names2 <- paste0("\\quad ",deco$names2)
  }
  
  statistics <- switch(stats[1],
                       inequality=c("Std. Dev.","Gini", "P90/P10", "P90/P50", "P50/P10"),
                       distribution=c("Mean","P10", "P50", "P90", "\\%obs.=0")
                       )
  digits_statistics <- switch(stats[1],
                              inequality= c(0,2,3,2,2,2),
                              distribution=c(0,2,2,2,2,3))
  seleffect <- switch(type[1],
                      aggregate=which(is.element(deco$cat2,c("aggregate","direct"))),
                      marginal=which(is.element(deco$cat1,c("Marginal")))[-c(1,sum(is.element(deco$cat1,c("Marginal"))))],
                      copula=which(is.element(deco$cat1,c("Copula")))[-c(1,sum(is.element(deco$cat1,c("Copula"))))])
  subtitle_pos <- switch(type[1],
                        aggregate=c(3,3+d,3+d+sum(1:(d-1)))+0:2,
                        marginal=d,
                        copula=sum(1:(d-1)))
  if(tex){
  subtitle <- switch(type[1],
                     aggregate=c("\\textbf{Direct marginal effects}","\\textbf{Direct dependence effects}","\\textbf{Interactions}"),
                     marginal=c("\\textbf{'Two-way' interactions}"),
                     copula=c("\\textbf{'Two-way' interactions}"))
  }else{
  subtitle <- switch(type[1],
                     aggregate=c("Direct marginal effects","Direct dependence effects","Interactions"),
                     marginal=c("'Two-way' interactions"),
                     copula=c("'Two-way' interactions"))
  }
  
  
  tab <- deco[seleffect,statistics]
  rownames(tab) <- deco[seleffect,"names2"]
  colnames(tab) <- statistics
  add <- ifelse(tex,"\\quad","")
  if(type[1]=="aggregate"&nrow(tab)>3){
    rownames(tab)[(nrow(tab)-1):nrow(tab)] <- paste0(add,c("1 Marginals", "2 Dependence"))
  }else if(nrow(tab)>3){
    rownames(tab)[nrow(tab)] <- paste0(add,"Higher order interactions")
  }
  #Round results
  cn <- colnames(tab)
  rn <- rownames(tab)
  tab <- do.call("cbind",lapply(1:ncol(tab), function(i) round(tab[,i], digits_statistics[-1][i])))
  colnames(tab) <- cn
  rownames(tab) <- rn
  
  # Add standard errors if availale
  if(is.null(VCopObj$VCdeco_se)==FALSE){
    tabse <- VCopObj$VCdeco_se$aggregate$stats_deco_se
    tabse <- tabse[match(deco$names1,tabse$names1),]
    tabse <- tabse[seleffect,statistics]
    # Round s.e.
    cn <- colnames(tabse)
    rn <- rownames(tabse)
    tabse <- do.call("cbind",lapply(1:ncol(tabse), function(i) round(tabse[,i], digits_statistics[-1][i])))
    colnames(tabse) <- cn
    rownames(tabse) <- rn
    tab0 <- tab
    tab <- NULL
    for(i in 1:ncol(tab0)){
      tab <- cbind(tab,tab0[,i],tabse[,i])
    }
    cn <-rep(colnames(tab0),each=2)
    cn[seq(2,length(cn),2)] <- unlist(lapply(mapply(rep, times=1:(length(cn)/2), x=rep(" ",length(cn)/2)),function(x) paste0(x,collapse="")))
    colnames(tab) <- cn 
    rownames(tab) <- rownames(tab0)
    rm("tab0")
    
  }
  
  # Add subtitle
  for(i in 1:length(subtitle)){
    tab <- rbind(tab[1:subtitle_pos[i],],
                 tab[subtitle_pos[i],],
                 tab[(subtitle_pos[i]+1):nrow(tab),])
    tab[subtitle_pos[i]+1,] <- c(rep(NA,ncol(tab)))
    rownames(tab)[subtitle_pos[i]+1] <- subtitle[i]
  }
  
  # Add brackets to se
  if(is.null(VCopObj$VCdeco_se)==FALSE){
    for(i in seq(2,ncol(tab),2)){
      tab[,i] <- paste0("(",tab[,i],")")
    }
    for(i in 1:nrow(tab)){
      if(is.na(tab[i,1])){tab[i,] <- rep(NA, ncol(tab))} 
    }
  }
  
  if(tex){
  tab <- xtable(tab,
                caption=caption,
                type="latex",
                booktabs = FALSE,
                auto=TRUE)
  #digits(tab) <-  matrix(rep(digits_statistics,nrow(tab)),nrow(tab),ncol(tab)+1,byrow=TRUE)
  align(tab) <- c("X",rep(paste0("R{",col.width,"cm}"),ncol(tab)))
  return(tab)
  }else{
  print(as.data.frame(tab))
}
}
############################################################################
## slightly adapted function getAnywhere("plot.vinecop")
require(gridExtra)
vinegraph <- function (x, tree = 1, var_names = "ignore", edge_labels = NULL, 
                       ...) 
{
  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop("The 'ggraph' package must be installed to plot.")
  }
  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("The 'grid' package must be installed to plot.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The 'ggplot2' package must be installed to plot.")
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("The 'igraph' package must be installed to plot.")
  }
  d <- dim(x)[1]
  trunc_lvl <- dim(x)[2]
  if (!identical(tree, "ALL") && any(tree > trunc_lvl)) {
    stop("Selected tree does not exist.")
  }
  if (any(tree == "ALL")) {
    if (d > 5) {
      warning(paste("tree = 'ALL' is not recommended for d > 5 and", 
                    " it is set as c(1,2), please use tree = 1:d"))
      tree <- c(1, 2)
    }
    else {
      tree <- seq_len(trunc_lvl)
    }
  }
  assertthat::assert_that(rvinecopulib:::in_set(var_names, c("ignore", "use", 
                                  "legend")))
  if (!is.null(edge_labels)) {
    assertthat::assert_that(rvinecopulib:::in_set(edge_labels, c("pair", "tau", 
                                      "family", "family_tau")))
  }
  if (is.null(x$names)) {
    x$names <- as.character(1:d)
  }
  if (var_names %in% c("ignore", "legend")) {
    names <- x$names
    x$names <- as.character(1:d)
  }
  g <- lapply(tree, rvinecopulib:::get_graph, vc = x, edge_labels = edge_labels)
                #,var_names = var_names)
  plots <- vector("list", length(tree))
  name <- NULL
  for (i in seq_along(tree)) {
    p <- ggraph::ggraph(g[[i]], "igraph", algorithm = "tree", 
                        circular = TRUE)
    if (!is.null(edge_labels)) {
      p <- p + ggraph::geom_edge_link(ggplot2::aes(label = name), 
                                      colour = "darkgrey",#"#000000",
                                      angle_calc = "along", 
                                      label_dodge = grid::unit(7, "points"))
    }
    else {
      p <- p + ggraph::geom_edge_link(colour ="darkgrey",#"#000000",
                                      )
    }
    p <- p + ggraph::geom_node_point(col = "#EB811B", # "#56B4E9", 
                                     size = 3) + ggraph::geom_node_text(ggplot2::aes(label = name), 
                                                                        #fontface = "bold", 
                                                                        repel = TRUE) + ggplot2::theme_void() + 
      ggplot2::labs(title = paste0("Tree ", tree[i]),
                    fontface = "bold")
    if (var_names == "legend") {
      p <- p + ggplot2::labs(caption = paste(x$names, names, 
                                             sep = " = ", collapse = ", "))
    }
    plots[[i]] <- p + ggplot2::theme(plot.margin = ggplot2::margin(1, 
                                                                   1, 1, 1, "pt"))
  }
  if (length(tree) > 3) {
    return(grid.arrange(grobs=plots, nrow=ceiling(length(plots)/2), ncol=2))
    #invisible(multiplot(plotlist = plots, cols = 2))
  }
  else {
    return(grid.arrange(grobs=plots, nrow=length(plots), ncol=1))
    #return(invisible(multiplot(plotlist = plots)))
  }
}