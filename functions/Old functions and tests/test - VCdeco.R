#############################################
### Test VCdeco with simulated data



#Simulate data 
n=10000

#bi0 <- bicop_dist(family = "gumbel",rotation = 0,  parameters = 1.2,  var_types = c("c", "c"))
bi0 <- bicop_dist(family = "gauss",rotation = 0,  parameters = 0.2,  var_types = c("c", "c"))
bi1 <- bicop_dist(family = "gauss",rotation = 0,  parameters = 0.5,  var_types = c("c", "c"))

plot(bi0)
plot(bi1)

set.seed(123)
ru0 <- rbicop(n,bi0)
ru1 <- rbicop(n,bi1)

curve(dlnorm(x,3,0.8),1,100)
curve(dlnorm(x,3.5,0.9),1,100, add=TRUE, col="red")

summary(qlnorm(ru0[,2],3,0.8))
summary(qlnorm(ru1[,2],3,0.9))

data0 <- data.frame(d=qbinom(ru0[,1],1,0.7),
                    d2=qbinom(ru0[,2],1,0.3),
                    h=qpois(ru0[,1],30),
                    x=qlnorm(ru0[,2],3,0.8))
data0$y <- data0$h*data0$x
data0$y1 <- data0$d*data0$x
data0$y2 <- data0$d2*data0$d
data0$group <- 0
data1 <- data.frame(d=qbinom(ru1[,1],1,0.7),
                    d2=qbinom(ru1[,2],1,0.3),
                    h=qpois(ru1[,1],32),
                    x=qlnorm(ru1[,2],3,0.8))  #,3.5,0.9) #exp(qlnorm(u1[,1],1.3,0.65))
data1$y <- data1$h*data1$x
data1$y1 <- data1$d*data1$x
data1$y2 <- data1$d2*data1$d
data1$group <- 1

plot(density(data1$d2))


df <- rbind(data0,data1)

#u0 <- cbind(rank.wt(data0$x),rank.wt(data0$d),rank.wt(data0$d,left_limit_CDF = TRUE))
#u1 <- cbind(rank.wt(data1$x),rank.wt(data1$d),rank.wt(data1$d,left_limit_CDF = TRUE))
#plot(vinecop(u0,var_types=c("c","d"), family_set = "gaussian", weights=rep(1,n)))

#
# Decompose
set.seed(324)
VCopObj<- VC_deco(y1 ~ d:x,df, 
                  weights=NULL, group=group, 
                  var_types = c("d","d"), #c("c","d","c","d"),
                  nonpar_method = "constant",
                  structure=NA,
                  family_set="tll",
                  tree_crit = "rho",
                  mult=0.1,
                  multiply=1,
                  tau=c(1,seq(5,95,5),99)/100,
                  deco=TRUE, n=10000, quasi = TRUE, detailed_deco = FALSE,
                  bootstrap=FALSE, it=100, cores=4,
                  print_info=TRUE)

head(VCopObj$VCfit$u0)
VC_deco_qqplot(VCopObj,formula=NULL, type=1, logscale=TRUE, tau=3:99/100)
VC_deco_qqplot(VCopObj,formula=NULL, type=2, logscale=FALSE, tau=3:99/100)
VC_deco_qqplot(VCopObj,formula=NULL, type=3, logscale=TRUE, tau=3:99/100)
VC_deco_qqplot(VCopObj,formula=y1 ~ d:x, type=3, logscale=TRUE, tau=35:99/100) #+ ylim(-0.05, 0.05)


plotlist <- get_density_plots(VCopObj,varnames=c("x","d"))
compare_plots(plotlist,edge=1, groups=c(1976,2016))
plot(VCopObj$VCfit$VCfit0, edge_labels = "tau")
plot(VCopObj$VCfit$VCfit1, edge_labels = "tau")


### KD estimates: kdecop vs. bicop
#bi0 <- bicop_dist(family = "gumbel",rotation = 0,  parameters = 1.2,  var_types = c("c", "c"))
n=5000
bi0 <- bicop_dist(family = "gauss",rotation = 0,  parameters = 0.5,  var_types = c("c", "c"))
set.seed(123)
ru0 <- rbicop(n,bi0)
data2 <- data.frame(d1=qbinom(ru0[,1],1,0.7),
                    d2=qbinom(ru0[,2],1,0.3))

set.seed(4003)
data2$d3 <- add_jitter(data2$d1, theta=0, nu=5, bd=0.5, quasi=TRUE)
set.seed(4005)
data2$d4 <- add_jitter(data2$d2, theta=0, nu=5, bd=0.5, quasi=TRUE)

u1 <- apply(data2, 2, function(x) rank.wt(x, ties_method="random"))[,1:2]
u2 <- apply(data2, 2, function(x) rank.wt(x, ties_method="max"))[,3:4]
kd1 <- kdecop(u1, mult=0.25, method="TLL2", knots=30)
kd2 <- kdecop(u1, mult=0.5, method="TLL2", knots=30)
set.seed(123)
usim1 <- rkdecop(5000,kd1)
set.seed(123)
usim2 <- rkdecop(5000,kd2)
dsim1 <- data.frame(d1=qbinom(usim1[,1],1,0.7),
                    d2=qbinom(usim1[,2],1,0.3)) 
dsim2 <- data.frame(d1=qbinom(usim2[,1],1,0.7),
                    d2=qbinom(usim2[,2],1,0.3)) 
obs <- table(data2[,1:2])/5000
sim1 <- table(dsim1)/5000
sim2 <- table(dsim2)/5000
obs
sim1
sim2
plot(kd1)
plot(kd2)


n=10000
bi0 <- bicop_dist(family = "gauss",rotation = 0,  parameters = 0.5,  var_types = c("c", "c"))
set.seed(4003)
ru0 <- rbicop(n,bi0)
data2 <- data.frame(d1=qbinom(ru0[,1],1,0.7),
                    d2=qbinom(ru0[,2],1,0.3))
u3 <- cbind(apply(data2[,1:2], 2, function(x) rank.wt(x, ties_method="max")),
            apply(data2[,1:2], 2, function(x) rank.wt(x, ties_method="max", left_limit_CDF = TRUE)))
head(u3)
#u4 <- apply(apply(data2[,1:2],2, add_jitter), 2, function(x) rank.wt(x, ties_method="max"))
#head(u4)
kd3 <- bicop(u3, var_types=c("d","d"), mult=1, family_set="tll" , nonpar_method = "quadratic")
plot(kd3)
kd4 <- bicop(u3, var_types=c("d","d"), mult=0.5, family_set="tll", nonpar_method = "constant")
plot(kd4)

AIC(kd3)
AIC(kd4)

n2=10000
set.seed(4003)
usim3 <- rvinecopulib:::bicop_sim_cpp(kd4, n2, TRUE, rvinecopulib:::get_seeds()) #rbicop(n,kd3)
dsim3 <- data.frame(d1=qbinom(usim3[,1],1,0.7),
                    d2=qbinom(usim3[,2],1,0.3)) 
obs <- table(data2[,1:2])/n
sim <- table(dsim3)/n2
obs
sim

library(cctools)
dat <- cbind(ordered(data2$d1),
      ordered(data2$d2))
fit <- cckde(dat)
d <- dcckde(cbind(ordered(c(0,1,0,1)),ordered(c(0,0,1,1))), fit)
d/sum(d)

?locfit


set.seed(4003)
data2$d5 <- add_jitter(data2$d1, theta=0, nu=5, bd=0.5, quasi=TRUE)
set.seed(4005)
data2$d6 <- add_jitter(data2$d2, theta=0, nu=5, bd=0.5, quasi=TRUE)


set.seed(4003)
data2$d7 <- cont_conv(ordered(data2$d1), theta=0, nu=5, quasi=TRUE)
set.seed(4005)
data2$d8 <- cont_conv(ordered(data2$d2), theta=0, nu=5, quasi=TRUE)

head(data2)

# Pseudo-Copula Obs
#u_1 <- apply(data2, 2, function(x) rank.wt(x,ties_method="average"))
u_1 <- do.call("cbind",lapply(1:8, function(x) rank.wt(data2[,x], ties_method="average")))
u_2 <- do.call("cbind",lapply(1:8, function(x) pseudo_obs(data2[,x])))

u_1 <- cbind(do.call("cbind",lapply(1:2, function(x) rank.wt(data2[,x], ties_method="max"))),
             do.call("cbind",lapply(1:2, function(x) rank.wt(data2[,x], ties_method="max", left_limit_CDF = TRUE))),
             do.call("cbind",lapply(1:2, function(x) rank.wt(data2[,x], ties_method="random"))))

u_1 <- cbind(apply(data2, 2, function(x) rank.wt(x, ties_method="max")),
             apply(data2, 2, function(x) rank.wt(x, ties_method="max", left_limit_CDF = TRUE)))
             
u_3 <- NULL
for(i in 1:8){u_3 <- cbind(u_3,rank.wt(data2[,i], ties_method="average"))}
u_4 <- NULL
for(i in 1:8){u_4 <- cbind(u_4,pseudo_obs(data2[,i], ties_method="average"))}

plot(density(data2$d1))
plot(density(data2$d3))
plot(density(data2$d5))
plot(density(data2$d7))

plot(bi0)
plot(bicop(u_1[,1:4], family_set="gauss", var_types = c("d","d")))
plot(bicop(u_2[,1:2], family_set="tll"))
plot(bicop(u_2[,3:4], family_set="tll", mult=1, nonpar_method = "quadratic"))
plot(bicop(u_2[,5:6], family_set="tll", mult=1))
plot(bicop(u_2[,7:8], family_set="tll", mult=1))

bi4 <- bicop(u_2[,5:6], family_set="tll", mult=1)
plot(bi4)
bi5 <- bicop(u_1[,1:4], var_types=c("d","d"), family_set="tll", mult=1)
plot(bi5)
bi6 <- bicop(u_1[,5:6], var_types=c("c","c"), family_set="tll", mult=1)
plot(bi6)
set.seed(123)
u_sim <- rbicop(n,bi6)
d_sim <- data.frame(d1=qbinom(u_sim[,1],1,0.7),
                    d2=qbinom(u_sim[,2],1,0.3))
table(d_sim)/nrow(d_sim)-table(data2[,1:2])/nrow(data2[,1:2])

#table(data2)/rowSums(table(data2))
#table(data2)/matrix(rep(colSums(table(data2)),2), byrow=TRUE,2,2)

n=5000
bi0 <- bicop_dist(family = "gauss",rotation = 0,  parameters = 0.5,  var_types = c("c", "c"))
ru0 <- rbicop(n,bi0)
data2 <- data.frame(d1=qbinom(ru0[,1],1,0.7),
                    d2=qbinom(ru0[,2],1,0.4))
est0 <- kdecop(uA, method="TLL2", knots=30)
est1 <- kdecop(uA, method="TLL2nn", knots=30)
set.seed(123)
u_sim0 <- rkdecop(100000,est0)
set.seed(123)
u_sim1 <- rkdecop(100000,est1)
table(data.frame(d1=qbinom(u_sim0[,1],1,0.7),
                d2=qbinom(u_sim0[,2],1,0.4)))/100000 - table(data2)/n
table(data.frame(d1=qbinom(u_sim1[,1],1,0.7),
                 d2=qbinom(u_sim1[,2],1,0.4)))/100000 - table(data2)/n
plot(est1)

n=20000
it=100
bi0 <- bicop_dist(family = "gauss",rotation = 0,  parameters = 0.5,  var_types = c("c", "c"))
mse <- diffAa <- diffA <- diffB <- diffC <- matrix(rep(NA,it*4),it,4)
for(i in 1:it){
  print(i)
  ru0 <- rbicop(n,bi0)
  data2 <- data.frame(d1=qbinom(ru0[,1],1,0.7),
                      d2=qbinom(ru0[,2],1,0.4))
  uA <- cbind(rank.wt(add_jitter(data2$d1, theta=0, nu=5, bd=0.5, quasi=TRUE)),
              rank.wt(add_jitter(data2$d2, theta=0, nu=5, bd=0.5, quasi=TRUE)))
  #uB <- cbind(rank.wt(data2$d1, ties_method = "random"),
  #            rank.wt(data2$d2, ties_method = "random"))
  uB <- cbind(rank.wt(add_jitter(data2$d1, theta=0, nu=5, bd=0.5, quasi=FALSE)),
              rank.wt(add_jitter(data2$d2, theta=0, nu=5, bd=0.5, quasi=FALSE)))
  uC <-  cbind(rank.wt(data2$d1, ties_method = "average"),
               rank.wt(data2$d2, ties_method = "average"),
               rank.wt(data2$d1, ties_method = "average", left_limit_CDF = TRUE),
               rank.wt(data2$d2, ties_method = "average", left_limit_CDF = TRUE))
  bicA <- bicop(uA, var_types=c("c","c"), family_set="tll", mult=1)
  #bicAa <- bicop(uA, var_types=c("c","c"), family_set="tll", mult=3)
  bicAa <- bicop(uA, var_types=c("c","c"), family_set="tll", nonpar_method="constant", mult=1)
  bicB <- bicop(uB, var_types=c("c","c"), family_set="tll", mult=1)
  bicC <- bicop(uC, var_types=c("d","d"), family_set="tll", mult=1)
  u_simA <- rbicop(n,  bicA)
  u_simAa <- rbicop(n,  bicAa)
  u_simB <- rbicop(n,  bicB)
  u_simC <- rbicop(n,  bicC)
  d_simA <- data.frame(d1=qbinom(u_simA[,1],1,0.7),
                      d2=qbinom(u_simA[,2],1,0.4))
  d_simAa <- data.frame(d1=qbinom(u_simAa[,1],1,0.7),
                        d2=qbinom(u_simAa[,2],1,0.4))
  d_simB <- data.frame(d1=qbinom(u_simB[,1],1,0.7),
                        d2=qbinom(u_simB[,2],1,0.4))
  d_simC <- data.frame(d1=qbinom(u_simC[,1],1,0.7),
                       d2=qbinom(u_simC[,2],1,0.4))
  obs <- as.vector(table(data2[,1:2]))/n
  diffA[i,] <- as.vector(table(d_simA))/n-obs
  diffAa[i,]<- as.vector(table(d_simAa))/n-obs
  diffB[i,] <- as.vector(table(d_simB))/n-obs
  diffC[i,] <- as.vector(table(d_simC))/n-obs
  mse[i,1] <- mean((as.vector(table(d_simA))/n-obs)^2)
  mse[i,2] <- mean((as.vector(table(d_simAa))/n-obs)^2)
  mse[i,3] <- mean((as.vector(table(d_simB))/n-obs)^2)
  mse[i,4] <- mean((as.vector(table(d_simC))/n-obs)^2)
}
rbind(colMeans(diffA),
      colMeans(diffAa),
      colMeans(diffB),
      colMeans(diffC))
colMeans(mse)*100

library(cctools)
set.seed(4007)
n=100
bi0 <- bicop_dist(family = "gauss",rotation = 0,  parameters = 0.5,  var_types = c("c", "c"))
ru0 <- rbicop(n,bi0)
X <- qbinom(ru0[,1],1,0.3)
X1 <- rank.wt(X, ties_method="random")
X2 <- rank.wt(add_jitter(X, quasi=FALSE)) # rank.wt(cont_conv(ordered(X), quasi=FALSE))
Y <- qbinom(ru0[,2],1,0.3)
Y1 <- rank.wt(Y, ties_method="random")
Y2 <- rank.wt(add_jitter(Y, quasi=FALSE)) # rank.wt(cont_conv(ordered(Y), quasi=FALSE))
pairs_copula_data(cbind(X1,Y1,X2,Y2))
pairs_copula_data(cbind(X2,Y2))
plot(X2,Y2)

table(as.data.frame(cbind(X,Y)))/100

set.seed(123)
A <- cont_conv(ordered(Y), quasi=FALSE)
set.seed(123)
B <- add_jitter(Y, bd=0.5, quasi=FALSE)
cbind(A,B)

sum(apply(data2, 2, function(x) rank.wt(x, ties_method="average"))-do.call("cbind", lapply(1:2, function(x) rank.wt(data2[,x],ties_method="average"))))

test <- kdecop(u_2[,5:6])
plot(test)
test$bw

head(apply(data2[,1:2],2,function(x) x+1))
head(apply(data2[,7:8],2,function(x) rank.wt(x, ties_method="average")))
head(apply(data2[,7:8],2,function(x) pseudo_obs(x)))







plotdf <- VCopObj$VCdeco$quants_deco
plotdf <- melt(plotdf, id.vars=names(plotdf)[1:4],
               measure.vars=names(plotdf)[5:ncol(plotdf)],
               variable.name = "tau",
               value.name = "delta")
plotdf$tau <- as.numeric(as.character(plotdf$tau))

plot_agg <- ggplot(plotdf[which(plotdf$cat2=="aggregate"),], aes(tau,delta,col=names2,fill=names2)) + 
  geom_line() + geom_point() + 
  geom_hline(yintercept=0, col="lightgrey") + 
  facet_wrap(~names2, ncol=3) + theme(legend.position = "none")+ 
  labs(x="Quantile",y="Log difference", color="Effect", fill="Effect") 
plot_agg