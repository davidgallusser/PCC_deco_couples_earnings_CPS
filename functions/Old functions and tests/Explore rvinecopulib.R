# Explore vinecop / vinecopulib
# https://vinecopulib.github.io/


library("rvinecopulib")
library("kdecopula")
data("wdbc")

library(sloop)
s3_methods_class("vinecop")
s3_get_method(plot.vinecop)

# Non-parametric vinc copula modelling 
u <- pobs(wdbc[,c(4:7)])
w <- runif(nrow(u),0.5,1.5)
fit <- vinecop(u, family_set = "tll", structure=NA,
               nonpar_method="quadratic",
               weights=w, tree_crit="rho", show_trace=TRUE)
summary(fit)
plot(fit, tree=c(1,2,3))
contour(fit)

structure <- dvine_structure(4:1)
structure <- cvine_structure(4:1)


fit <- vinecop(u, family_set = "tll", structure=structure,
               nonpar_method="quadratic",
               weights=w, tree_crit="rho", show_trace=TRUE)
plot(fit, tree=c(1:3),edge_labels = "tau")


M <- as_rvine_matrix(cvine_structure(1:4, trunc_lvl = 1))
rvinecopulib:::get_graph()


require("rvinecopulib")

setwd(paste0(sav_d,datafiles_d))
load(paste0("CPScouples",2016,".rda"))
jit <- c(2,4)
data1 <- cdata[sample(1:nrow(cdata),1000),c("EARNINGS","WAGE_m","ANNUALHOURS_m", "WAGE_f","ANNUALHOURS_f","ASECWTH_m")]
w1 <- data1[,"ASECWTH_m"]
names(data1)[6] <- ("(weights)")
data1[,jit+1] <- apply(data1[,jit+1], 2, function(x) add_jitter(x))
u1  <- apply(data1[,2:5], 2, function(x) rank.wt(x,w1))

load(paste0("CPScouples",1976,".rda"))
jit <- c(2,4)
data0 <- cdata[sample(1:nrow(cdata),1000),c("EARNINGS","WAGE_m","ANNUALHOURS_m", "WAGE_f","ANNUALHOURS_f","ASECWTH_m")]
w0 <- data0[,"ASECWTH_m"]
names(data0)[6] <- ("(weights)")
data0[,jit+1] <- apply(data0[,jit+1], 2, function(x) add_jitter(x))
u0  <- apply(data0[,2:5], 2, function(x) rank.wt(x,w0))


VCfit0 <- vinecop(u0, family_set = "tll", nonpar_method="quadratic", weights=w0, tree_crit="rho", show_trace=TRUE)
summary(VCfit0)

VCfit1 <- vinecop(u1,family_set = "tll", nonpar_method="quadratic", weights=w1, tree_crit="rho", show_trace=TRUE)
summary(VCfit1)


rvinecop(10, VCfit1, qrng = FALSE, cores = 1)



plot(bicop(u0[,c(1,3)], family_set = "tll", nonpar_method="linear", weights=w0))
plot(bicop(u1[,c(2,1)], family_set = "tll", nonpar_method="quadratic", weights=w1))
plot(bicop(u1[,c(1,2)], family_set = "tll", nonpar_method="quadratic", weights=w1))

str <- get_structure(VCfit1)
as_rvine_matrix(str)

m <- get_matrix(VCfit1)
m
as_rvine_structure(m)

c <- get_pair_copula(VCfit1, tree=1, edge=1)
plot(c)
su <- summary(VCfit1)
su

fit <- VCfit1



str$struct_array
predict(VCfit1)
test <- as_rvine_matrix(str)

c1 <- bicop(u1[,c(2,1)], family_set = "tll", nonpar_method="quadratic")
c2 <- bicop(u1[,c(1,2)], family_set = "tll", nonpar_method="quadratic")
plot(c1)
plot(c2)
c2$parameters <- t(c1$parameters)
plot(c2)



setwd(paste0(sav_d,datafiles_d))
load(paste0("CPScouples",2016,".rda"))
jit <- c(2,4)
data1 <- cdata[sample(1:nrow(cdata),1000),c("EARNINGS","WAGE_m","ANNUALHOURS_m", "WAGE_f","ANNUALHOURS_f","ASECWTH_m")]
w1 <- data1[,"ASECWTH_m"]
names(data1)[6] <- ("(weights)")
data1[,jit+1] <- apply(data1[,jit+1], 2, function(x) add_jitter(x))
u1  <- apply(data1[,2:5], 2, function(x) rank.wt(x,w))

VCfit1 <- vinecop(u1, weights=w1,  treecrit="rho", info=TRUE)


load(paste0("CPScouples",1976,".rda"))
jit <- c(2,4)
data0 <- cdata[sample(1:nrow(cdata),1000),c("EARNINGS","WAGE_m","ANNUALHOURS_m", "WAGE_f","ANNUALHOURS_f","ASECWTH_m")]
w0 <- data0[,"ASECWTH_m"]
names(data0)[6] <- ("(weights)")
data0[,jit+1] <- apply(data0[,jit+1], 2, function(x) add_jitter(x))
u0  <- apply(data0[,2:5], 2, function(x) rank.wt(x,w))

VCfit0 <- wtd.kdevinecop(u0, weights=w0,   treecrit="rho", info=TRUE)


VineCopula:::findMaxTree()
