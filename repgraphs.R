library(statnet)
data(faux.mesa.high)
data(sampson)

Yfmh <- as.sociomatrix(faux.mesa.high)
diag(Yfmh) <- NA
Ysamp <- as.sociomatrix(samplike)
diag(Ysamp) <- NA

fmh <- faux.mesa.high
n <- nrow(Yfmh)
fmh_sdd <- sd(degree(fmh, gmode="graph"))
fmh_dens <- gden(fmh, mode="graph")
fmh_close <- centralization(fmh, closeness, mode="graph", cmode="undirected")
fmh_ingroup <- sum(diag(mixingmatrix(fmh, "Sex")[[2]]))
fmh_mm <- mixingmatrix(fmh,"Sex")[[2]][[2,2]]
fmh_ff <- mixingmatrix(fmh,"Sex")[[2]][[1,1]]
fmh_tri <- summary(fmh ~ triangle)

samp_sdd <- sd(degree(samplike, gmode="digraph"))
samp_sdi <- sd(colSums(Ysamp, na.rm=T))
samp_sdo <- sd(rowSums(Ysamp, na.rm=T))
samp_dens <- gden(samplike, mode="digraph")
samp_close <- centralization(samplike, closeness, mode="digraph", cmode="directed")
samp_mut <- sum(Ysamp*t(Ysamp),na.rm=TRUE)/2
ns <- nrow(Ysamp)

## for histograms in workflow image
brg_dens <- NULL
brg_sdd <- NULL
for(s in 1:500){
  Ysim <- matrix(rbinom(n^2,1,fmh_dens),n,n)
  diag(Ysim)<- NA
  brg_dens <- rbind(brg_dens,mean(Ysim,na.rm=T))
  brg_sdd <- rbind(brg_sdd, sd(rowSums(Ysim, na.rm=TRUE)))
}
BRGcol <- adjustcolor("firebrick")
CUGcol <- adjustcolor("orangered")
obsblue <- "royalblue2"
tgray <- adjustcolor("gray",alpha.f=0.4)
par(mfrow=c(1,2))
hist(brg_dens,col=tgray,border=BRGcol,xlab="density",ylab=NULL,main=NULL)
abline(v=fmh_dens,col=obsblue,lwd=2)
hist(brg_sdd,col=tgray,border=BRGcol,xlab="sd(degree)",ylab=NULL,main=NULL,
     xlim=range(brg_sdd,fmh_sdd))
abline(v=fmh_sdd,col=obsblue,lwd=2)

tsim <- function(Y, gmode, cmode){
  n <- nrow(Y)
  if(gmode=="graph"){
    s <- sum(Y, na.rm=T)/2
  } else {
    s <- sum(Y, na.rm=T)
  }
  obs_dens <- mean(Y, na.rm=TRUE)
  
  brg_dens <- NULL
  brg_close <- NULL
  brg_sdd <- NULL
  cug_sdd <- NULL
  cug_close <- NULL
  cug_sdi <- NULL
  cug_sdo <- NULL
  brg_sdi <- NULL
  brg_sdo <- NULL
  brg_mut <- NULL
  cug_mut <- NULL
  brg_ingroupsex <- NULL
  cug_ingroupsex <- NULL
  
  for(i in 1:500){
    Ycug <- matrix(0,n,n)
    diag(Ycug)<- NA
    if(gmode=="graph"){
      Ycug[upper.tri(Ycug)] <- sample(c(rep(1,s),rep(0,n*(n-1)/2-s)))
      Ycug <- symmetrize(Ycug, rule="upper")
    } else {
      Ycug[!is.na(Ycug)] <- sample(c(rep(1,s), rep(0,n*(n-1)-s)))
    }
    cug_sdd <- rbind(cug_sdd, sd(degree(Ycug, gmode=gmode)))
    cug_close <- rbind(cug_close, centralization(Ycug, closeness,mode=gmode,cmode=cmode))
    cug_sdo <- rbind(cug_sdo, sd(rowSums(Ycug, na.rm=TRUE)))
    cug_sdi <- rbind(cug_sdi, sd(colSums(Ycug, na.rm=TRUE)))
    cug_mut <- rbind(cug_mut, sum(Ycug*t(Ycug),na.rm=T)/2)
    cug_ingroupsex <- rbind(cug_ingroupsex, sum(Ycug[sex01==1,sex01==1],na.rm=T)+sum(Ycug[sex01==0,sex01==0],na.rm=T))
    
    Ybrg<-matrix(rbinom(n^2,1,obs_dens),n,n)
    diag(Ybrg)<-NA
    if(gmode=="graph"){Ybrg <- symmetrize(Ybrg,rule="upper")}
    brg_sdd <- rbind(brg_sdd, sd(degree(Ybrg, gmode=gmode)))
    brg_dens <- rbind(brg_dens, gden(Ybrg, mode=gmode))
    brg_close <- rbind(brg_close, centralization(Ybrg, closeness,mode=gmode,cmode=cmode))
    brg_sdo <- rbind(brg_sdo, sd(rowSums(Ybrg, na.rm=TRUE)))
    brg_sdi <- rbind(brg_sdi, sd(colSums(Ybrg, na.rm=TRUE)))
    brg_mut <- rbind(brg_mut, sum(Ybrg*t(Ybrg),na.rm=T)/2)
    brg_ingroupsex <- rbind(brg_ingroupsex, sum(Ybrg[sex01==1,sex01==1],na.rm=T)+sum(Ybrg[sex01==0,sex01==0],na.rm=T))
  }
  
  dat <- list()
  dat$cug <- cbind(cug_sdd, cug_close, cug_sdo, cug_sdi, cug_mut, cug_ingroupsex)
  dat$brg <- cbind(brg_sdd, brg_close, brg_sdo, brg_sdi, brg_mut, brg_ingroupsex, brg_dens)
  return(dat)
}


## histograms in null distribution section
fmhsim <- tsim(Yfmh, gmode="graph", cmode="undirected")

BRGcol <- adjustcolor("firebrick")
CUGcol <- adjustcolor("orangered")
obsblue <- "royalblue2"
tgray <- adjustcolor("gray",alpha.f=0.4)

par(mfrow=c(1,1), xpd=FALSE)
layout(rbind(c(1,2,3,4)), widths=c(2,2,2,1.1))
hist(fmhsim$brg[,7], col=tgray,border=BRGcol,ylab=NULL,main=NULL,xlab="density")
abline(v=fmh_dens, col=obsblue, lwd=2)

hist(fmhsim$cug[,1], col=tgray,border=CUGcol, ylab=NULL, main=NULL, xlab="sd(degree)", xlim=range(fmhsim$cug[,1],fmhsim$brg[,1],fmh_sdd+.2), ylim=c(0,150))
hist(fmhsim$brg[,1], col=tgray,border=BRGcol, ylab=NULL, main=NULL, xlab="sd(degree)", add=T)
abline(v=fmh_sdd, col=obsblue, lwd=2)

hist(fmhsim$cug[,6], col=tgray,border=CUGcol, ylab=NULL, main=NULL, xlab="within group ties (sex)", xlim=range(fmhsim$cug[,6],fmhsim$brg[,6],fmh_ingroup))
hist(fmhsim$brg[,6], col=tgray,border=BRGcol, ylab=NULL, main=NULL, xlab="within group ties (sex)", add=T)
abline(v=fmh_ingroup, col=obsblue, lwd=2)

plot.new(); par(xpd=TRUE)
legend(x="top",legend=c("observed","BRG(rho)","CUG(s)"), fill=c(0,tgray,tgray), border=c(0,BRGcol,CUGcol),
       col=c(obsblue,0,0), lty=c(1,0,0),lwd=c(2,0,0), merge=TRUE, bty="n", cex=1.5)


par(mfrow=c(1,1), xpd=FALSE)
#hist(brg_close, col=BRGcol, ylab=NULL, main=NULL, xlab="closeness centrality")
#all zeros because graph is disconnected


# Blockmodels -------------------------------------------------------------


sex <- get.vertex.attribute(fmh,"Sex")
p00 <- mean(Yfmh[sex=="F",sex=="F"], na.rm=T)
p01 <- mean(Yfmh[sex=="F",sex=="M"], na.rm=T)
p10 <- mean(Yfmh[sex=="M",sex=="F"], na.rm=T)
p11 <- mean(Yfmh[sex=="M",sex=="M"], na.rm=T)




# Regression --------------------------------------------------------------


b0 <- 0.0009
par(mfrow=c(2,1))
b1 <- 2
x <- seq(-2.5,3,0.05)
mu <- exp(b0 + b1*x)/(1+exp(b0 + b1*x))
plot(x,mu,type="l",xaxt="n", main="beta1 > 0", ylab="Pr(Y=1 | x)")

n <- nrow(Yfmh)
sex <- get.vertex.attribute(fmh,"Sex")
sex01 <- 1*(sex=="M")
sexR <- matrix(sex01,nrow=n,ncol=n)
sexC <- t(sexR)
fit <- glm(c(Yfmh) ~ c(sexR) + c(sexC) + c(sexR*sexC), family=binomial)
summary(fit) #same
fit2 <- ergm(fmh ~ edges + nodefactor("Sex") + nodematch("Sex",keep=2))
summary(fit2) #same
fit3 <- glm(c(Yfmh) ~ c(sexR) + c(sexR*sexC), family=binomial)
summary(fit3) #different
fit4 <- glm(c(Yfmh) ~ c(sexR) + c(sexR)*c(sexC), family=binomial)
summary(fit4) #same

exp(fit$coef)
P <- fit$fitted.values

sim_sdd <- NULL
sim_ingroupsex <- NULL
sim_mm <- NULL
sim_ff <- NULL
sim_tri <- NULL

#from glm fitted values
for(j in 1:1000){
  Ysim <- matrix(0,n,n)
  diag(Ysim)<- NA
  Ysim[!is.na(Ysim)] <- rbinom(n*(n-1),1,P)
  Ysim <- symmetrize(Ysim,rule="upper")
  sim_sdd <- rbind(sim_sdd, sd(rowSums(Ysim,na.rm=T)))
  sim_ingroupsex <- rbind(sim_ingroupsex, sum(Ysim[sex01==1,sex01==1],na.rm=T)+sum(Ysim[sex01==0,sex01==0],na.rm=T))
  sim_mm <- rbind(sim_mm, sum(Ysim[sex01==1,sex01==1],na.rm=T))
  sim_ff <- rbind(sim_ff, sum(Ysim[sex01==0,sex01==0],na.rm=T))
  diag(Ysim) <- 0
  sim_tri <- rbind(sim_tri,sum(diag(Ysim%*%Ysim%*%Ysim)))
}
# from ergm
for(j in 1:500){
  Ysim <- simulate(fit2)
  sim_sdd <- rbind(sim_sdd, sd(summary(Ysim ~ degree(0:15))))
  sim_ingroupsex <- rbind(sim_ingroupsex, summary(Ysim ~ nodematch("Sex")))
  sim_mm <- rbind(sim_mm, summary(Ysim ~ nodematch("Sex",diff=TRUE,keep=1)))
  sim_ff <- rbind(sim_ff, summary(Ysim ~ nodematch("Sex",diff=TRUE,keep=2)))
  sim_tri <- rbind(sim_tri,summary(Ysim ~ triangles))
}

par(mfrow=c(1,5))
hist(sim_sdd,col=tgray,xlim=range(sim_sdd,fmh_sdd), border=BRGcol,xlab="sd(degree)",
     ylab=NULL,main=NULL)
abline(v=fmh_sdd,col=obsblue, lwd=2)
hist(sim_ingroupsex,col=tgray,xlim=range(sim_ingroupsex,fmh_ingroup), border=BRGcol,
     xlab="within group ties (sex)",ylab=NULL,main=NULL)
abline(v=fmh_ingroup,col=obsblue, lwd=2)
hist(sim_mm,col=tgray,xlim=range(sim_mm,fmh_mm), border=BRGcol,
     xlab="male-male ties",ylab=NULL,main=NULL)
abline(v=fmh_mm,col=obsblue, lwd=2)
hist(sim_ff,col=tgray,xlim=range(sim_ff,fmh_ff), border=BRGcol,
     xlab="female-female ties",ylab=NULL,main=NULL)
abline(v=fmh_ff,col=obsblue, lwd=2)
hist(sim_tri, col=tgray,xlim=range(sim_tri,fmh_tri), border=BRGcol,xlab="triangles",
     ylab=NULL,main=NULL)
abline(v=fmh_tri, col=obsblue, lwd=2)

BRGcol <- adjustcolor("firebrick")
CUGcol <- adjustcolor("orangered")
obsblue <- "royalblue2"
tgray <- adjustcolor("gray",alpha.f=0.4)



# Row Column Effects ------------------------------------------------------


sampsimdat <- tsim(Ysamp, gmode="digraph", cmode="directed")

samp_sdo <- sd(rowSums(Ysamp, na.rm=TRUE))
samp_sdi <- sd(colSums(Ysamp, na.rm=TRUE))
samp_mut <- sum(Ysamp*t(Ysamp),na.rm=TRUE)/2
par(mfrow=c(1,3))
hist(sampsimdat$cug[,3], border=CUGcol,col=tgray,xlab="sd(out-degree)", 
     ylab=NULL,main=NULL,
     xlim=c(0,max(sampsimdat$cug[,3],samp_sdo)))
hist(sampsimdat$brg[,3], border=BRGcol,col=tgray, add=T)
abline(v=samp_sdo,col=obsblue, lwd=2)

hist(sampsimdat$cug[,4], border=CUGcol,col=tgray,xlab="sd(in-degree)",
     ylab=NULL,main=NULL,
     xlim=c(0,max(sampsimdat$cug[,4],samp_sdi)))
hist(sampsimdat$brg[,4], border=BRGcol,col=tgray, add=T)
abline(v=samp_sdi,col=obsblue, lwd=2)

hist(sampsimdat$cug[,5], border=CUGcol,col=tgray, xlab="mutual dyads",
     ylab=NULL, main=NULL,
     breaks=c(0:24),xlim=c(0,samp_mut))
hist(sampsimdat$brg[,5], border=BRGcol,col=tgray, breaks=c(0:24), add=T)
abline(v=samp_mut, col=obsblue, lwd=2)
par(mfrow=c(1,1)); title(main="sampson data"); 

Ridx <- matrix(1:ns,ns,ns)
Cidx <- t(Ridx)
rce.fit <- glm(c(Ysamp) ~ factor(c(Ridx))+factor(c(Cidx)), family=binomial)
summary(rce.fit)
P.rce <- rce.fit$fitted.values
sim_sdi <- NULL
sim_sdo <- NULL
sim_mut <- NULL
for(j in 1:1000){
  Ysim <- matrix(0,ns,ns)
  diag(Ysim)<- NA
  Ysim[!is.na(Ysim)] <- rbinom(ns*(ns-1),1,P.rce)
  sim_sdo <- rbind(sim_sdo, sd(rowSums(Ysim,na.rm=T)))
  sim_sdi <- rbind(sim_sdi, sd(colSums(Ysim,na.rm=T)))
  sim_mut <- rbind(sim_mut, sum(Ysim*t(Ysim),na.rm=T)/2)
}
par(mfrow=c(1,3))
hist(sim_sdo, col=tgray,border="darkred", xlab="sd(out-degree)",ylab=NULL, main=NULL,
     xlim=range(sim_sdo,samp_sdo))
abline(v=samp_sdo,col=obsblue,lwd=2)

hist(sim_sdi, col=tgray,border="darkred", xlab="sd(in-degree)",ylab=NULL, main=NULL,
     xlim=range(sim_sdi,samp_sdi))
abline(v=samp_sdi,col=obsblue,lwd=2)

hist(sim_mut, col=tgray,border="darkred", xlab="mutual dyads",ylab=NULL, main=NULL,
     xlim=range(sim_mut,samp_mut))
abline(v=samp_mut,col=obsblue,lwd=2)

### p1

fitp1 <- ergm(samplike ~ edges + sender + receiver + mutual)
sim_sdi <- NULL
sim_sdo <- NULL
sim_mut <- NULL
for(k in 1:1000){
  Ysim <- simulate(fitp1)
  sim_sdi <- rbind(sim_sdi,sd(summary(Ysim ~ idegree(0:15))))
  sim_sdo <- rbind(sim_sdo,sd(summary(Ysim ~ odegree(0:15))))
  sim_mut <- rbind(sim_mut,summary(Ysim ~ mutual))
}

par(mfrow=c(1,3))
hist(sim_sdo, col=tgray,border="darkred", xlab="sd(out-degree)",ylab=NULL, main=NULL,
     xlim=range(sim_sdo,samp_sdo))
abline(v=samp_sdo,col=obsblue,lwd=2)

hist(sim_sdi, col=tgray,border="darkred", xlab="sd(in-degree)",ylab=NULL, main=NULL,
     xlim=range(sim_sdi,samp_sdi))
abline(v=samp_sdi,col=obsblue,lwd=2)

hist(sim_mut, col=tgray,border="darkred", xlab="mutual dyads",ylab=NULL, main=NULL,
     xlim=range(sim_mut,samp_mut))
abline(v=samp_mut,col=obsblue,lwd=2)
par(mfrow=c(2,2))
plot(gof(fitp1))
