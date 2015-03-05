library(statnet)
data(faux.mesa.high)
data(sampson)

Yfmh <- as.sociomatrix(faux.mesa.high)
diag(Yfmh) <- NA
Ysamp <- as.sociomatrix(samplike)
diag(Ysamp) <- NA

fmh <- faux.mesa.high

fmh_sdd <- sd(degree(fmh, gmode="graph"))
fmh_dens <- gden(fmh, mode="graph")
fmh_close <- centralization(fmh, closeness, mode="graph", cmode="undirected")

samp_sdd <- sd(degree(samplike, gmode="digraph"))
samp_sdideg <- sd(colSums(Ysamp, na.rm=T))
samp_sdodeg <- sd(rowSums(Ysamp, na.rm=T))
samp_dens <- gden(samplike, mode="digraph")
samp_close <- centralization(samplike, closeness, mode="digraph", cmode="directed")
samp_mut <- sum(Ysamp*t(Ysamp),na.rm=TRUE)/2
ns <- nrow(Ysamp)

## for histograms in workflow image
brg_dens <- NULL
brg_sdi <- NULL
brg_sdo <- NULL
for(s in 1:500){
  Ysim <- matrix(rbinom(ns^2,1,samp_dens),ns,ns)
  diag(Ysim)<- NA
  brg_dens <- rbind(brg_dens,mean(Ysim,na.rm=T))
  brg_sdo <- rbind(brg_sdo, sd(rowSums(Ysim, na.rm=TRUE)))
  brg_sdi <- rbind(brg_sdi, sd(colSums(Ysim, na.rm=TRUE)))
}


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
    
    Ybrg<-matrix(rbinom(n^2,1,obs_dens),n,n)
    diag(Ybrg)<-NA
    if(gmode=="graph"){Ybrg <- symmetrize(Ybrg,rule="upper")}
    brg_sdd <- rbind(brg_sdd, sd(degree(Ybrg, gmode=gmode)))
    brg_dens <- rbind(brg_dens, gden(Ybrg, mode=gmode))
    brg_close <- rbind(brg_close, centralization(Ybrg, closeness,mode=gmode,cmode=cmode))
    brg_sdo <- rbind(brg_sdo, sd(rowSums(Ybrg, na.rm=TRUE)))
    brg_sdi <- rbind(brg_sdi, sd(colSums(Ybrg, na.rm=TRUE)))
    brg_mut <- rbind(brg_mut, sum(Ybrg*t(Ybrg),na.rm=T)/2)
  }
  
  dat <- list()
  dat$cug <- cbind(cug_sdd, cug_close, cug_sdo, cug_sdi, cug_mut)
  dat$brg <- cbind(brg_sdd, brg_close, brg_sdo, brg_sdi, brg_mut, brg_dens)
  return(dat)
}


## histograms in null distribution section
fmhsim <- tsim(Yfmh, gmode="graph", cmode="undirected")

BRGcol <- adjustcolor("darkred", alpha.f=0.7)
CUGcol <- adjustcolor("sienna1", alpha.f=0.9)
obsblue <- "#4D94FF"

par(mfrow=c(1,1), mar=mar, xpd=FALSE)
layout(rbind(c(1,2,3)), widths=c(2,2,1.1))
hist(fmhsim$brg[,6], col=BRGcol, ylab=NULL, main=NULL, xlab="density")
abline(v=fmh_dens, col=obsblue, lwd=2)

hist(fmhsim$cug[,1], col=CUGcol, ylab=NULL, main=NULL, xlab="sd(degree)", xlim=range(fmhsim$cug[,1],fmhsim$brg[,1],fmh_sdd+.2), ylim=c(0,150))
hist(fmhsim$brg[,1], col=BRGcol, ylab=NULL, main=NULL, xlab="sd(degree)", add=T)
abline(v=fmh_sdd, col=obsblue, lwd=2)
plot.new(); par(xpd=TRUE)
legend(x="top",legend=c("observed","BRG(rho)","CUG(s)"), fill=c(0,"darkred","sienna1"), border=c(0,"black","black"),
       col=c(obsblue,0,0), lty=c(1,0,0),lwd=c(2,0,0), merge=TRUE, bty="n", cex=1.5)


par(mfrow=c(1,1), xpd=FALSE)
#hist(brg_close, col=BRGcol, ylab=NULL, main=NULL, xlab="closeness centrality")
#all zeros because graph is disconnected

### Blockmodels
sex <- get.vertex.attribute(fmh,"Sex")
p00 <- mean(Yfmh[sex=="F",sex=="F"], na.rm=T)
p01 <- mean(Yfmh[sex=="F",sex=="M"], na.rm=T)
p10 <- mean(Yfmh[sex=="M",sex=="F"], na.rm=T)
p11 <- mean(Yfmh[sex=="M",sex=="M"], na.rm=T)



### Regression
n <- nrow(Yfmh)
sex01 <- 1*(sex=="M")
sexR <- matrix(sex01,nrow=n,ncol=n)
sexC <- t(sexR)
fit <- glm(c(Yfmh) ~ c(sexR) + c(sexC) + c(sexR)*c(sexC), family=binomial)
summary(fit)
exp(fit$coef)
P <- fit$fitted.values

sim_sdd <- NULL
sim_sdi <- NULL

for(j in 1:500){
  Ysim <- matrix(0,n,n)
  diag(Ysim)<- NA
  Ysim[!is.na(Ysim)] <- rbinom(n*(n-1),1,P)
  Ysim <- symmetrize(Ysim,rule="upper")
  sim_sdd <- rbind(sim_sdd, sd(rowSums(Ysim,na.rm=T)))
}




### Row Column Effects
sampsimdat <- tsim(Ysamp, gmode="digraph", cmode="directed")

samp_sdo <- sd(rowSums(Ysamp, na.rm=TRUE))
samp_sdi <- sd(colSums(Ysamp, na.rm=TRUE))
samp_mut <- sum(Ysamp*t(Ysamp),na.rm=TRUE)/2
par(mfrow=c(1,3))
hist(sampsimdat$cug[,3], col=CUGcol, xlab="sd(out-degree)", ylab=NULL,main=NULL,
     xlim=c(0,max(sampsimdat$cug[,3],samp_sdo)))
hist(sampsimdat$brg[,3], col=BRGcol, add=T)
abline(v=samp_sdo,col=obsblue, lwd=2)

hist(sampsimdat$cug[,4], col=CUGcol, xlab="sd(in-degree)", ylab=NULL,main=NULL,
     xlim=c(0,max(sampsimdat$cug[,4],samp_sdi)))
hist(sampsimdat$brg[,4], col=BRGcol, add=T)
abline(v=samp_sdi,col=obsblue, lwd=2)

hist(sampsimdat$cug[,5], col=CUGcol, xlab="mutual dyads", ylab=NULL, main=NULL,
     breaks=c(0:24),xlim=c(0,samp_mut))
hist(sampsimdat$brg[,5], col=BRGcol, breaks=c(0:24), add=T)
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
hist(sim_sdo, col="gray", xlab="sd(out-degree)",ylab=NULL, main=NULL,
     xlim=range(sim_sdo,samp_sdo))
abline(v=samp_sdo,col=obsblue,lwd=2)
hist(sim_sdi, col="gray", xlab="sd(in-degree)",ylab=NULL, main=NULL,
     xlim=range(sim_sdi,samp_sdi))
abline(v=samp_sdi,col=obsblue,lwd=2)
hist(sim_mut, col="gray", xlab="mutual dyads",ylab=NULL, main=NULL,
     xlim=range(sim_mut,samp_mut))
abline(v=samp_mut,col=obsblue,lwd=2)
