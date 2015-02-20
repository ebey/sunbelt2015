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
samp_dens <- gden(samplike, mode="digraph")
samp_close <- centralization(samplike, closeness, mode="digraph", cmode="directed")

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
    
    Ybrg<-matrix(rbinom(n^2,1,obs_dens),n,n)
    diag(Ybrg)<-NA
    Ybrg <- symmetrize(Ybrg,rule="upper")
    brg_sdd <- rbind(brg_sdd, sd(degree(Ybrg, gmode=gmode)))
    brg_dens <- rbind(brg_dens, gden(Ybrg, mode=gmode))
    brg_close <- rbind(brg_close, centralization(Ybrg, closeness,mode=gmode,cmode=cmode))
    brg_sdo <- rbind(brg_sdo, sd(rowSums(Ybrg, na.rm=TRUE)))
    brg_sdi <- rbind(brg_sdi, sd(colSums(Ybrg, na.rm=TRUE)))
  }
  
  dat <- list()
  dat$cug <- cbind(cug_sdd, cug_close, cug_sdo, cug_sdi)
  dat$brg <- cbind(brg_sdd, brg_close, brg_sdo, brg_sdi, brg_dens)
  return(dat)
}



fmhsim <- tsim(Yfmh, gmode="graph", cmode="undirected")


CUGcol <- adjustcolor("darkred", alpha.f=0.7)
BRGcol <- adjustcolor("sienna1", alpha.f=0.9)

par(mfrow=c(1,1), mar=mar, xpd=FALSE)
layout(rbind(c(1,2,3)), widths=c(2,2,1.1))
hist(brg_theta, col=BRGcol, ylab=NULL, main=NULL, xlab="density")
abline(v=obs_theta, col="blue", lty=2, lwd=2)

hist(brg_sdd, col=BRGcol, ylab=NULL, main=NULL, xlab="sd(degree)", xlim=range(brg_sdd,cug_sdd,obs_sdd+.2), ylim=c(0,150))
hist(cug_sdd, col=CUGcol, ylab=NULL, main=NULL, xlab="sd(degree)", add=T)
abline(v=obs_sdd, col="blue", lty=2, lwd=2)
plot.new(); par(xpd=TRUE)
legend(x="top",legend=c("observed","BRG","CUG(n,s)"), fill=c(0,"sienna1","darkred"), border=c(0,"black","black"),
       col=c("blue",0,0), lty=c(2,0,0),lwd=c(2,0,0), merge=TRUE, bty="n")


par(mfrow=c(1,1), xpd=FALSE)
#hist(brg_close, col=BRGcol, ylab=NULL, main=NULL, xlab="closeness centrality")
#all zeros because graph is disconnected

### Blockmodels
sex <- get.vertex.attribute(fmh,"Sex")
p00 <- mean(Yfmh[sex=="F",sex=="F"])
p01 <- mean(Yfmh[sex=="F",sex=="M"])
p10 <- mean(Yfmh[sex=="M",sex=="F"])
p11 <- mean(Yfmh[sex=="M",sex=="M"])



### Regression
sex01 <- 1*(sex=="M")
sexR <- matrix(sex01,nrow=n,ncol=n)
sexC <- t(sexR)
fit <- glm(c(Yobs) ~ c(sexR) + c(sexR)*c(sexC), family=binomial)
summary(fit)
exp(fit$coef)


sampsimdat <- tsim(Ysamp, gmode="digraph", cmode="directed")

### Row Column Effects
samp_sdo <- sd(rowSums(Ysamp, na.rm=TRUE))
samp_sdi <- sd(colSums(Ysamp, na.rm=TRUE))
par(mfrow=c(1,3))
hist(sampsimdat$brg[,3], col=BRGcol, xlab="sd(out-degree)", ylab=NULL,main=NULL,
     xlim=range(c(sampsimdat$cug[,3],samp_sdo)))
hist(sampsimdat$cug[,3], col=CUGcol, add=T)
abline(v=samp_sdo,col="blue", lty=2, lwd=2)

hist(sampsimdat$brg[,4], col=BRGcol, xlab="sd(in-degree)", ylab=NULL,main=NULL,
     xlim=range(c(sampsimdat$cug[,4],samp_sdi)))
hist(sampsimdat$cug[,4], col=CUGcol, add=T)
abline(v=samp_sdi,col="blue", lty=2, lwd=2)
plot.new()
par(mar=c(5.1,0,4.1,2.1))
legend(x="topleft",legend=c("observed","BRG","CUG(n,s)"), fill=c(0,"sienna1","darkred"), border=c(0,"black","black"),
       col=c("blue",0,0), lty=c(2,0,0),lwd=c(2,0,0), merge=TRUE, bty="n")
par(mfrow=c(1,1)); title(main="sampson data"); 
