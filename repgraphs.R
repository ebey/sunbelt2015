library(statnet)
data(faux.mesa.high)
fmh <- faux.mesa.high
n <- fmh$gal$n

obs_sdd <- sd(degree(fmh, gmode="graph"))
obs_theta <- gden(fmh, mode="graph")
obs_close <- centralization(fmh, closeness, mode="graph", cmode="undirected")

brg_theta <- NULL
brg_close <- NULL
brg_sdd <- NULL
cug_sdd <- NULL
cug_close <- NULL

for(i in 1:500){
  Ycug <- matrix(0,n,n)
  diag(Ycug)<- NA
  Ycug[upper.tri(Ycug)] <- sample(c(rep(1,203),rep(0,n*(n-1)/2-203)))
  Ycug <- symmetrize(Ycug, rule="upper")
  cug_sdd <- rbind(cug_sdd, sd(degree(Ycug, gmode="graph")))
  cug_close <- rbind(cug_close, centralization(Ycug, closeness,mode="graph",cmode="undirected"))
  
  Ybrg<-matrix(rbinom(n^2,1,obs_theta),n,n)
  diag(Ybrg)<-NA
  Ybrg <- symmetrize(Ybrg,rule="upper")
  brg_sdd <- rbind(brg_sdd, sd(degree(Ybrg, gmode="graph")))
  brg_theta <- rbind(brg_theta, gden(Ybrg, mode="graph"))
  brg_close <- rbind(brg_close, centralization(Ybrg, closeness,mode="graph",cmode="undirected"))
}

CUGcol <- adjustcolor("firebrick", alpha.f=0.55)
BRGcol <- adjustcolor("orangered", alpha.f=0.55)

hist(brg_sdd, col=BRGcol, ylab=NULL, main=NULL, xlab="sd(degree)", xlim=range(brg_sdd,cug_sdd,obs_sdd+.2), ylim=c(0,150))
hist(cug_sdd, col=CUGcol, ylab=NULL, main=NULL, xlab="sd(degree)", add=T)
abline(v=obs_sdd, col="blue", lty=2, lwd=2)

hist(brg_theta, col=BRGcol, ylab=NULL, main=NULL, xlab="density")
abline(v=obs_theta, col="blue", lty=2, lwd=2)

hist(brg_close, col=BRGcol, ylab=NULL, main=NULL, xlab="closeness centrality")
#all zeros because graph is disconnected

### Blockmodels
Yobs <- as.sociomatrix(fmh)
sex <- get.vertex.attribute(fmh,"Sex")
p00 <- mean(Yobs[sex=="F",sex=="F"])
p01 <- mean(Yobs[sex=="F",sex=="M"])
p10 <- mean(Yobs[sex=="M",sex=="F"])
p11 <- mean(Yobs[sex=="M",sex=="M"])
