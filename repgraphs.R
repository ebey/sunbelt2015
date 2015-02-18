library(statnet)
data(faux.mesa.high)
n <- faux.mesa.high$gal$n

obs_sdd <- sd(degree(faux.mesa.high, gmode="graph"))
theta <- gden(faux.mesa.high, mode="graph")

brg_sdd <- NULL
cug_sdd <- NULL

for(i in 1:500){
  Ycug <- matrix(0,n,n)
  diag(Ycug)<- NA
  Ycug[upper.tri(Ycug)] <- sample(c(rep(1,203),rep(0,n*(n-1)/2-203)))
  Ycug <- symmetrize(Ycug, rule="upper")
  cug_sdd <- rbind(cug_sdd, sd(degree(Ycug, gmode="graph")))
  
  Ybrg<-matrix(rbinom(n^2,1,theta),n,n)
  diag(Ybrg)<-NA
  Ybrg <- symmetrize(Ybrg,rule="upper")
  brg_sdd <- rbind(brg_sdd, sd(degree(Ybrg, gmode="graph")))
}

CUGcol <- adjustcolor("firebrick", alpha.f=0.7)
BRGcol <- adjustcolor("orangered", alpha.f=0.7)

hist(brg_sdd, col=BRGcol, ylab=NULL, main=NULL, xlab="sd(degree)", xlim=range(brg_sdd,cug_sdd,obs_sdd+.2), ylim=c(0,150))
hist(cug_sdd, col=CUGcol, ylab=NULL, main=NULL, xlab="sd(degree)", add=T)
abline(v=obs_sdd, col="#79AED4", lty=2, lwd=2)
