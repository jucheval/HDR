require(spatstat)

routine = function(pp, title){
  par(mar=c(0,3,2,0))
  plot(pp, main=title)
  
  par(mar=c(1,4,1,1))
  K = Kest(pp, correction = c("isotropic"))
  plot(K, main="", legend=F, xaxt="n")
  grid()
  
  par(mar=c(4,4,0,1))
  g = pcf(pp, correction = c("isotropic"))
  plot(g, main="", legend=F)
  grid()
}

set.seed(2)
par(mfcol = c(3,3))
# Determinantal Bessel
model = dppBessel(lambda=100, alpha=.05, sigma=0, d=2)
pp = simulate(model, 1)
routine(pp, "repulsive")

# Poisson
pp = rpoispp(100)
routine(pp, "Poisson")

# Clustering
pp = rMatClust(10, 0.2, 10)
routine(pp, "clustered")


