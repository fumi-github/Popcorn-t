

### SNP-wise ML
# single norm model
nll = function(ag) {
  a = ag[1]; # intercept for stratification
  g = ag[2]; # heritability
  v = a + (Nmax/M)*g*xkeep;
  #  0.5 * sum(wkeep * (log(v) + ykeep/v)) # ykeep=Z^2, simplified
  #  0.5 * sum(wkeep * (log(v) + ykeep/v + log(2*pi))) # ykeep=Z^2
  #  0.5 * sum(wkeep * (log(v) + ykeep^2/v + log(2*pi))) # ykeep=Z
  sum(wkeep * -dnorm(x=ykeep, sd=sqrt(v), log=TRUE)) # ykeep=Z
}
#optim(c(1,0.5), nll)
optim(c(1,0.5), nll, lower=c(0,0.01), upper=c(2,0.99), method="L-BFGS-B")
### RA (EAS vs EUR)
# score1 0.1806062
# score2 0.1286437
# scoreX 0.11246498
0.11246498/sqrt(0.1806062*0.1286437) #0.7378307
### SBP (EAS vs EUR)
# score1 0.08438481
# score2 0.06834101
# scoreX 0.104603500
0.104603500/sqrt(0.08438481*0.06834101) #1.377442
### SBP (JP vs EASwoJP; one population)
# score1 0.09184674
# score2 0.1828195
# scoreX 0.1117416
0.1117416/sqrt(0.09184674*0.1828195) #0.8623265
### SBP (JP vs EASwoJP; one population) gene_effect
# score1 0.08809854
# score2 0.1565064
# scoreX 0.1106435
0.1106435/sqrt(0.08809854*0.1565064) #0.9422703

# single t model
nll = function(agn) {
  a = agn[1]; # intercept for stratification
  g = agn[2]; # heritability
  nu = agn[3];
  v = a + (Nmax/M)*g*xkeep;
  s = v*(1-2/nu); # divide by the variance of standard t
  sum(wkeep * -log( dt(x=ykeep/sqrt(s), df=nu)/sqrt(s) ))
}
optim(c(1,0.5,2.5), nll, lower=c(0,0.01,2.01), upper=c(2,0.99,10), method="L-BFGS-B")
# EASSBP 1.117701  0.096621 10.000000
# EASRA  1.012467  0.201768 10.000000
# EASRA  1.288153  0.236272 fix=4.2
# EURSBP 1.026371  0.254781 10.107019(upper=20)

library(mvtnorm)

# model two distributions, and pre-integrate one
#dx=0.5; normx=seq(-10,10,dx);
dx=0.25; normx=seq(-4,4,dx); #better
normddx = dnorm(normx)*dx;
sum(normddx)
#
sum(dnorm(3 - sqrt(0.7)*normx, sd=sqrt(0.3))*normddx)
dnorm(3)
#
# two dimensinal
normx1 = as.numeric(matrix(normx, nrow=length(normx), ncol=length(normx), byrow=T))
normx2 = as.numeric(matrix(normx, nrow=length(normx), ncol=length(normx), byrow=F))
normddx1dx2 = dnorm(normx1)*dx*dnorm(normx2)*dx;
sum(normddx1dx2)




# norm + norm model (not modelling stratification)
nll = function(g) { #heritability
  n = length(xkeep);
  m = length(normx);
  v = g*(Nmax/M*xkeep + 1);
  sum(wkeep *
        -log( 
          # dnorm(x=(matrix(ykeep,nrow=n,ncol=m,byrow=F)
          #          -sqrt(1-g)*
          #            matrix(normx,nrow=n,ncol=m,byrow=T)),
          #       sd=matrix(sqrt(v),nrow=n,ncol=m,byrow=F))
          (dnorm(x=(matrix(ykeep,nrow=n,ncol=m,byrow=F)
                    -sqrt(1-g)*
                      matrix(normx,nrow=n,ncol=m,byrow=T))/
                   matrix(sqrt(v),nrow=n,ncol=m,byrow=F)) /
             matrix(sqrt(v),nrow=n,ncol=m,byrow=F))
          %*%
            matrix(normddx)))
}
optim(0.5, nll, lower=0.01, upper=0.99, method="L-BFGS-B")



# # t + norm model (not modelling stratification)
# nll = function(gn) {
#   n = length(xkeep);
#   m = length(normx);
#   g = gn[1]; #heritability
#   nu = gn[2];
#   v = g*(Nmax/M*xkeep + 1);
#   s = v*(1-2/nu); # divide by the variance of standard t
#   sum(wkeep *
#         -log(
#           colSums(t(
#           (dt(x=(matrix(ykeep,nrow=n,ncol=m,byrow=F)
#                  -sqrt(1-g)*
#                    matrix(normx,nrow=n,ncol=m,byrow=T))/
# #                matrix(sqrt(s),nrow=n,ncol=m,byrow=F),
#                 sqrt(s),
#               df=nu) /
# #             matrix(sqrt(s),nrow=n,ncol=m,byrow=F))
#               sqrt(s))
# )*normddx)
# #%*%
# #            matrix(normddx)
# ))
# }
# optim(c(0.5,2.5), nll, lower=c(0.01,2.01), upper=c(0.99,10), method="L-BFGS-B")
#optim(0.5, nll, lower=0.01, upper=0.99, method="L-BFGS-B")
#optim(2.5, nll, lower=2.01, upper=10, method="L-BFGS-B")
# EASSBP Zsim1 nu2.5_trial0 0.07161543 3.95857305
# EURSBP Zsim2 nu2.5_trial0 0.07886822 4.50401533
# EASSBP Zsim1 nu3_trial0   0.07240888 5.09570104
# EURSBP Zsim2 nu3_trial0   0.07653706 5.96124694
# EASSBP Zsim1 nu4_trial0   0.07216395 10.00000000
# EURSBP Zsim2 nu4_trial0   0.07808539 10.00000000


# t + norm model (not modelling stratification)
weightedlogliknu =
  function(x) { # global variables: normx normddx
    g = x[1];
    sqrts = x[2];
    y = x[3];
    nu = x[4];
    z = (y - sqrt(1-g)*normx);
    log( 
      sum(
        dt(x=z/sqrts, df=nu)/sqrts *
          normddx))
  }
nll = function(gn) {
  g = gn[1]; #heritability
  nu = gn[2];
  v = g*(Nmax/M*xkeep + 1);
  s = v*(1-2/nu); # divide by the variance of standard t
  sum(wkeep *
        -parApply(cl,   # parallel
                  cbind(
                    g,
                    sqrt(s),
                    ykeep,
                    nu),
                  1,
                  weightedlogliknu)
  )
}
clusterExport(cl, c("normx", "normddx"))
optim(c(0.5,2.5), nll, lower=c(0.01,2.01), upper=c(0.99,10), method="L-BFGS-B")

# t + norm model (modelling constant)
weightedlogliknu =
  function(x) { # global variables: normx normddx
    g = x[1]; # not used
    sqrts = x[2];
    y = x[3];
    nu = x[4];
    c = x[5]; # constant
#    z = (y - c*normx);
    z = (y - c*sqrt(1-g)*normx);
    log( 
      sum(
        dt(x=z/sqrts, df=nu)/sqrts *
          normddx))
  }
nll = function(gnc) {
  g = gnc[1]; #heritability
  nu = gnc[2];
  c = gnc[3];
  v = g*(Nmax/M*xkeep + 1);
  s = v*(1-2/nu); # divide by the variance of standard t
  sum(wkeep *
        -parApply(cl,   # parallel
                  cbind(
                    g,
                    sqrt(s),
                    ykeep,
                    nu,
                    c),
                  1,
                  weightedlogliknu)
  )
}
clusterExport(cl, c("normx", "normddx"))
optim(c(0.5,2.5,1), nll, lower=c(0.01,2.01,0.5), upper=c(0.99,10,1.5), method="L-BFGS-B")



# t + norm model (modelling GC correction)
nll = function(gnl) {
  n = length(xkeep);
  m = length(normx);
  g = gnl[1]; #heritability
  nu = gnl[2];
  lamdaGC = gnl[3];
  v = g*(Nmax/M*xkeep + 1);
  s = v*(1-2/nu); # divide by the variance of standard t
  sum(wkeep *
        -log(
          (dt(x=(matrix(ykeep,nrow=n,ncol=m,byrow=F)*sqrt(lamdaGC)
                 -sqrt(1-g)*
                   matrix(normx,nrow=n,ncol=m,byrow=T))/
                matrix(sqrt(s),nrow=n,ncol=m,byrow=F),
              df=nu) /
             matrix(sqrt(s),nrow=n,ncol=m,byrow=F))
          %*%
            matrix(normddx)))
}
optim(c(0.5,2.5,1), nll, lower=c(0.01,2.01,0.5), upper=c(0.99,10,1.5), method="L-BFGS-B")



# # fit all parameters at once
# nlls = function(ag) {
#   a1 = ag[1]; a2 = ag[2]; aX = ag[3];
#   h1 = ag[4]; h2 = ag[5]; hX = ag[6]; #heritability
#   v1 = a1 + (Ns[[1]]/M)*h1*xkeeps[[1]];
#   v2 = a2 + (Ns[[2]]/M)*h2*xkeeps[[2]];
#   vX = aX + (Ns[[3]]/M)*sqrt(h1*h2)*hX*xkeeps[[3]];
#   0.5 * sum(wkeeps[[1]] * (log(v1) + ykeeps[[1]]/v1)) +
#     0.5 * sum(wkeeps[[2]] * (log(v2) + ykeeps[[2]]/v2)) +
#     0.5 * sum(wkeeps[[3]] * (log(vX) + ykeeps[[3]]/vX))
# }
# optim(c(0.1,0.1,0.1, 0.1,0.1,0.5), nlls) #not converge well
# optim(c(0.1,0.1,0.1, 0.1,0.1,0.5), nlls,
#       lower=c(0,0,0, 0.01,0.01,0.01),
#       upper=c(2,2,2, 1,1,1), method="L-BFGS-B")
# #0.09225807 0.18312765 0.85959073




#### corr coeff

# single norm model; underestimates
a1=1; h1=0.1;
a2=1; h2=0.1;
sqrtv1 = sqrt(a1 + (N1max/M)*h1*x1keep);
sqrtv2 = sqrt(a2 + (N2max/M)*h2*x2keep);
nll = function(r) {
  #  v = (NXmax/M)*r*sqrt(h1)*sqrt(h2)*xXkeep/sqrtv1/sqrtv2; # genetic-correlation
  v = r; # correlation of explained variance model
  sum(wkeep *
        #    apply(cbind(y1keep/sqrtv1, y2keep/sqrtv2, v),
        #          1,
        #          function(x){-dmvnorm(x[1:2],
        #                            sigma=matrix(c(1,x[3],x[3],1), ncol=2),
        #                            log=TRUE)})
        (0.5*log(1-v^2) +
           0.25/(1+v)*(y1keep/sqrtv1 + y2keep/sqrtv2)^2 +
           0.25/(1-v)*(y1keep/sqrtv1 - y2keep/sqrtv2)^2)
  )
}
optim(0, nll, lower=-0.99, upper=0.99, method="L-BFGS-B")

# single t model; underestimates
a1=1; h1=0.1;
a2=1; h2=0.1;
nu=4.5
sqrtv1 = sqrt(a1 + (N1max/M)*h1*x1keep);
sqrtv2 = sqrt(a2 + (N2max/M)*h2*x2keep);
nll = function(r) {
  #  v = (NXmax/M)*r*sqrt(h1)*sqrt(h2)*xXkeep/sqrtv1/sqrtv2; # genetic-correlation
  v = r; # correlation of explained variance model
  sum(wkeep *
        # apply(cbind(y1keep/sqrtv1, y2keep/sqrtv2, v),
        #       1,
        #       function(x){-dmvt(x[1:2],
        #         sigma=matrix(c(1,x[3],x[3],1), ncol=2)*(1-2/nu),
        #         df=nu,
        #         log=TRUE)})
        (0.5*log(1-v^2) +
           0.5*(nu+2)*log(1 + 1/(nu-2)*(
             0.5/(1+v)*(y1keep/sqrtv1 + y2keep/sqrtv2)^2 +
               0.5/(1-v)*(y1keep/sqrtv1 - y2keep/sqrtv2)^2 ) ) )
  )
}
optim(0, nll, lower=-0.99, upper=0.99, method="L-BFGS-B")

# unique nu as for dmvt
weightedloglik =
  function(x) { # global variables: h1 h2 nu normx1 normx2 normddx1dx2
    sqrtv1 = x[1];
    sqrtv2 = x[2];
    v  = x[3]; # correlation coefficient
    y1 = x[4];
    y2 = x[5];
    log(
      sum(
        exp( # -log(2*pi)
          -log((nu-2)/nu) -
          0.5*log(1-v^2) -
              0.5*(nu+2)*log(1 + 1/(nu-2)*(
                0.5/(1+v)*((y1 - sqrt(1-h1)*normx1)/sqrtv1 +
                             (y2 - sqrt(1-h2)*normx2)/sqrtv2)^2 +
                  0.5/(1-v)*((y1 - sqrt(1-h1)*normx1)/sqrtv1 - 
                               (y2 - sqrt(1-h2)*normx2)/sqrtv2)^2 ) ) ) *
          normddx1dx2))
  }
# not modeling stratification
weightedlogliknu =
  function(x) { # global variables: h1 h2 normx1 normx2 normddx1dx2
    sqrtv1 = x[1];
    sqrtv2 = x[2];
    v  = x[3]; # correlation coefficient
    y1 = x[4];
    y2 = x[5];
    nu = x[6];
    log(
      sum(
        exp( # -log(2*pi)
          -log((nu-2)/nu) -
            0.5*log(1-v^2) -
            0.5*(nu+2)*log(1 + 1/(nu-2)*(
              0.5/(1+v)*((y1 - sqrt(1-h1)*normx1)/sqrtv1 +
                           (y2 - sqrt(1-h2)*normx2)/sqrtv2)^2 +
                0.5/(1-v)*((y1 - sqrt(1-h1)*normx1)/sqrtv1 - 
                             (y2 - sqrt(1-h2)*normx2)/sqrtv2)^2 ) ) ) *
          normddx1dx2))
  }
# modelling constant
weightedlogliknu =
  function(x) { # global variables: h1 h2 c1 c2 normx1 normx2 normddx1dx2
    sqrtv1 = x[1];
    sqrtv2 = x[2];
    v  = x[3]; # correlation coefficient
    y1 = x[4];
    y2 = x[5];
    nu = x[6];
    log(
      sum(
        exp( # -log(2*pi)
          -log((nu-2)/nu) -
            0.5*log(1-v^2) -
            0.5*(nu+2)*log(1 + 1/(nu-2)*(
              0.5/(1+v)*((y1 - c1*sqrt(1-h1)*normx1)/sqrtv1 +
                           (y2 - c2*sqrt(1-h2)*normx2)/sqrtv2)^2 +
                0.5/(1-v)*((y1 - c1*sqrt(1-h1)*normx1)/sqrtv1 - 
                             (y2 - c2*sqrt(1-h2)*normx2)/sqrtv2)^2 ) ) ) *
          normddx1dx2))
  }

source("~/Documents/R/bivt/bivt.R", chdir=TRUE)
# two separte nu's (nu1 nu2) as for dbivt
weightedloglik =
  function(x) { # global variables: h1 h2 nu1 nu2 normx1 normx2 normddx1dx2
    sqrtv1 = x[1];
    sqrtv2 = x[2];
    v  = x[3]; # correlation coefficient
    y1 = x[4];
    y2 = x[5];
    z1 = 
        (y1 - sqrt(1-h1)*normx1)/sqrtv1;
    z2 = 
        (y2 - sqrt(1-h2)*normx2)/sqrtv2;
    log( 
      sum(
        dbivt(cbind(z1,z2), v, nu1, nu2) *
        normddx1dx2))
  }

# two separte NORMALs
weightedloglik =
  function(x) { # global variables: h1 h2 normx1 normx2 normddx1dx2
    sqrtv1 = x[1];
    sqrtv2 = x[2];
    v  = x[3]; # correlation coefficient
    y1 = x[4];
    y2 = x[5];
    z1 = (y1 - sqrt(1-h1)*normx1)/sqrtv1;
    z2 = (y2 - sqrt(1-h2)*normx2)/sqrtv2;
    log( 
      sum(
        exp(-(0.5*log(1-v^2) +
                0.25/(1+v)*(z1 + z2)^2 +
                0.25/(1-v)*(z1 - z2)^2)) *
          normddx1dx2))
  }



# OPTION: cache function
# library(memoise)
# forget(weightedloglik);
# weightedloglik = memoise(weightedloglik)

# OPTION: parallel
library(parallel)
cl = makeCluster(detectCores()-1)
# stopCluster(cl)

# t + norm (unique nu or two nu's)
lambdaGC1=1; lambdaGC2=1;
nll = function(r) {
  sqrtv1 = sqrt(h1*((N1max/M)*x1keep + 1));
  sqrtv2 = sqrt(h2*((N2max/M)*x2keep + 1));
  v = (NXmax/M)*r*sqrt(h1)*sqrt(h2)*xXkeep/sqrtv1/sqrtv2; # r as genetic-correlation
  #v = r; # correlation of explained variance model
  sum(wkeep *
        #      -apply(
        -parApply(cl,   # parallel
                  cbind(
                    sqrtv1, # round(sqrtv1, digits=2),   # round for cache
                    sqrtv2,
                    v,
                    y1keep*sqrt(lambdaGC1),
                    y2keep*sqrt(lambdaGC2)),
                  1,
                  weightedloglik)
  )
}
clusterExport(cl, c("h1", "h2", "nu", "normx1", "normx2", "normddx1dx2"))
clusterExport(cl, c("h1", "h2", "nu1", "nu2", "normx1", "normx2", "normddx1dx2"))
clusterExport(cl, c("h1", "h2", "nu1", "nu2", "normx1", "normx2", "normddx1dx2",
                    "dbivt", "nu1nucombtonu2"))
#needs finite values of 'fn'
optim(0, nll, lower=-0.99, upper=0.99, method="L-BFGS-B") 
# parabolic interpolation; Brent of optim; not fast; sometimes LESS ACCURATE THAN optim
optimize(nll, lower=-0.99, upper=0.99) 
# Nelder-Mead; allows fn=Inf; robust but slow
# optim(0, nll) # AVOID; sometimes caught in local optima near 0

# t + norm; infer unique nu
lambdaGC1=1; lambdaGC2=1;
nll = function(nu) {
  sqrtv1 = sqrt(h1*((N1max/M)*x1keep + 1));
  sqrtv2 = sqrt(h2*((N2max/M)*x2keep + 1));
  v = (NXmax/M)*gencor*sqrt(h1)*sqrt(h2)*xXkeep/sqrtv1/sqrtv2; # gencor as genetic-correlation
  #v = r; # correlation of explained variance model
  sum(wkeep *
        #      -apply(
        -parApply(cl,   # parallel
                  cbind(
                    sqrtv1, # round(sqrtv1, digits=2),   # round for cache
                    sqrtv2,
                    v,
                    y1keep*sqrt(lambdaGC1),
                    y2keep*sqrt(lambdaGC2),
                    nu),
                  1,
                  weightedlogliknu)
  )
}
clusterExport(cl, c("h1", "h2", "gencor", "normx1", "normx2", "normddx1dx2"))
optimize(nll, lower=2.01, upper=10)

# t + norm; infer unique nu and gencor
lambdaGC1=1; lambdaGC2=1;
nll = function(nugencor) {
  nu = nugencor[1];
  gencor = nugencor[2];
  sqrtv1 = sqrt(h1*((N1max/M)*x1keep + 1));
  sqrtv2 = sqrt(h2*((N2max/M)*x2keep + 1));
  v = (NXmax/M)*gencor*sqrt(h1)*sqrt(h2)*xXkeep/sqrtv1/sqrtv2; # gencor as genetic-correlation
  #v = r; # correlation of explained variance model
  sum(wkeep *
        #      -apply(
        -parApply(cl,   # parallel
                  cbind(
                    sqrtv1, # round(sqrtv1, digits=2),   # round for cache
                    sqrtv2,
                    v,
                    y1keep*sqrt(lambdaGC1),
                    y2keep*sqrt(lambdaGC2),
                    nu),
                  1,
                  weightedlogliknu)
  )
}
clusterExport(cl, c("h1", "h2", "normx1", "normx2", "normddx1dx2"))
clusterExport(cl, c("h1", "h2", "c1", "c2", "normx1", "normx2", "normddx1dx2"))
#needs finite values of 'fn'
optim(c(5,0), nll, lower=c(2.01,-0.99), upper=c(10,0.99), method="L-BFGS-B") 
# Nelder-Mead; allows fn=Inf; robust but slow
# optim(0, nll) # AVOID; sometimes caught in local optima near 0




### jackknife single
# t + norm model (modelling constant)
weightedlogliknu =
  function(x) { # global variables: normx normddx
    g = x[1]; # not used
    sqrts = x[2];
    y = x[3];
    nu = x[4];
    c = x[5]; # constant
    #    z = (y - c*normx);
    z = (y - c*sqrt(1-g)*normx);
    log( 
      sum(
        dt(x=z/sqrts, df=nu)/sqrts *
          normddx))
  }

jackknifeML = function () {
  njackknife = 200
  hlist = c()
  for (i in 1:njackknife) {
    print(i);
    omit = seq(floor((i-1)/njackknife*length(xkeep)) + 1,
               floor(i/njackknife*length(xkeep)));
    xsampled = xkeep[-omit];
    ysampled = ykeep[-omit];
    wsampled = wkeep[-omit];
    
    nll = function(gnc) {
      g = gnc[1]; #heritability
      nu = gnc[2];
      c = gnc[3];
      v = g*(Nmax/M*xsampled + 1);
      s = v*(1-2/nu); # divide by the variance of standard t
      sum(wsampled *
            -parApply(cl,   # parallel
                      cbind(
                        g,
                        sqrt(s),
                        ysampled,
                        nu,
                        c),
                      1,
                      weightedlogliknu)
      )
    }
    a1= optim(c(0.5,2.5,1), nll, lower=c(0.01,2.01,0.5), upper=c(0.99,10,1.5), method="L-BFGS-B")
    hlist = c(hlist, a1$par)
  }
  matrix(hlist, ncol=3, byrow=TRUE)
}

clusterExport(cl, c("normx", "normddx"))

Nmax=N1max; xkeep=x1keep; wkeep=1/pmax(1,x1unweightedkeep); ykeep=y1keep
h12jk = jackknifeML()
#
Nmax=N2max; xkeep=x2keep; wkeep=1/pmax(1,x2unweightedkeep); ykeep=y2keep
h22jk = jackknifeML()


# modelling constant
weightedlogliknu =
  function(x) { # global variables: h1 h2 c1 c2 normx1 normx2 normddx1dx2
    sqrtv1 = x[1];
    sqrtv2 = x[2];
    v  = x[3]; # correlation coefficient
    y1 = x[4];
    y2 = x[5];
    nu = x[6];
    log(
      sum(
        exp( # -log(2*pi)
          -log((nu-2)/nu) -
            0.5*log(1-v^2) -
            0.5*(nu+2)*log(1 + 1/(nu-2)*(
              0.5/(1+v)*((y1 - c1*sqrt(1-h1)*normx1)/sqrtv1 +
                           (y2 - c2*sqrt(1-h2)*normx2)/sqrtv2)^2 +
                0.5/(1-v)*((y1 - c1*sqrt(1-h1)*normx1)/sqrtv1 - 
                             (y2 - c2*sqrt(1-h2)*normx2)/sqrtv2)^2 ) ) ) *
          normddx1dx2))
  }

jackknifeML = function () {
  njackknife = 200
  hlist = c()
  for (i in 1:njackknife) {
    print(i);
    omit = seq(floor((i-1)/njackknife*length(x1keep)) + 1,
               floor(i/njackknife*length(x1keep)));
    x1sampled = x1keep[-omit];
    x2sampled = x2keep[-omit];
    xXsampled = xXkeep[-omit];
    y1sampled = y1keep[-omit];
    y2sampled = y2keep[-omit];
    wsampled = wkeep[-omit];

    # t + norm; infer unique nu and gencor
    lambdaGC1=1; lambdaGC2=1;
    nll = function(nugencor) {
      nu = nugencor[1];
      gencor = nugencor[2];
      sqrtv1 = sqrt(h1*((N1max/M)*x1sampled + 1));
      sqrtv2 = sqrt(h2*((N2max/M)*x2sampled + 1));
      v = (NXmax/M)*gencor*sqrt(h1)*sqrt(h2)*xXsampled/sqrtv1/sqrtv2; # gencor as genetic-correlation
      #v = r; # correlation of explained variance model
      sum(wsampled *
            #      -apply(
            -parApply(cl,   # parallel
                      cbind(
                        sqrtv1, # round(sqrtv1, digits=2),   # round for cache
                        sqrtv2,
                        v,
                        y1sampled*sqrt(lambdaGC1),
                        y2sampled*sqrt(lambdaGC2),
                        nu),
                      1,
                      weightedlogliknu)
      )
    }
    a1 = optim(c(5,0), nll, lower=c(2.01,-0.99), upper=c(10,0.99), method="L-BFGS-B") 
    
    hlist = c(hlist, a1$par)
  }
  matrix(hlist, ncol=2, byrow=TRUE)
}

h1=0.1069329; h2=0.08556117; c1=1.0349866; c2=0.97422558; #SBP
h1=0.1048264; h2=0.07534258; c1=1.0235582; c2=0.99869954; #T2D
clusterExport(cl, c("h1", "h2", "c1", "c2", "normx1", "normx2", "normddx1dx2"))

wkeep=1/pmax(1,xXunweightedkeep); 
hxjk = jackknifeML()
#
#output = cbind(h12jk, h22jk)
output = hxjk;
write.table(output,
            jackknifefile,
            sep="\t", quote=F, row.names=F)

