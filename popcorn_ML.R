### SNP-wise ML

cat('
Estimate heritability parameters by
optim(initial_value, nll, lower, upper, method="L-BFGS-B")

initial lower upper parameter
0.5     0.01  0.99  heritability
0       -0.99 0.99  genetic corelation
2.5     2.01  10    nu for t-distribution (for one population; historical & irrelevant?)
5       2.01  10    nu for t-distribution (for two populations)
1       0.5   1.5   intercept of LD score regression
1       0.5   1.5   lambdaGC of genomic control
')


###
### For one population
### 

# negative log-likelihood
# assuming association Z statistics ~ normal distribution
# modelling intercept of LD score regression
nll_onepop_Znorm_intercept = function(ha) {
  h = ha[1]; # heritability
  a = ha[2]; # intercept
  v = a + (Nmax/M)*h*xkeep;
  #  0.5 * sum(wkeep * (log(v) + ykeep/v)) # ykeep=Z^2, simplified
  #  0.5 * sum(wkeep * (log(v) + ykeep/v + log(2*pi))) # ykeep=Z^2
  #  0.5 * sum(wkeep * (log(v) + ykeep^2/v + log(2*pi))) # ykeep=Z
  sum(wkeep * -dnorm(x=ykeep, sd=sqrt(v), log=TRUE)) # ykeep=Z
}

# negative log-likelihood
# assuming association Z statistics ~ t-distribution
# modelling intercept of LD score regression
# This model converges to nu=10(upper), suggesting normality
nll_onepop_Zt_intercept = function(hna) {
  h = hna[1]; # heritability
  nu = hna[2];
  a = hna[3]; # intercept
  v = a + (Nmax/M)*h*xkeep;
  s = v*(1-2/nu); # divide by the variance of standard t
  sum(wkeep * -log( dt(x=ykeep/sqrt(s), df=nu)/sqrt(s) ))
}

# negative log-likelihood
# assuming allele substitution effect of SNPs ~ normal distribution
# assuming environmental effect ~ normal distribution
# [Not good] nll value can become infinite near boundaries of 0 < h < 1
nll_onepop_snpnorm_envnorm = function(h) { #heritability
  n = length(xkeep);
  m = length(normx);
  v = h*(Nmax/M*xkeep + 1);
  sum(wkeep *
        -log( 
          # dnorm(x=(matrix(ykeep,nrow=n,ncol=m,byrow=F)
          #          -sqrt(1-h)*
          #            matrix(normx,nrow=n,ncol=m,byrow=T)),
          #       sd=matrix(sqrt(v),nrow=n,ncol=m,byrow=F))
          (dnorm(x=(matrix(ykeep,nrow=n,ncol=m,byrow=F)
                    -sqrt(1-h)*
                      matrix(normx,nrow=n,ncol=m,byrow=T))/
                   matrix(sqrt(v),nrow=n,ncol=m,byrow=F)) /
             matrix(sqrt(v),nrow=n,ncol=m,byrow=F))
          %*%
            matrix(normddx)))
}

# negative log-likelihood
# assuming allele substitution effect of SNPs ~ t distribution
# assuming environmental effect ~ normal distribution
weightedloglik_onepop_snpt_envnorm =
  function(x) {
    h = x[1];
    sqrts = x[2];
    y = x[3];
    nu = x[4];
    z = (y - sqrt(1-h)*normx);
    log( 
      sum(
        dt(x=z/sqrts, df=nu)/sqrts *
          normddx))
  }
nll_onepop_snpt_envnorm = function(hn) {
  h = hn[1]; #heritability
  nu = hn[2];
  v = h*(Nmax/M*xkeep + 1);
  s = v*(1-2/nu); # divide by the variance of standard t
  sum(wkeep *
        -parApply(cl,   # parallel
                  cbind(
                    h,
                    sqrt(s),
                    ykeep,
                    nu),
                  1,
                  weightedloglik_onepop_snpt_envnorm)
  )
}

# negative log-likelihood
# assuming allele substitution effect of SNPs ~ t distribution
# assuming environmental effect ~ normal distribution
# modelling intercept of LD score regression
weightedloglik_onepop_snpt_envnorm_intercept =
  function(x) {
    h = x[1];
    sqrts = x[2];
    y = x[3];
    nu = x[4];
    a = x[5]; # intercept
    z = (y - a*sqrt(1-h)*normx);
    log( 
      sum(
        dt(x=z/sqrts, df=nu)/sqrts *
          normddx))
  }
nll_onepop_snpt_envnorm_intercept = function(hna) {
  h = hna[1]; #heritability
  nu = hna[2];
  a = hna[3];
  v = h*(Nmax/M*xkeep + 1);
  s = v*(1-2/nu); # divide by the variance of standard t
  sum(wkeep *
        -parApply(cl,   # parallel
                  cbind(
                    h,
                    sqrt(s),
                    ykeep,
                    nu,
                    a),
                  1,
                  weightedloglik_onepop_snpt_envnorm_intercept)
  )
}

# negative log-likelihood
# assuming allele substitution effect of SNPs ~ t distribution
# assuming environmental effect ~ normal distribution
# modelling GC correction
# [Not good] does not converge as expected
nll_onepop_snpt_envnorm_GC = function(hnl) {
  n = length(xkeep);
  m = length(normx);
  h = hnl[1]; #heritability
  nu = hnl[2];
  lamdaGC = hnl[3];
  v = h*(Nmax/M*xkeep + 1);
  s = v*(1-2/nu); # divide by the variance of standard t
  sum(wkeep *
        -log(
          (dt(x=(matrix(ykeep,nrow=n,ncol=m,byrow=F)*sqrt(lamdaGC)
                 -sqrt(1-h)*
                   matrix(normx,nrow=n,ncol=m,byrow=T))/
                matrix(sqrt(s),nrow=n,ncol=m,byrow=F),
              df=nu) /
             matrix(sqrt(s),nrow=n,ncol=m,byrow=F))
          %*%
            matrix(normddx)))
}


###
### For two populations
### 

# negative log-likelihood
# assuming association Z statistics ~ normal distribution
# modelling intercept of LD score regression
# [Not good] underestimates
nll_twopop_Znorm_intercept = function(gencor, h1, h2, a1, a2) {
  sqrtv1 = sqrt(a1 + (N1max/M)*h1*x1keep);
  sqrtv2 = sqrt(a2 + (N2max/M)*h2*x2keep);
  v = (NXmax/M)*gencor*sqrt(h1)*sqrt(h2)*xXkeep/sqrtv1/sqrtv2; # gencor as genetic-correlation
  # v = gencor; # correlation of explained variance model
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

# negative log-likelihood
# assuming association Z statistics ~ t-distribution
# modelling intercept of LD score regression
# This model converges to nu=10(upper), suggesting normality
nll_twopop_Zt_intercept = function(gencornu, h1, h2, a1, a2) {
  gencor = gencornu[1];
  nu = gencornu[2];
  sqrtv1 = sqrt(a1 + (N1max/M)*h1*x1keep);
  sqrtv2 = sqrt(a2 + (N2max/M)*h2*x2keep);
  v = (NXmax/M)*gencor*sqrt(h1)*sqrt(h2)*xXkeep/sqrtv1/sqrtv2; # gencor as genetic-correlation
  # v = gencor; # correlation of explained variance model
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

# negative log-likelihood
# assuming allele substitution effect of SNPs ~ t distribution
# assuming environmental effect ~ normal distribution
weightedloglik_twopop_snpt_envnorm =
  function(x) {
    sqrtv1 = x[1];
    sqrtv2 = x[2];
    v  = x[3]; # correlation coefficient
    y1 = x[4];
    y2 = x[5];
    nu = x[6];
    h1 = x[7];
    h2 = x[8];
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
# estimate gencor, given unique nu
nll_twopop_snpt_envnorm_givennu = function(gencor, h1, h2, nu) {
  sqrtv1 = sqrt(h1*((N1max/M)*x1keep + 1));
  sqrtv2 = sqrt(h2*((N2max/M)*x2keep + 1));
  v = (NXmax/M)*gencor*sqrt(h1)*sqrt(h2)*xXkeep/sqrtv1/sqrtv2;
  sum(wkeep *
        -parApply(cl,
                  cbind(
                    sqrtv1,
                    sqrtv2,
                    v,
                    y1keep,
                    y2keep,
                    nu, h1, h2 # fixed
                    ),
                  1,
                  weightedloglik_twopop_snpt_envnorm)
  )
}
# estimate unique nu, given gencor
nll_twopop_snpt_envnorm_givengencor = function(nu, h1, h2, gencor) {
  sqrtv1 = sqrt(h1*((N1max/M)*x1keep + 1));
  sqrtv2 = sqrt(h2*((N2max/M)*x2keep + 1));
  v = (NXmax/M)*gencor*sqrt(h1)*sqrt(h2)*xXkeep/sqrtv1/sqrtv2; # fixed
  sum(wkeep *
        -parApply(cl,
                  cbind(
                    sqrtv1,
                    sqrtv2,
                    v, # fixed
                    y1keep,
                    y2keep,
                    nu,
                    h1, h2 # fixed
                    ),
                  1,
                  weightedloglik_twopop_snpt_envnorm)
  )
}
# estimate unique nu and gencor
nll_twopop_snpt_envnorm = function(gencornu, h1, h2) {
  gencor = gencornu[1];
  nu = gencornu[2];
  sqrtv1 = sqrt(h1*((N1max/M)*x1keep + 1));
  sqrtv2 = sqrt(h2*((N2max/M)*x2keep + 1));
  v = (NXmax/M)*gencor*sqrt(h1)*sqrt(h2)*xXkeep/sqrtv1/sqrtv2;
  sum(wkeep *
        -parApply(cl,
                  cbind(
                    sqrtv1,
                    sqrtv2,
                    v,
                    y1keep,
                    y2keep,
                    nu,
                    h1, h2 # fixed
                  ),
                  1,
                  weightedloglik_twopop_snpt_envnorm)
  )
}

# negative log-likelihood
# assuming allele substitution effect of SNPs ~ t distribution
# assuming environmental effect ~ normal distribution
# modelling intercept of LD score regression
weightedloglik_twopop_snpt_envnorm_intercept =
  function(x) {
    sqrtv1 = x[1];
    sqrtv2 = x[2];
    v      = x[3]; # correlation coefficient
    y1     = x[4];
    y2     = x[5];
    nu     = x[6];
    h1     = x[7];
    h2     = x[8];
    a1     = x[9];
    a2     = x[10];
    log(
      sum(
        exp( # -log(2*pi)
          -log((nu-2)/nu) -
            0.5*log(1-v^2) -
            0.5*(nu+2)*log(1 + 1/(nu-2)*(
              0.5/(1+v)*((y1 - a1*sqrt(1-h1)*normx1)/sqrtv1 +
                           (y2 - a2*sqrt(1-h2)*normx2)/sqrtv2)^2 +
                0.5/(1-v)*((y1 - a1*sqrt(1-h1)*normx1)/sqrtv1 - 
                             (y2 - a2*sqrt(1-h2)*normx2)/sqrtv2)^2 ) ) ) *
          normddx1dx2))
  }
# estimate unique nu and gencor
nll_twopop_snpt_envnorm_intercept =
  function(gencornu, h1, h2, a1, a2) {
    gencor = gencornu[1];
    nu = gencornu[2];
    sqrtv1 = sqrt(h1*((N1max/M)*x1keep + 1));
    sqrtv2 = sqrt(h2*((N2max/M)*x2keep + 1));
    v = (NXmax/M)*gencor*sqrt(h1)*sqrt(h2)*xXkeep/sqrtv1/sqrtv2; # gencor as genetic-correlation
    #v = gencor; # correlation of explained variance model
    sum(wkeep *
        -parApply(cl,   # parallel
                  cbind(
                    sqrtv1,
                    sqrtv2,
                    v,
                    y1keep,
                    y2keep,
                    nu,
                    h1, h2, a1, a2 # fixed
                    ),
                  1,
                  weightedloglik_twopop_snpt_envnorm_intercept)
  )
}




### jackknife single
# t + norm model (modelling intercept)
weightedlogliknu =
  function(x) { # global variables: normx normddx
    g = x[1]; # not used
    sqrts = x[2];
    y = x[3];
    nu = x[4];
    c = x[5]; # intercept
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


# modelling intercept
weightedlogliknu =
  function(x) { # global variables: h1 h2 a1 a2 normx1 normx2 normddx1dx2
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
              0.5/(1+v)*((y1 - a1*sqrt(1-h1)*normx1)/sqrtv1 +
                           (y2 - a2*sqrt(1-h2)*normx2)/sqrtv2)^2 +
                0.5/(1-v)*((y1 - a1*sqrt(1-h1)*normx1)/sqrtv1 - 
                             (y2 - a2*sqrt(1-h2)*normx2)/sqrtv2)^2 ) ) ) *
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

h1=0.1069329; h2=0.08556117; a1=1.0349866; a2=0.97422558; #SBP
h1=0.1048264; h2=0.07534258; a1=1.0235582; a2=0.99869954; #T2D
clusterExport(cl, c("h1", "h2", "a1", "a2", "normx1", "normx2", "normddx1dx2"))

wkeep=1/pmax(1,xXunweightedkeep); 
hxjk = jackknifeML()
#
#output = cbind(h12jk, h22jk)
output = hxjk;
write.table(output,
            jackknifefile,
            sep="\t", quote=F, row.names=F)

