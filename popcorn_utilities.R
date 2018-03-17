# randomize the effect allele before computing var, cov, cor
library(moments);
myvar = function(x, ...) {
  set.seed(123);
  d=rbinom(length(x), 1, 0.5)*2-1;
  var(d*x, ...)
}
mycov = function(x, y, ...) {
  set.seed(123);
  d=rbinom(length(x), 1, 0.5)*2-1;
  cov(d*x, d*y, ...)
}
mycor = function(x, y, ...) {
  set.seed(123);
  d=rbinom(length(x), 1, 0.5)*2-1;
  cor(d*x, d*y, ...)
}
mykurtosis = function(x, ...) {
  set.seed(123);
  d=rbinom(length(x), 1, 0.5)*2-1;
  kurtosis(d*x, ...)
}

# To count the number of loci
snpstowindow = function(chr, pos, pow, distance=5e5) {
  if (length(chr)==0) {
    return(data.frame(row.names=c("chr", "start", "end", "pow")))
  }
  # reorder by physical position
  x = order(chr, pos);
  chr=chr[x];
  pos=pos[x];
  pow=pow[x];
  rm(x);

  result = c(chr[1], pos[1]);
  cprev = chr[1];
  pprev = pos[1];
  powmax  = pow[1];
  if (length(chr) > 1) {
    for (i in 2:length(chr)) {
      if (chr[i]==cprev & (pos[i]-pprev <= distance)) {
        pprev = pos[i];
        powmax = max(powmax, pow[i]);
      } else {
        result = c(result, pprev, powmax, chr[i], pos[i]);
        cprev = chr[i];
        pprev = pos[i];
        powmax = pow[i];
      }
    }
  }
  result = c(result, pprev, powmax);
  result = data.frame(matrix(result, ncol=4, byrow=TRUE));
  names(result) = c("chr", "start", "end", "pow");
  result
}
windowoverlapcount = function(w1, w2) {
  c = 0;
  for (i in 1:nrow(w1)) {
    x = (w2$chr==w1$chr[i]) &
        !((w2$start > w1$end[i]) | (w2$end < w1$start[i]));
    c = c + w1$pow[i] * sum(w2$pow[x]);
  }
  c
}

# n is sample size
gwaspowerncp = 
  Vectorize(
    function (n, ncp, alpha) {
      if (ncp==Inf) {
        1
      } else {
        threshold = qf(alpha, 1, n-2, lower.tail=F);
        pf(threshold, 1, n-2, ncp=ncp, lower.tail=F);
      }
    }
  )

# For numerical integration
dx = 0.25;
normx = seq(-4, 4, dx);
normddx = dnorm(normx)*dx;
normx1 = as.numeric(matrix(normx, nrow=length(normx), ncol=length(normx), byrow=T))
normx2 = as.numeric(matrix(normx, nrow=length(normx), ncol=length(normx), byrow=F))
normddx1dx2 = dnorm(normx1)*dx*dnorm(normx2)*dx;
rm(dx);
lockBinding(c("normx", "normddx", "normx1", "normx2", "normddx1dx2"), .GlobalEnv);

# Transform from a bivariate normal distribution with covariance matrix
# ((1, r), (r, 1)) to one with the covariance matrix ((1, s), (s, 1))
# by quantile normalization.
# Absolute values of r, s should be smaller than one.
# TODO simplify the code by using math
bivnormchangecor =
  function (x, r, s) {
    cosX1 = cos(0.25*pi-0.5*acos(r));
    cosX2 = cos(0.25*pi+0.5*acos(r));
    sinX1 = sin(0.25*pi-0.5*acos(r));
    sinX2 = sin(0.25*pi+0.5*acos(r));
    cosY1 = cos(0.25*pi-0.5*acos(s));
    cosY2 = cos(0.25*pi+0.5*acos(s));
    sinY1 = sin(0.25*pi-0.5*acos(s));
    sinY2 = sin(0.25*pi+0.5*acos(s));
    if (! is.matrix(x)) {
      x = matrix(x, nrow=1)
    }
    U = x %*%
      matrix(c( sinX2, -cosX2,
                -sinX1,  cosX1) /
               (cosX1*sinX2-sinX1*cosX2),
             nrow=2, byrow=T);
    U %*%
      matrix(c(cosY1, cosY2,
               sinY1, sinY2),
             nrow=2, byrow=T)
  }

# For no correlation; ie var = ((1,0),(0,1))
bivnormtomvt =
  function (x, df) {
    if (! is.matrix(x)) {
      x = matrix(x, nrow=1)
    }
    rpre = sqrt(apply(x^2, 1, sum));
    rpost = sqrt((df-2)*(exp(rpre^2/df)-1));
    x/rpre*rpost
  }
