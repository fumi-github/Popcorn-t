source("popcorn_utilities.R")
source("popcorn_ML.R")

# We need to use multiple CPU cores for parallel computation
library(parallel)
cl = makeCluster(detectCores()-1)
clusterExport(cl, c("normx", "normddx", "normx1", "normx2", "normddx1dx2"));


###
### Load example data for LDL
###

# Download processed data from
# http://103.253.147.127/PUBLICATIONS/popcorn_EASLDL_EURLDL_dumpdata.zip
# http://www.fumihiko.takeuchi.name//PUBLICATIONS/popcorn_EASLDL_EURLDL_dumpdata.zip
# The original data sources are
# East Asian lipid GWAS by AGEN (PubMed ID 28334899)
# http://blog.nus.edu.sg/agen/summary-statistics/lipids/
# European ancestry lipid GWAS by GLGC (PubMed ID 24097068)
# http://csg.sph.umich.edu/abecasis/public/lipids2013/

inputfilebase = "lipid/EASLDL_EURLDL.out.dumpdata";
dataunweighted = read.csv(inputfilebase,
                          header=T, stringsAsFactors=F);
dataunweighted = dataunweighted[order(dataunweighted$chr, dataunweighted$pos), ];

inputfile = sub("out", "out_h2weight", inputfilebase);
data = read.csv(inputfile,
                header=T, stringsAsFactors=F);
data = data[order(data$chr, data$pos), ];
data = data[data$N1 >= 0.5*max(data$N1) & data$N2 >= 0.5*max(data$N2), ];
data = data[!is.na(data$Z1), ];

M = 4155885; #Total number of common SNPs in 1000G EAS_EUR shared dataset


###
### Estimate heritability and genetic correlation (Fig. 2)
###

# EAS
Nmax = max(data$N1);
keep = (data$score1<=200 & data$score2<=200 & data$scoreX<=200 &
        data$score1>=1   & data$score2>=1   & data$scoreX>=1);
xkeep = (data$score1*(data$N1/Nmax))[keep];
ykeep = data$Z1[keep];
wkeep = 1/pmax(1,
               (dataunweighted$score1[match(data$id, dataunweighted$id)])[keep]);
# The computation below took 34 minutes on a 24 core machine.
# Estimated values should be
# heritability                      = 0.06630128
# nu for t-distribution             = 3.38048012
# intercept for LD score regression = 0.98506231
optim(c(0.5, 2.5, 1),
      nll_onepop_snpt_envnorm_intercept,
      lower=c(0.01, 2.01, 0.5),
      upper=c(0.99, 10, 1.5),
      method="L-BFGS-B");

# EUR
Nmax = max(data$N2);
keep = (data$score1<=200 & data$score2<=200 & data$scoreX<=200 &
        data$score1>=1   & data$score2>=1   & data$scoreX>=1);
xkeep = (data$score2*(data$N2/Nmax))[keep];
ykeep = data$Z2[keep];
wkeep = 1/pmax(1,
               (dataunweighted$score2[match(data$id, dataunweighted$id)])[keep]);
# The computation below took 59 minutes on a 24 core machine.
# Estimated values should be
# heritability                      = 0.1180501
# nu for t-distribution             = 2.8696062
# intercept for LD score regression = 0.9958312
optim(c(0.5, 2.5, 1),
      nll_onepop_snpt_envnorm_intercept,
      lower=c(0.01, 2.01, 0.5),
      upper=c(0.99, 10, 1.5),
      method="L-BFGS-B");

# EAS vs EUR
N1max = max(data$N1);
N2max = max(data$N2);
NXmax = max(sqrt(data$N1*as.numeric(data$N2))); # avoid integer overflow when multiplying
keep = (data$score1<=200 & data$score2<=200 & data$scoreX<=200 &
        data$score1>=1   & data$score2>=1   & data$scoreX>=1);
x1keep = (data$score1*(data$N1/N1max))[keep];
x2keep = (data$score2*(data$N2/N2max))[keep];
xXkeep = (data$scoreX*(sqrt(data$N1*as.numeric(data$N2))/NXmax))[keep];
y1keep = data$Z1[keep];
y2keep = data$Z2[keep];
wkeep = 1/pmax(1,
               (dataunweighted$scoreX[match(data$id, dataunweighted$id)])[keep]);
# The computation below took 19 minutes on a 24 core machine.
# Estimated values should be
# genetic correlation   = 0.7581449
# nu for t-distribution = 3.0750883
optim(c(0, 5),
      h1 = 0.06630128, # heritability
      h2 = 0.1180501,
      a1 = 0.98506231, # intercept for LD score regression
      a2 = 0.9958312,
      nll_twopop_snpt_envnorm_intercept,
      lower=c(-0.99, 2.01),
      upper=c(0.99, 10),
      method="L-BFGS-B");

# EAS (jackknife)
# The computation below takes very long
Nmax = max(data$N1);
keep = (data$score1<=200 & data$score2<=200 & data$scoreX<=200 &
        data$score1>=1   & data$score2>=1   & data$scoreX>=1);
xkeepall = (data$score1*(data$N1/Nmax))[keep];
ykeepall = data$Z1[keep];
wkeepall = 1/pmax(1,
                  (dataunweighted$score1[match(data$id, dataunweighted$id)])[keep]);
njackknife = 200;
hlist = c();
for (i in 1:njackknife) {
  print(i);
  omit = seq(floor((i-1)/njackknife*length(xkeepall)) + 1,
             floor(i/njackknife*length(xkeepall)));
  xkeep = xkeepall[-omit];
  ykeep = ykeepall[-omit];
  wkeep = wkeepall[-omit];
  a1 = 
    optim(c(0.5, 2.5, 1),
          nll_onepop_snpt_envnorm_intercept,
          lower=c(0.01, 2.01, 0.5),
          upper=c(0.99, 10, 1.5),
          method="L-BFGS-B");
  hlist = c(hlist, a1$par);
}
h1jk = matrix(hlist, ncol=3, byrow=TRUE);

# EUR (jackknife)
# The computation below takes very long
Nmax = max(data$N2);
keep = (data$score1<=200 & data$score2<=200 & data$scoreX<=200 &
        data$score1>=1   & data$score2>=1   & data$scoreX>=1);
xkeepall = (data$score2*(data$N2/Nmax))[keep];
ykeepall = data$Z2[keep];
wkeepall = 1/pmax(1,
                  (dataunweighted$score2[match(data$id, dataunweighted$id)])[keep]);
njackknife = 200;
hlist = c();
for (i in 1:njackknife) {
  print(i);
  omit = seq(floor((i-1)/njackknife*length(xkeepall)) + 1,
             floor(i/njackknife*length(xkeepall)));
  xkeep = xkeepall[-omit];
  ykeep = ykeepall[-omit];
  wkeep = wkeepall[-omit];
  a1 = 
    optim(c(0.5, 2.5, 1),
          nll_onepop_snpt_envnorm_intercept,
          lower=c(0.01, 2.01, 0.5),
          upper=c(0.99, 10, 1.5),
          method="L-BFGS-B");
  hlist = c(hlist, a1$par);
}
h2jk = matrix(hlist, ncol=3, byrow=TRUE);

# EAS vs EUR (jackknife)
# The computation below takes very long
N1max = max(data$N1);
N2max = max(data$N2);
NXmax = max(sqrt(data$N1*as.numeric(data$N2))); # avoid integer overflow when multiplying
keep = (data$score1<=200 & data$score2<=200 & data$scoreX<=200 &
          data$score1>=1   & data$score2>=1   & data$scoreX>=1);
x1keepall = (data$score1*(data$N1/N1max))[keep];
x2keepall = (data$score2*(data$N2/N2max))[keep];
xXkeepall = (data$scoreX*(sqrt(data$N1*as.numeric(data$N2))/NXmax))[keep];
y1keepall = data$Z1[keep];
y2keepall = data$Z2[keep];
wkeepall = 1/pmax(1,
                  (dataunweighted$scoreX[match(data$id, dataunweighted$id)])[keep]);
njackknife = 200;
hlist = c();
for (i in 1:njackknife) {
  print(i);
  omit = seq(floor((i-1)/njackknife*length(x1keepall)) + 1,
             floor(i/njackknife*length(x1keepall)));
  x1keep = x1keepall[-omit];
  x2keep = x2keepall[-omit];
  xXkeep = xXkeepall[-omit];
  y1keep = y1keepall[-omit];
  y2keep = y2keepall[-omit];
  wkeep = wkeepall[-omit];
  a1 = 
    optim(c(0, 5),
          h1 = 0.06630128, # heritability
          h2 = 0.1180501,
          a1 = 0.98506231, # intercept for LD score regression
          a2 = 0.9958312,
          nll_twopop_snpt_envnorm_intercept,
          lower=c(-0.99, 2.01),
          upper=c(0.99, 10),
          method="L-BFGS-B");
  hlist = c(hlist, a1$par);
}
hxjk = matrix(hlist, ncol=2, byrow=TRUE);


###
### Power calculation of GWAS under given heritability parameters (Fig. 3)
###

