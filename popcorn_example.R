source("popcorn_utilities.R")

# We need to use multiple CPU cores for parallele computation
library(parallel)
cl = makeCluster(detectCores()-1)
clusterExport(cl, c("normx", "normddx", "normx1", "normx2", "normddx1dx2"));


inputfilebase = "~/human/popcorn/lipid/EASLDL_EURLDL.out.dumpdata";
dataunweighted = read.csv(inputfilebase,
                          header=T, stringsAsFactors=F);
dataunweighted = dataunweighted[order(dataunweighted$chr, dataunweighted$pos), ]

inputfile = sub("out", "out_h2weight", inputfilebase);
data = read.csv(inputfile,
                header=T, stringsAsFactors=F);
data = data[order(data$chr, data$pos), ];
data = data[data$N1 >= 0.5*max(data$N1) & data$N2 >= 0.5*max(data$N2), ];
data = data[!is.na(data$Z1), ];

M = 4155885 #Total number of common SNPs in 1000G EAS_EUR shared dataset


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
      c1 = 0.98506231, # intercept for LD score regression
      c2 = 0.9958312,
      nll_twopop_snpt_envnorm_intercept,
      lower=c(-0.99, 2.01),
      upper=c(0.99, 10),
      method="L-BFGS-B");
