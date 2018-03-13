source("popcorn_utilities.R")

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
