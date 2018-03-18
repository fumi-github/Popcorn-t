twopop_Zsim = function(h1, h2, gencor, nu) {
  for(trial in 0:99) {
    print(paste0("Trial ",trial));
    foo = read.table(
      paste0("cscore/EAS_EUR.cscore_h2weightbetasim_h_EAS1_EUR1_x0_nuInf_trial", trial),
      header=F, stringsAsFactors=F);
    data$Zpartsim1_0 = NA;
    data$Zpartsim2_0 = NA;
    for (c in 1:22) {
      data$Zpartsim1_0[data$chr==c] = foo$V8[match(data$id[data$chr==c], foo$V3)];
      data$Zpartsim2_0[data$chr==c] = foo$V9[match(data$id[data$chr==c], foo$V3)];
    }
    # S^(2) in Methods:
    data$Zpartsim1_0 =
      qnorm((rank(data$Zpartsim1_0/sqrt(data$score1))-0.5)/nrow(data)) *sqrt(data$score1);
    data$Zpartsim2_0 =
      qnorm((rank(data$Zpartsim2_0/sqrt(data$score2))-0.5)/nrow(data)) *sqrt(data$score2);
    
    foo =
      apply(
        #parApply(cl,   # parallel
        cbind(
          # S^(3) in Methods:
          data$Zpartsim1_0 / sqrt(data$score1),
          data$Zpartsim2_0 / sqrt(data$score2),
          0      * data$scoreX / sqrt(data$score1*data$score2),
          gencor * data$scoreX / sqrt(data$score1*data$score2),
          nu),
        1,
        function (y) {
          # S^(5) in Methods:
          bivnormchangecor(
            # S^(4) in Methods:
            bivnormtomvt(
              y[1:2],
              y[5]),
            y[3], y[4])
        } );
    data$Zpartsim1_0ct = foo[1, ] * sqrt(data$score1);
    data$Zpartsim2_0ct = foo[2, ] * sqrt(data$score2);
    rm(foo);
    
    # S^(6) in Methods:
    data[, paste0("Zsim1_",trial)] = data$Zpartsim1_0ct * sqrt(h1/M);
    data[, paste0("Zsim2_",trial)] = data$Zpartsim2_0ct * sqrt(h2/M);
  }
  data[, c(paste0("Zsim1_",0:99), paste0("Zsim2_",0:99))];
}

twopop_expected_number_of_gwsignificant_loci = function () {
  ### Number of detectable locus
  distance = 500 * 1000;
  result = c();
  for (trial in 0:99) {
    print(paste0("Trial ",trial));
    x = data[, paste0("Zsim1_",trial)];
    x = x / pmax(abs(x),1);
    n = 1*10^5;
    p = gwaspowerncp(n, n*x^2/(1-x^2), 5*10^-8);
    y = (!is.na(p) & p >= 0.05);
    w100K1 = snpstowindow(data$chr[y], data$pos[y], p[y]);
    n = 2*10^5;
    p = gwaspowerncp(n, n*x^2/(1-x^2), 5*10^-8);
    y = (!is.na(p) & p >= 0.05);
    w200K1 = snpstowindow(data$chr[y], data$pos[y], p[y]);
    n = 5*10^5;
    p = gwaspowerncp(n, n*x^2/(1-x^2), 5*10^-8);
    y = (!is.na(p) & p >= 0.05);
    w500K1 = snpstowindow(data$chr[y], data$pos[y], p[y]);
    n = max(data$N1)
    p = gwaspowerncp(n, n*x^2/(1-x^2), 5*10^-8);
    y = (!is.na(p) & p >= 0.05);
    wmaxN1 = snpstowindow(data$chr[y], data$pos[y], p[y]);
    #
    x = data[, paste0("Zsim2_",trial)];
    x = x / pmax(abs(x),1);
    n = 1*10^5;
    p = gwaspowerncp(n, n*x^2/(1-x^2), 5*10^-8);
    y = (!is.na(p) & p >= 0.05);
    w100K2 = snpstowindow(data$chr[y], data$pos[y], p[y]);
    n = 2*10^5;
    p = gwaspowerncp(n, n*x^2/(1-x^2), 5*10^-8);
    y = (!is.na(p) & p >= 0.05);
    w200K2 = snpstowindow(data$chr[y], data$pos[y], p[y]);
    n = 5*10^5;
    p = gwaspowerncp(n, n*x^2/(1-x^2), 5*10^-8);
    y = (!is.na(p) & p >= 0.05);
    w500K2 = snpstowindow(data$chr[y], data$pos[y], p[y]);
    n = max(data$N2);
    p = gwaspowerncp(n, n*x^2/(1-x^2), 5*10^-8);
    y = (!is.na(p) & p >= 0.05);
    wmaxN2 = snpstowindow(data$chr[y], data$pos[y], p[y]);
    #
    result = c(result,
               sum(w100K1$pow), sum(w200K1$pow), sum(w500K1$pow),
               sum(w100K2$pow), sum(w200K2$pow), sum(w500K2$pow),
               windowoverlapcount(w100K1, w100K2),
               windowoverlapcount(w100K1, w200K2),
               windowoverlapcount(w100K1, w500K2),
               windowoverlapcount(w200K1, w100K2),
               windowoverlapcount(w200K1, w200K2),
               windowoverlapcount(w200K1, w500K2),
               windowoverlapcount(w500K1, w100K2),
               windowoverlapcount(w500K1, w200K2),
               windowoverlapcount(w500K1, w500K2),
               sum(wmaxN1$pow), sum(wmaxN2$pow), windowoverlapcount(wmaxN1, wmaxN2));
  }
  result = matrix(result, ncol=18, byrow=TRUE);
  apply(result[,1:15], 2, mean) # no. of loci does not need scaling for M
}
