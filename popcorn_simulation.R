# generates data$Zsim1 data$Zsim2

# xpre = x used in cscore_h2weightbetasim simulation
# gencor = target genetic correlation
# !!! Mall = M for .cscore; scale as if all 1000G SNPs were included in GWAS
# Version 1; h1, h2, gencor computed by simple LD score regression
h1=0.1108609; h2=0.08472976; xpre=0; gencor=0.694; nu=4.3; M=4155885; Mall=4174819 #SBPEAS_SBPEUR
h1=0.1128007; h2=h1;         xpre=0; gencor=0;     nu=4.3; M=5148923; Mall=5510924 #SBPEAS_self
h1=0.09498145; h2=h1;        xpre=0; gencor=0;     nu=4.3; M=5942349; Mall=5976091 #SBPEUR_self
lx=0.04; ly=0.02

h1=0.09028171; h2=0.08525658; xpre=0; gencor=0.661; nu=3.9; M=4155885; Mall=4174819 #DBPEAS_DBPEUR
h1=0.08642003; h2=h1;         xpre=0; gencor=0;     nu=3.9; M=5148923; Mall=5510924 #DBPEAS_self
h1=0.09043663; h2=h1;        xpre=0; gencor=0;      nu=3.9; M=5942349; Mall=5976091 #DBPEUR_self
lx=0.03; ly=0.03

h1=0.2623927; h2=0.2113607; xpre=0; gencor=0.787; nu=4.3; M=4155885; Mall=4174819 #HEIGHTEAS_HEIGHTEURMEN
h1=0.27436; h2=h1;          xpre=0; gencor=0;     nu=4.3; M=5148923; Mall=5510924 #HEIGHTEAS_self
h1=0.2551615; h2=h1;        xpre=0; gencor=0;     nu=4.3; M=5942349; Mall=5976091 #HEIGHTEURMEN_self
lx=0.04; ly=0.03

h1=0.0970393; h2=0.439972; xpre=0; gencor=0.852; nu=2.8; M=4155885; Mall=4174819 #HDLEAS_HDLEUR
h1=0.07659177; h2=h1;         xpre=0; gencor=0;  nu=2.8; M=5148923; Mall=5510924 #HDLEAS_self
h1=0.5614308; h2=h1;        xpre=0; gencor=0;    nu=2.8; M=5942349; Mall=5976091 #HDLEUR_self
lx=0.04; ly=0.07

h1=0.04721333; h2=0.2913462; xpre=0; gencor=0.426; nu=2.8; M=4155885; Mall=4174819 #LDLEAS_LDLEUR
h1=0.05704457; h2=h1;         xpre=0; gencor=0;  nu=2.8; M=5148923; Mall=5510924 #LDLEAS_self
h1=0.4047947; h2=h1;        xpre=0; gencor=0;    nu=2.8; M=5942349; Mall=5976091 #LDLEUR_self
lx=0.03; ly=0.05

h1=0.07661929; h2=0.349348; xpre=0; gencor=0.488; nu=2.7; M=4155885; Mall=4174819 #TCEAS_TCEUR
h1=0.0887034;  h2=h1;       xpre=0; gencor=0;     nu=2.7; M=5148923; Mall=5510924 #TCEAS_self
h1=0.4910303;  h2=h1;       xpre=0; gencor=0;     nu=2.7; M=5942349; Mall=5976091 #TCEUR_self
lx=0.04; ly=0.06

h1=0.0900469; h2=0.5052914; xpre=0; gencor=1.034; nu=2.7; M=4155885; Mall=4174819 #TGEAS_TGEUR
h1=0.1083463;  h2=h1;       xpre=0; gencor=0;     nu=2.7; M=5148923; Mall=5510924 #TGEAS_self
h1=0.5985967;  h2=h1;       xpre=0; gencor=0;     nu=2.7; M=5942349; Mall=5976091 #TGEUR_self
lx=0.04; ly=0.07

h1=0.4933613; h2=0.0816908; xpre=0; gencor=0.629; nu=8.6; M=4155885; Mall=4174819 #T2DEAS_T2DEUR
h1=0.4893035;  h2=h1;       xpre=0; gencor=0;     nu=8.6; M=5148923; Mall=5510924 #T2DEAS_self
# h1=0.08736977;  h2=h1;       xpre=0; gencor=0;     nu=8.6; M=5942349; Mall=5976091 #T2DEUR_self; overwritten by nu=4.9 below
lx=0.04; ly=0.03

h1=0.1562091; h2=0.07742799; xpre=0; gencor=0.913; nu=4.9; M=4155885; Mall=4174819 #T2DEAS-2_T2DEUR
h1=0.147656;  h2=h1;       xpre=0; gencor=0;       nu=4.9; M=5148923; Mall=5510924 #T2DEAS-2_self
h1=0.08736977;  h2=h1;       xpre=0; gencor=0;     nu=4.9; M=5942349; Mall=5976091 #T2DEUR_self
lx=0.03; ly=0.04


h1=0.1687944; h2=0.1504502; xpre=0; gencor=0.943; nu=5.4; M=4155885; Mall=4174819 #BMIEAS_BMIMEN
h1=0.1751875; h2=h1;          xpre=0; gencor=0;     nu=5.4; M=5148923; Mall=5510924 #BMIEAS_self
h1=0.1786419; h2=h1;        xpre=0; gencor=0;     nu=5.4; M=5942349; Mall=5976091 #BMIEURMEN_self
lx=0.04; ly=0.05


# Version 2; h1, h2 prefit by t+norm in each pop
R --no-save >> `hostname`.log << EOF
h1=0.1297765; h2=0.06321855; xpre=0; gencor=0.9409613; nu=4.2; M=4155885; Mall=4174819 #SBPEAS_SBPEUR
inputfilebase="~/human/popcorn/BP/EASSBP_EURSBP.out.dumpdata"; source("~/Downloads/template_EAS_EUR.R");
.04; ly=.03;
EOF
R --no-save >> `hostname`.log << EOF
h1=0.1311458; h2=h1;         xpre=0; gencor=0;     nu=4.7; M=5148923; Mall=5510924 #SBPEAS_self
inputfilebase="~/human/popcorn/BP/EASSBP_self.out.dumpdata"; source("~/Downloads/template_EAS_self.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.06357095; h2=h1;        xpre=0; gencor=0;     nu=3.5; M=5942349; Mall=5976091 #SBPEUR_self
inputfilebase="~/human/popcorn/BP/EURSBP_self.out.dumpdata"; source("~/Downloads/template_EUR_self.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.1053324; h2=0.06394449; xpre=0; gencor=0.9331359; nu=3.8; M=4155885; Mall=4174819 #DBPEAS_DBPEUR
inputfilebase="~/human/popcorn/BP/EASDBP_EURDBP.out.dumpdata"; source("~/Downloads/template_EAS_EUR.R");
.04; ly=.04;
EOF
R --no-save >> `hostname`.log << EOF
h1=0.1041399; h2=h1;         xpre=0; gencor=0;     nu=4.4; M=5148923; Mall=5510924 #DBPEAS_self
inputfilebase="~/human/popcorn/BP/EASDBP_self.out.dumpdata"; source("~/Downloads/template_EAS_self.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.06270845; h2=h1;        xpre=0; gencor=0;      nu=3.3; M=5942349; Mall=5976091 #DBPEUR_self
inputfilebase="~/human/popcorn/BP/EURDBP_self.out.dumpdata"; source("~/Downloads/template_EUR_self.R");
EOF

R --no-save >> `hostname`.log << EOF
h1=0.2375312; h2=0.1432211; xpre=0; gencor=0.99; nu=4.4; M=4155885; Mall=4174819 #HEIGHTEAS_HEIGHTEURMEN
inputfilebase="~/human/popcorn/GIANT/EASHEIGHT_EURHEIGHTMEN.out.dumpdata"; source("~/Downloads/template_EAS_EUR.R");
.04; ly=.06
EOF
R --no-save >> `hostname`.log << EOF
h1=0.2441181; h2=h1;          xpre=0; gencor=0;     nu=6.6; M=5148923; Mall=5510924 #HEIGHTEAS_self
inputfilebase="~/human/popcorn/GIANT/EASHEIGHT_self.out.dumpdata"; source("~/Downloads/template_EAS_self.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.1469808; h2=h1;        xpre=0; gencor=0;     nu=2.5; M=5942349; Mall=5976091 #HEIGHTEURMEN_self
inputfilebase="~/human/popcorn/GIANT/EURHEIGHTMEN_self.out.dumpdata"; source("~/Downloads/template_EUR_self.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.09227932; h2=0.1593731; xpre=0; gencor=0.99; nu=3.2; M=4155885; Mall=4174819 #HDLEAS_HDLEUR
inputfilebase="~/human/popcorn/lipid/EASHDL_EURHDL.out.dumpdata"; source("~/Downloads/template_EAS_EUR.R");
.03; ly=.05;
EOF
R --no-save >> `hostname`.log << EOF
h1=0.09180711; h2=h1;         xpre=0; gencor=0;  nu=2.6; M=5148923; Mall=5510924 #HDLEAS_self
inputfilebase="~/human/popcorn/lipid/EASHDL_self.out.dumpdata"; source("~/Downloads/template_EAS_self.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.1770591; h2=h1;        xpre=0; gencor=0;    nu=3.2; M=5942349; Mall=5976091 #HDLEUR_self
inputfilebase="~/human/popcorn/lipid/EURHDL_self.out.dumpdata"; source("~/Downloads/template_EUR_self.R");
EOF

R --no-save >> `hostname`.log << EOF
h1=0.06190259; h2=0.1169615; xpre=0; gencor=0.8315008; nu=3.0; M=4155885; Mall=4174819 #LDLEAS_LDLEUR
inputfilebase="~/human/popcorn/lipid/EASLDL_EURLDL.out.dumpdata"; source("~/Downloads/template_EAS_EUR.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.06250484; h2=h1;         xpre=0; gencor=0;  nu=3.1; M=5148923; Mall=5510924 #LDLEAS_self
inputfilebase="~/human/popcorn/lipid/EASLDL_self.out.dumpdata"; source("~/Downloads/template_EAS_self.R");
.03; ly=.06;
EOF
R --no-save >> `hostname`.log << EOF
h1=0.1238256; h2=h1;        xpre=0; gencor=0;    nu=2.8; M=5942349; Mall=5976091 #LDLEUR_self
inputfilebase="~/human/popcorn/lipid/EURLDL_self.out.dumpdata"; source("~/Downloads/template_EUR_self.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.08135579; h2=0.1473798; xpre=0; gencor=0.8184684; nu=2.9; M=4155885; Mall=4174819 #TCEAS_TCEUR
inputfilebase="~/human/popcorn/lipid/EASTC_EURTC.out.dumpdata"; source("~/Downloads/template_EAS_EUR.R");
.03; ly==.06;
EOF
R --no-save >> `hostname`.log << EOF
h1=0.08014178;  h2=h1;       xpre=0; gencor=0;     nu=2.8; M=5148923; Mall=5510924 #TCEAS_self
inputfilebase="~/human/popcorn/lipid/EASTC_self.out.dumpdata"; source("~/Downloads/template_EAS_self.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.1581742;  h2=h1;       xpre=0; gencor=0;     nu=2.8; M=5942349; Mall=5976091 #TCEUR_self
inputfilebase="~/human/popcorn/lipid/EURTC_self.out.dumpdata"; source("~/Downloads/template_EUR_self.R");
EOF

R --no-save >> `hostname`.log << EOF
h1=0.07174658; h2=0.1211473; xpre=0; gencor=0.99; nu=3.2; M=4155885; Mall=4174819 #TGEAS_TGEUR
inputfilebase="~/human/popcorn/lipid/EASTG_EURTG.out.dumpdata"; source("~/Downloads/template_EAS_EUR.R");
.04; ly=.05;
EOF
R --no-save >> `hostname`.log << EOF
h1=0.07071624;  h2=h1;       xpre=0; gencor=0;     nu=2.9; M=5148923; Mall=5510924 #TGEAS_self
inputfilebase="~/human/popcorn/lipid/EASTG_self.out.dumpdata"; source("~/Downloads/template_EAS_self.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.1337642;  h2=h1;       xpre=0; gencor=0;     nu=3.1; M=5942349; Mall=5976091 #TGEUR_self
inputfilebase="~/human/popcorn/lipid/EURTG_self.out.dumpdata"; source("~/Downloads/template_EUR_self.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.1998018; h2=0.07448747; xpre=0; gencor=0.7770905; nu=5.3; M=4155885; Mall=4174819 #T2DEAS-2_T2DEUR
inputfilebase="~/human/popcorn/T2D/EAST2D-2_EURT2D.out.dumpdata"; source("~/Downloads/template_EAS_EUR.R");
.03; ly=.04;
EOF
R --no-save >> `hostname`.log << EOF
h1=0.2013361;  h2=h1;       xpre=0; gencor=0;       nu=7.1; M=5148923; Mall=5510924 #T2DEAS-2_self
inputfilebase="~/human/popcorn/T2D/EAST2D-2_self.out.dumpdata"; source("~/Downloads/template_EAS_self.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.07563739;  h2=h1;       xpre=0; gencor=0;     nu=3.6; M=5942349; Mall=5976091 #T2DEUR_self
inputfilebase="~/human/popcorn/T2D/EURT2D_self.out.dumpdata"; source("~/Downloads/template_EUR_self.R");
EOF

R --no-save >> `hostname`.log << EOF
h1=0.1992564; h2=0.09055735; xpre=0; gencor=0.99; nu=5.5; M=4155885; Mall=4174819 #BMIEAS_BMIMEN
inputfilebase="~/human/popcorn/GIANT/EASBMI_EURBMIMEN.out.dumpdata"; source("~/Downloads/template_EAS_EUR.R");
.04; ly=.06;
EOF
R --no-save >> `hostname`.log << EOF
h1=0.2039724; h2=h1;          xpre=0; gencor=0;     nu=6.2; M=5148923; Mall=5510924 #BMIEAS_self
inputfilebase="~/human/popcorn/GIANT/EASBMI_self.out.dumpdata"; source("~/Downloads/template_EAS_self.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.1013778; h2=h1;        xpre=0; gencor=0;     nu=2.4; M=5942349; Mall=5976091 #BMIEURMEN_self
inputfilebase="~/human/popcorn/GIANT/EURBMIMEN_self.out.dumpdata"; source("~/Downloads/template_EUR_self.R");
EOF




# Version 3; h1, h2, c1, c2 prefit by t+norm in each pop
R --no-save >> `hostname`.log << EOF
h1=0.1069329; h2=0.08556117; c1=1.0349866; c2=0.97422558; xpre=0; gencor=0.898; nu=4.0; M=4155885; Mall=4174819 #SBPEAS_SBPEUR
lx=.03; ly=.04;
inputfilebase="~/human/popcorn/BP/EASSBP_EURSBP.out.dumpdata"; source("~/Downloads/t_loadsim_EAS_EUR.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.1092447; c1=1.0283512; h2=h1; c2=c1;        xpre=0; gencor=0;     nu=4.1; M=5148923; Mall=5510924 #SBPEAS_self
lx=.03; ly=.04;
inputfilebase="~/human/popcorn/BP/EASSBP_self.out.dumpdata"; source("~/Downloads/t_loadsim_EAS_self.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.09709136; c1=0.96724951; h2=h1; c2=c1;        xpre=0; gencor=0;     nu=4.7; M=5942349; Mall=5976091 #SBPEUR_self
lx=.03; ly=.04;
inputfilebase="~/human/popcorn/BP/EURSBP_self.out.dumpdata"; source("~/Downloads/t_loadsim_EUR_self.R");
EOF

R --no-save >> `hostname`.log << EOF
h1=0.08966287; h2=0.09098336; c1=1.02455281; c2=0.9675149; xpre=0; gencor=0.851; nu=3.9; M=4155885; Mall=4174819 #DBPEAS_DBPEUR
lx=.03; ly=.04;
inputfilebase="~/human/popcorn/BP/EASDBP_EURDBP.out.dumpdata"; source("~/Downloads/t_loadsim_EAS_EUR.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.0871617; c1=1.0221491; h2=h1; c2=c1;         xpre=0; gencor=0;     nu=3.8; M=5148923; Mall=5510924 #DBPEAS_self
lx=.03; ly=.04;
inputfilebase="~/human/popcorn/BP/EASDBP_self.out.dumpdata"; source("~/Downloads/t_loadsim_EAS_self.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.09660772; c1=0.96581063; h2=h1; c2=c1;        xpre=0; gencor=0;      nu=4.5; M=5942349; Mall=5976091 #DBPEUR_self
lx=.03; ly=.04;
inputfilebase="~/human/popcorn/BP/EURDBP_self.out.dumpdata"; source("~/Downloads/t_loadsim_EUR_self.R");
EOF

R --no-save >> `hostname`.log << EOF
h1=0.08504717; h2=0.1548046; c1=0.99138531; c2=1.0107311; xpre=0; gencor=0.99; nu=3.1; M=4155885; Mall=4174819 #HDLEAS_HDLEUR
lx=.03; ly=.06;
inputfilebase="~/human/popcorn/lipid/EASHDL_EURHDL.out.dumpdata"; source("~/Downloads/t_loadsim_EAS_EUR.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.08484241; c1=0.99099684; h2=h1; c2=c1;         xpre=0; gencor=0;  nu=2.8; M=5148923; Mall=5510924 #HDLEAS_self
lx=.03; ly=.06;
inputfilebase="~/human/popcorn/lipid/EASHDL_self.out.dumpdata"; source("~/Downloads/t_loadsim_EAS_self.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.1699546; c1=1.0144683; h2=h1; c2=c1;        xpre=0; gencor=0;    nu=3.0; M=5942349; Mall=5976091 #HDLEUR_self
lx=.03; ly=.06;
inputfilebase="~/human/popcorn/lipid/EURHDL_self.out.dumpdata"; source("~/Downloads/t_loadsim_EUR_self.R");
EOF

R --no-save >> `hostname`.log << EOF
h1=0.06630128; h2=0.1180501; c1=0.98506231; c2=0.9958312; xpre=0; gencor=0.758; nu=3.1; M=4155885; Mall=4174819 #LDLEAS_LDLEUR
lx=.03; ly=.05;
inputfilebase="~/human/popcorn/lipid/EASLDL_EURLDL.out.dumpdata"; source("~/Downloads/t_loadsim_EAS_EUR.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.06703234; c1=0.9835395; h2=h1; c2=c1;         xpre=0; gencor=0;  nu=3.4; M=5148923; Mall=5510924 #LDLEAS_self
lx=.03; ly=.05;
inputfilebase="~/human/popcorn/lipid/EASLDL_self.out.dumpdata"; source("~/Downloads/t_loadsim_EAS_self.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.1243268; c1=0.9985278; h2=h1; c2=c1;        xpre=0; gencor=0;    nu=2.9; M=5942349; Mall=5976091 #LDLEUR_self
lx=.03; ly=.05;
inputfilebase="~/human/popcorn/lipid/EURLDL_self.out.dumpdata"; source("~/Downloads/t_loadsim_EUR_self.R");
EOF

R --no-save >> `hostname`.log << EOF
h1=0.08702871; h2=0.1476327; c1=0.98754701; c2=1.0044978; xpre=0; gencor=0.773; nu=2.9; M=4155885; Mall=4174819 #TCEAS_TCEUR
lx=.03; ly=.06;
inputfilebase="~/human/popcorn/lipid/EASTC_EURTC.out.dumpdata"; source("~/Downloads/t_loadsim_EAS_EUR.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.08323712; c1=0.98371431; h2=h1; c2=c1;       xpre=0; gencor=0;     nu=3.3; M=5148923; Mall=5510924 #TCEAS_self
lx=.03; ly=.06;
inputfilebase="~/human/popcorn/lipid/EASTC_self.out.dumpdata"; source("~/Downloads/t_loadsim_EAS_self.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.1580177; c1=1.0075465; h2=h1; c2=c1;       xpre=0; gencor=0;     nu=2.7; M=5942349; Mall=5976091 #TCEUR_self
lx=.03; ly=.06;
inputfilebase="~/human/popcorn/lipid/EURTC_self.out.dumpdata"; source("~/Downloads/t_loadsim_EUR_self.R");
EOF

R --no-save >> `hostname`.log << EOF
h1=0.069759; h2=0.1253272; c1=0.98974; c2=0.9927937; xpre=0; gencor=0.99; nu=3.4; M=4155885; Mall=4174819 #TGEAS_TGEUR
lx=.03; ly=.05;
inputfilebase="~/human/popcorn/lipid/EASTG_EURTG.out.dumpdata"; source("~/Downloads/t_loadsim_EAS_EUR.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.06932404; c1=0.98777252;  h2=h1; c2=c1;       xpre=0; gencor=0;     nu=3.2; M=5148923; Mall=5510924 #TGEAS_self
lx=.03; ly=.05;
inputfilebase="~/human/popcorn/lipid/EASTG_self.out.dumpdata"; source("~/Downloads/t_loadsim_EAS_self.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.1398853; c1=0.9910802; h2=h1; c2=c1;       xpre=0; gencor=0;     nu=3.3; M=5942349; Mall=5976091 #TGEUR_self
lx=.03; ly=.05;
inputfilebase="~/human/popcorn/lipid/EURTG_self.out.dumpdata"; source("~/Downloads/t_loadsim_EUR_self.R");
EOF

R --no-save >> `hostname`.log << EOF
h1=0.1048264; h2=0.07534258; c1=1.0235582; c2=0.99869954; xpre=0; gencor=0.99; nu=4.1; M=4155885; Mall=4174819 #T2DEAS-2_T2DEUR
lx=.03; ly=.04;
inputfilebase="~/human/popcorn/T2D/EAST2D-2_EURT2D.out.dumpdata"; source("~/Downloads/t_loadsim_EAS_EUR.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.08964393; c1=1.02533038;  h2=h1; c2=c1;       xpre=0; gencor=0;       nu=4.1; M=5148923; Mall=5510924 #T2DEAS-2_self
lx=.03; ly=.04;
inputfilebase="~/human/popcorn/T2D/EAST2D-2_self.out.dumpdata"; source("~/Downloads/t_loadsim_EAS_self.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.07903342; c1=0.99525918;  h2=h1; c2=c1;       xpre=0; gencor=0;     nu=3.7; M=5942349; Mall=5976091 #T2DEUR_self
lx=.03; ly=.04;
inputfilebase="~/human/popcorn/T2D/EURT2D_self.out.dumpdata"; source("~/Downloads/t_loadsim_EUR_self.R");
EOF

# R --no-save >> `hostname`.log << EOF
# h1=0.2162307; h2=0.1894674; c1=1.0091078; c2=0.9224786; xpre=0; gencor=0.929; nu=5.2; M=4155885; Mall=4174819 #HEIGHTEAS_HEIGHTEURMEN
# lx=.04; ly=.05;
# inputfilebase="~/human/popcorn/GIANT/EASHEIGHT_EURHEIGHTMEN.out.dumpdata"; source("~/Downloads/t_loadsim_EAS_EUR.R");
# EOF
R --no-save >> `hostname`.log << EOF
h1=0.2247666; c1=1.0077092; h2=h1; c2=c1;          xpre=0; gencor=0;     nu=6.1; M=5148923; Mall=5510924 #HEIGHTEAS_self
lx=.04; ly=.07;
inputfilebase="~/human/popcorn/GIANT/EASHEIGHT_self.out.dumpdata"; source("~/Downloads/t_loadsim_EAS_self.R");
EOF
# R --no-save >> `hostname`.log << EOF
# h1=0.2151893; c1=0.9095801; h2=h1; c2=c1;        xpre=0; gencor=0;     nu=5.3; M=5942349; Mall=5976091 #HEIGHTEURMEN_self
# lx=.04; ly=.05;
# inputfilebase="~/human/popcorn/GIANT/EURHEIGHTMEN_self.out.dumpdata"; source("~/Downloads/t_loadsim_EUR_self.R");
# EOF
R --no-save >> `hostname`.log << EOF
h1=0.2172695; h2=0.3060659; c1=1.0087438; c2=1.2066581; xpre=0; gencor=0.935; nu=3.7; M=4155885; Mall=4174819 #HEIGHTEAS_HEIGHTEUR2014
lx=.04; ly=.07
inputfilebase="~/human/popcorn/GIANT/EASHEIGHT_EURHEIGHT2014.out.dumpdata"; source("~/Downloads/t_loadsim_EAS_EUR.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.3332631; c1=1.1753883; h2=h1; c2=c1;        xpre=0; gencor=0;     nu=3.6; M=5942349; Mall=5976091 #HEIGHTEUR2014_self
lx=.04; ly=.07;
inputfilebase="~/human/popcorn/GIANT/EURHEIGHT2014_self.out.dumpdata"; source("~/Downloads/t_loadsim_EUR_self.R");
EOF


# R --no-save >> `hostname`.log << EOF
# h1=0.161419; h2=0.1286135; c1=1.053353; c2=0.9358289; xpre=0; gencor=0.99; nu=5.2; M=4155885; Mall=4174819 #BMIEAS_BMIMEN
# lx=.04; ly=.04;
# inputfilebase="~/human/popcorn/GIANT/EASBMI_EURBMIMEN.out.dumpdata"; source("~/Downloads/t_loadsim_EAS_EUR.R");
# EOF
R --no-save >> `hostname`.log << EOF
h1=0.1670337; c1=1.0501513; h2=h1; c2=c1;          xpre=0; gencor=0;     nu=5.2; M=5148923; Mall=5510924 #BMIEAS_self
lx=.03; ly=.04;
inputfilebase="~/human/popcorn/GIANT/EASBMI_self.out.dumpdata"; source("~/Downloads/t_loadsim_EAS_self.R");
EOF
# R --no-save >> `hostname`.log << EOF
# h1=0.1469824; c1=0.9300321; h2=h1; c2=c1;        xpre=0; gencor=0;     nu=5.8; M=5942349; Mall=5976091 #BMIEURMEN_self
# lx=.04; ly=.04;
# inputfilebase="~/human/popcorn/GIANT/EURBMIMEN_self.out.dumpdata"; source("~/Downloads/t_loadsim_EUR_self.R");
# EOF
R --no-save >> `hostname`.log << EOF
h1=0.1609931; h2=0.1143785; c1=1.0539514; c2=0.8073608; xpre=0; gencor=0.941; nu=5.1; M=4155885; Mall=4174819 #BMIEAS_BMI2015
lx=.03; ly=.04
inputfilebase="~/human/popcorn/GIANT/EASBMI_EURBMI2015.out.dumpdata"; source("~/Downloads/t_loadsim_EAS_EUR.R");
EOF
R --no-save >> `hostname`.log << EOF
h1=0.1274812; c1=0.792683; h2=h1; c2=c1;        xpre=0; gencor=0;     nu=5.6; M=5942349; Mall=5976091 #BMIEUR2015_self
lx=.03; ly=.04;
inputfilebase="~/human/popcorn/GIANT/EURBMI2015_self.out.dumpdata"; source("~/Downloads/t_loadsim_EUR_self.R");
EOF


### dull theoretical simulation; one GWAS per one SNP
library(mvtnorm) # r is not genetic correlation
r=0.65; nu=3.2; trial=0;
r=0.65; nu=5; trial=0;
hxv = r*data$scoreX/sqrt(data$score1+M/data$N1)/sqrt(data$score2+M/data$N2)
set.seed(trial)
source("~/Documents/R/bivt/bivt.R", chdir=TRUE)
clusterExport(cl, c("nu", "nu1nucombtonu2", "rbivt"))
betasim = t(
#  sapply(
  parSapply(cl, 
  hxv,
  function(x){
  # mvt formulation in mvtnorm library
    rmvt(1,
         sigma=matrix(c(1,x,x,1)*(1-2/nu),nrow=2),
         df=nu)
  # # bivt formulation in bivt.R
  #   rbivt(1,
  #         x,
  #         df1=nu, df2=nu)
  }))

# betasim = rmvt(nrow(data),
#                sigma=matrix(c(1,r,r,1)*(1-2/nu),nrow=2),
#                df=nu)
# betasim = rmvnorm(nrow(data),
#                sigma=matrix(c(1,r,hx,1),nrow=2))
data$Zpartsim1 = betasim[,1]
data$Zpartsim2 = betasim[,2]
data$rcg1_0 = rnorm(nrow(data))
data$rcg2_0 = rnorm(nrow(data))
h1=0.1; h2=0.1;
# dull theoretical single t
data$Zsim1 = sqrt(h1/M*data$N1*data$score1 + 1)*
                data$Zpartsim1
data$Zsim2 = sqrt(h2/M*data$N2*data$score2 + 1)*
                data$Zpartsim2
# dull theoretical t + norm
data$Zsim1 = (sqrt(h1/M*(data$N1*data$score1 + M))*
                data$Zpartsim1 +
                sqrt(1-h1)*data$rcg1_0)
data$Zsim2 = (sqrt(h2/M*(data$N2*data$score2 + M))*
                data$Zpartsim2 +
                sqrt(1-h2)*data$rcg2_0)


### 1000G based simluation

### randomlycombinedgenotype
for (pn in 1:2) {
  p = c("EAS","EUR")[pn];
  for (c in 1:22) {
    print(paste(p, c))
    bim = read.table(
      paste0("~/human/1000G/phase3_shapeit2/",p,".chr",c,".maf001.snv.bim"),
      header=F, stringsAsFactors=F)
    for (trial in 0:5){
      foo = read.table(
        paste0("~/human/1000G/phase3_shapeit2/",p,
               ".chr",c,
               ".maf001.snv.randomlycombinedgenotype_trial",trial),
        header=F)
      if(nrow(bim) != nrow(foo)) { print("ERROR"); break }
      data[data$chr==c, paste0("rcg", pn, "_", trial)] =
        foo$V1[match(data[data$chr==c, "id"], bim$V2)]
    }
  }
}
cov(data[,c("rcg1_0","rcg1_1","rcg2_0","rcg2_1")])
cor(data$rcg1_0[-1], data$rcg1_0[-nrow(data)]) # neighbouring SNPs correlated
cor(data$rcg1_0[1:1000000], data$rcg1_0[100001:1100000]) # far SNPs not correlated

trial=0;
## two populations; EAS_EUR
foo = read.table(
  paste0("~/human/popcorn/cscore/EAS_EUR.cscore_h2weightbetasim_h_EAS1_EUR1_x0_nuInf_trial",trial),
  header=F, stringsAsFactors=F)
data$Zpartsim1_0 = NA;
data$Zpartsim2_0 = NA;
for (c in 1:22) {
  data$Zpartsim1_0[data$chr==c] = foo$V8[match(data$id[data$chr==c], foo$V3)]
  data$Zpartsim2_0[data$chr==c] = foo$V9[match(data$id[data$chr==c], foo$V3)]
}
# foo = read.table(
#   paste0("~/human/popcorn/cscore/EAS_EUR.cscore_h2weightbetasim_h_EAS1_EUR1_x0_nuInf_trial",trial+1),
#   header=F, stringsAsFactors=F)
# data$Zpartsim1_1 = NA;
# data$Zpartsim2_1 = NA;
# for (c in 1:22) {
#   data$Zpartsim1_1[data$chr==c] = foo$V8[match(data$id[data$chr==c], foo$V3)]
#   data$Zpartsim2_1[data$chr==c] = foo$V9[match(data$id[data$chr==c], foo$V3)]
# }
## one population; EAS_self, EUR_self
foo = read.table(
  # without ldcutoff05, noise from SNPs in weak LD is added, resulting in lighter
  #"~/human/popcorn/cscore/EAS.cscore_h2weightbetasim_ldcutoff05_h_EAS1_EUR1_x0.65_nu2.3_trial0",
  #"~/human/popcorn/cscore/EAS.cscore_h2weightbetasim_ldcutoff05_h_EAS1_EUR1_x0.65_nu2.5_trial0",
  #"~/human/popcorn/cscore/EAS.cscore_h2weightbetasim_ldcutoff05_h_EAS1_EUR1_x0.65_nu2.7_trial0",
  #"~/human/popcorn/cscore/EAS.cscore_h2weightbetasim_ldcutoff05_h_EAS1_EUR1_x0.65_nu3_trial0",
  #"~/human/popcorn/cscore/EAS.cscore_h2weightbetasim_ldcutoff05_h_EAS1_EUR1_x0.65_nu4_trial0",
  #"~/human/popcorn/cscore/EAS.cscore_h2weightbetasim_ldcutoff05_h_EAS1_EUR1_x0.65_nuInf_trial0",
  #"~/human/popcorn/cscore/EAS.cscore_h2weightbetasim_ldcutoff0316_h_EAS1_EUR1_x0.65_nuInf_trial0",
  #"~/human/popcorn/cscore/EAS.cscore_h2weightbetasim_h_EAS1_EUR1_x0.65_nuInf_trial0",
  #"~/human/popcorn/cscore/EAS.cscore_h2weightbetasim_h_EAS1_EUR1_x0.75_nuInf_trial0",
  #"~/human/popcorn/cscore/EAS.cscore_h2weightbetasim_h_EAS1_EUR1_x0.9_nuInf_trial0",
  paste0("~/human/popcorn/cscore/EAS.cscore_h2weightbetasim_h_EAS1_EUR1_x0_nuInf_trial",trial),
  #"~/human/popcorn/cscore/EUR.cscore_h2weightbetasim_ldcutoff05_h_EAS1_EUR1_x0.65_nu2.3_trial0",
  #"~/human/popcorn/cscore/EUR.cscore_h2weightbetasim_ldcutoff05_h_EAS1_EUR1_x0.65_nu2.5_trial0",
  #"~/human/popcorn/cscore/EUR.cscore_h2weightbetasim_ldcutoff05_h_EAS1_EUR1_x0.65_nu2.7_trial0",
  #"~/human/popcorn/cscore/EUR.cscore_h2weightbetasim_ldcutoff05_h_EAS1_EUR1_x0.65_nu3_trial0",
  #"~/human/popcorn/cscore/EUR.cscore_h2weightbetasim_ldcutoff05_h_EAS1_EUR1_x0.65_nu4_trial0",
  #"~/human/popcorn/cscore/EUR.cscore_h2weightbetasim_ldcutoff05_h_EAS1_EUR1_x0.65_nuInf_trial0",
  #"~/human/popcorn/cscore/EUR.cscore_h2weightbetasim_ldcutoff0316_h_EAS1_EUR1_x0.65_nuInf_trial0",
  #"~/human/popcorn/cscore/EUR.cscore_h2weightbetasim_h_EAS1_EUR1_x0.65_nuInf_trial0",
  #"~/human/popcorn/cscore/EUR.cscore_h2weightbetasim_h_EAS1_EUR1_x0.75_nuInf_trial0",
  #"~/human/popcorn/cscore/EUR.cscore_h2weightbetasim_h_EAS1_EUR1_x0.9_nuInf_trial0",
  #paste0("~/human/popcorn/cscore/EUR.cscore_h2weightbetasim_h_EAS1_EUR1_x0_nuInf_trial",trial),
  header=F, stringsAsFactors=F)
data$Zpartsim1_0 = NA
for (c in 1:22) {
  data$Zpartsim1_0[data$chr==c] = foo$V7[match(data$id[data$chr==c], foo$V3)]
}
data$Zpartsim2_0 = data$Zpartsim1_0 #dummy
# foo = read.table(
#   #paste0("~/human/popcorn/cscore/EAS.cscore_h2weightbetasim_h_EAS1_EUR1_x0_nuInf_trial",trial+1),
#   paste0("~/human/popcorn/cscore/EUR.cscore_h2weightbetasim_h_EAS1_EUR1_x0_nuInf_trial",trial+1),
#   header=F, stringsAsFactors=F)
# data$Zpartsim1_1 = NA
# for (c in 1:22) {
#   data$Zpartsim1_1[data$chr==c] = foo$V7[match(data$id[data$chr==c], foo$V3)]
# }
# data$Zpartsim2_1 = data$Zpartsim1_1 #dummy


# INSPECT
cov(data[, c("Zpartsim1","Zpartsim2")])
cor(data$Zpartsim1[-1], data$Zpartsim1[-nrow(data)]) # neighbouring SNPs correlated
cor(data$Zpartsim1[1:1000000], data$Zpartsim1[100001:1100000]) # far SNPs not correlated
sd(data$Zpartsim1/sqrt(data$score1)) # ~1
sd(data$Zpartsim2/sqrt(data$score2)) # ~1
summary(lm((Zpartsim1)^2 ~ score1, data=data))
summary(lm((Zpartsim2)^2 ~ score2, data=data))
summary(lm(Zpartsim1*Zpartsim2 ~ scoreX, data=data))
summary(lm(log((data$Zpartsim1)^2) ~ log(data$score1), data=data)) #coeff=1 => linear
summary(lm(log((data$Zpartsim2)^2) ~ log(data$score2), data=data)) #coeff=1 => linear
summary(lm(log(abs(Zpartsim1*Zpartsim2)) ~ log(abs(scoreX)), data=data))
#
library(ggplot2)
ggplot(data, aes(Zpartsim1/sqrt(score1))) + geom_density()
ggplot(data, aes(Zpartsim2/sqrt(score2))) + geom_density()
ggplot(data, aes(factor(floor(score1/200)), Zpartsim1)) + geom_violin()
ggplot(data, aes(factor(floor(score2/200)), Zpartsim2)) + geom_violin()
ggplot(data, aes(factor(floor(scoreX/200)), Zpartsim1*Zpartsim2)) + geom_violin()
ggplot(data, aes(log10(score1), log10(Zpartsim1^2))) + geom_bin2d()
ggplot(data, aes(log10(score2), log10(Zpartsim2^2))) + geom_bin2d()
ggplot(data, aes(log10(abs(scoreX)), log10(abs(Zpartsim1*Zpartsim2)))) + geom_bin2d()
#
y1 = data$Zpartsim1;
y2 = data$Zpartsim2;
sqrtv1 = sqrt(data$score1);
sqrtv2 = sqrt(data$score2);
nll = function(r) {
  v = r*data$scoreX/(sqrtv1*sqrtv2);
  sum(
        #    apply(cbind(y1/sqrtv1, y2/sqrtv2, v),
        #          1,
        #          function(x){-dmvnorm(x[1:2],
        #                            sigma=matrix(c(1,x[3],x[3],1), ncol=2),
        #                            log=TRUE)})
        (0.5*log(1-v^2) +
           0.25/(1+v)*(y1/sqrtv1 + y2/sqrtv2)^2 +
           0.25/(1-v)*(y1/sqrtv1 - y2/sqrtv2)^2)
  )
}
optim(0, nll, lower=-0.99, upper=0.99, method="L-BFGS-B")
# Drastically diminishes (due to accumulation of weak LD pairs)
# h_EAS1_EUR1_x0.9_nuInf_trial0  0.3854008
# h_EAS1_EUR1_x0.65_nuInf_trial0 0.2668497
# ldcutoff0316_h_EAS1_EUR1_x0.65_nuInf_trial0  0.3927572
# ldcutoff05_h_EAS1_EUR1_x0.65_nuInf_trial0    0.5301982
optimize(nll, lower=-0.99, upper=0.99)
cor(y1/sqrtv1, y2/sqrtv2)
cor(y1/sqrt(data$score1 -3), y2/sqrt(data$score2 -3))


# data$Zpartsim1 = data$Zpartsim1/sd(data$Zpartsim1/sqrt(data$score1))
# data$Zpartsim2 = data$Zpartsim2/sd(data$Zpartsim2/sqrt(data$score2))
#
# score1predict = predict(lm((Zpartsim1)^2 ~ score1, data=data))
# data$Zpartsim1 = data$Zpartsim1*sqrt(data$score1/pmax(data$score1, score1predict))
# score2predict = predict(lm((Zpartsim2)^2 ~ score2, data=data))
# data$Zpartsim2 = data$Zpartsim2*sqrt(data$score2/pmax(data$score2, score2predict))
# quantile normalize; change is suttle
data$Zpartsim1_0 =
  qnorm((rank(data$Zpartsim1_0/sqrt(data$score1))-0.5)/nrow(data)) *sqrt(data$score1);
data$Zpartsim2_0 =
  qnorm((rank(data$Zpartsim2_0/sqrt(data$score2))-0.5)/nrow(data)) *sqrt(data$score2);
# data$Zpartsim1_1 =
#   qnorm((rank(data$Zpartsim1_1/sqrt(data$score1))-0.5)/nrow(data)) *sqrt(data$score1);
# data$Zpartsim2_1 =
#   qnorm((rank(data$Zpartsim2_1/sqrt(data$score2))-0.5)/nrow(data)) *sqrt(data$score2);


if (gencor==0) {
  # This part is for one population
  data$Zpartsim1_0c = data$Zpartsim1_0;
  data$Zpartsim2_0c = data$Zpartsim2_0;
  data$Zpartsim1_0ct =
    qt(pnorm(data$Zpartsim1_0c/sqrt(data$score1)), df=nu) *
    sqrt(1-2/nu)*sqrt(data$score1);
  data$Zpartsim2_0ct =
    qt(pnorm(data$Zpartsim2_0c/sqrt(data$score2)), df=nu) *
    sqrt(1-2/nu)*sqrt(data$score2);
  # data$Zpartsim1_1c = data$Zpartsim1_1;
  # data$Zpartsim2_1c = data$Zpartsim2_1;
} else {
  source("~/Documents/R/bivt/bivt.R", chdir=TRUE);
  clusterExport(cl, c("bivnormchangecor","bivnormtomvt","nu"));
  foo =
    #apply(
    parApply(cl,   # parallel
             cbind(
               data$Zpartsim1_0 / sqrt(data$score1),
               data$Zpartsim2_0 / sqrt(data$score2),
               xpre * data$scoreX / sqrt(data$score1*data$score2),
               gencor * data$scoreX / sqrt(data$score1*data$score2)),
             1,
             function (y) {
               bivnormchangecor(
                 bivnormtomvt(
                   y[1:2],
                 nu),
                 y[3], y[4])
             } );
  data$Zpartsim1_0ct = foo[1, ] * sqrt(data$score1);
  data$Zpartsim2_0ct = foo[2, ] * sqrt(data$score2);
  # foo =
  #   #apply(
  #   parApply(cl,   # parallel
  #            cbind(
  #              data$Zpartsim1_1 / sqrt(data$score1),
  #              data$Zpartsim2_1 / sqrt(data$score2),
  #              xpre * data$scoreX / sqrt(data$score1*data$score2),
  #              gencor * data$scoreX / sqrt(data$score1*data$score2)),
  #            1,
  #            function (y) {
  #              bivnormchangecor(
  #                y[1:2], y[3], y[4])
  #            } );
  # data$Zpartsim1_1c = foo[1, ] * sqrt(data$score1);
  # data$Zpartsim2_1c = foo[2, ] * sqrt(data$score2);
  rm(foo);
}


# construct t-distribution by dividing with chisq
# transtot1 =
#   sqrt(nu-2)/
#   sqrt(qchisq((rank((data$Zpartsim1_1c/sqrt(data$score1))^2)-0.5)/nrow(data), df=nu));
# transtot2 =
#   sqrt(nu-2)/
#   sqrt(qchisq((rank((data$Zpartsim2_1c/sqrt(data$score2))^2)-0.5)/nrow(data), df=nu));
# transtot12 =
#   sqrt(nu-2)/
#   sqrt(qchisq((rank( (
#       (data$Zpartsim1_1/sqrt(data$score1)+
#          data$Zpartsim2_1/sqrt(data$score2))/sqrt(2))^2 )-0.5)/nrow(data), df=nu));
# set.seed(trial)
# transtot12 =
#   sqrt(nu-2)/
#   sqrt(rchisq(nrow(data), df=nu));


# data$Zpartsim1_0ct = 
#   data$Zpartsim1_0c *transtot12
# data$Zpartsim2_0ct = 
#   data$Zpartsim2_0c *transtot12


### Escape here to record true effect
data[, paste0("Zsim1_",trial)] = data$Zpartsim1_0ct * sqrt(h1/M);
data[, paste0("Zsim2_",trial)] = data$Zpartsim2_0ct * sqrt(h2/M);
# EAS_EUR
output = data[, c("id","chr","pos","a1isREF",
                  paste0("Zsim1_",trial),
                  paste0("Zsim2_",trial))]
# EAS
output = data[, c("id","chr","pos","EURMAF","a1isREF",
                  paste0("Zsim1_",trial))]
# EUR
output = data[, c("id","chr","pos","EASMAF","a1isREF",
                  paste0("Zsim1_",trial))]
#
write.table(output,
  paste0(sub("dumpdata","trueeffect_trial",inputfile), trial),
  row.names=FALSE, quote=FALSE, sep="\t")
# load
# foo = read.table(paste0(sub("dumpdata","trueeffect_trial",inputfile), trial),
#                  header=T, sep="\t");
# data[, names(foo)[5:6]] = foo[, 5:6];


# # transform to t using bivnormtobivt
# sigma1diag = h1           /M*(data$N1*data$score1 + M);
# sigma2diag = h2           /M*(data$N2*data$score2 + M);
# sigmaXdiag = sqrt(h1*h2)*x/M*sqrt(data$N1*data$N2)*(data$scoreX);
# r = sigmaXdiag/sqrt(sigma1diag*sigma2diag);
# summary(r);
# rmax=1; r[r > rmax]=rmax;   # SET MANUALLY if bivnormtobivt makes error below
# source("~/Documents/R/bivt/bivt.R", chdir=TRUE)
# clusterExport(cl, c("nu", "nu1nucombtonu2", "bivnormtobivt"))
# Zpartsimtrans =
#   #apply(
#   parApply(cl,   # parallel
#            cbind(
#              (sqrt(h1*(data$N1)/M)*data$Zpartsim1_0c + sqrt(h1)*data$rcg1_0) /
#                sqrt(sigma1diag),
#              (sqrt(h2*(data$N2)/M)*data$Zpartsim2_0c + sqrt(h2)*data$rcg2_0) /
#                sqrt(sigma2diag),
#              r ),
#            1,
#            function (y) {
#               ## TODO REWRITE TO bivnormchangecor bivnormtomvt bivnormchangecor
#               bivnormtobivt(y[1:2], y[3], nu, nu) } );
# data$Zsim1 = 
#   Zpartsimtrans[1,]*sqrt(sigma1diag) +
#   sqrt(1-h1)*data$rcg1_1;   # trial 0 already used above
# data$Zsim2 = 
#   Zpartsimtrans[2,]*sqrt(sigma2diag) +
#   sqrt(1-h2)*data$rcg2_1;

# transform to t using transtot12
# data$Zsim1 = 
#   (sqrt(h1*(data$N1)/M)*data$Zpartsim1_0c + sqrt(h1)*data$rcg1_0)*transtot12 +
#   sqrt(1-h1)*data$rcg1_1;   # trial 0 already used above
# data$Zsim2 = 
#   (sqrt(h2*(data$N2)/M)*data$Zpartsim2_0c + sqrt(h2)*data$rcg2_0)*transtot12 +
#   sqrt(1-h2)*data$rcg2_1;


## correlation scale
observed = (1/M)*(data$Zpartsim1)^2
## Z scale (covariance between SNPs is inaccurate)
observed = (1/M)*(data$Zpartsim1)^2*data$N1
observed = (1/M)*(data$Zpartsim1)^2*
  (data$N1*data$score1 + M)/data$score1
observed = (data$Zsim1)^2;

data[head(order((data$Zsim1)^2, decreasing=TRUE), 100),
     c("id","chr","pos","af1","score1","N1","Zpartsim1","observed")]

observed = sort(observed, decreasing=TRUE)
observedsample = observed[seq(1,length(observed),ceiling(length(observed)/np))]
observedsamplesim = observedsample
observed = (data$Z1)^2
observed = sort(observed, decreasing=TRUE)
observedsample = observed[seq(1,length(observed),ceiling(length(observed)/np))]
plot(observedsamplesim,observedsample)
lines(c(0,50),c(0,50),col="Red")


### Number of detectable loci
# common SNPs in two populations (eg. EAS_EUR)
result = c();
for (trial in 0:99) {
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
output = c(apply(result[,1:15], 2, mean), # no. of loci does not need scaling for M
           mean(result[,ncol(result)-2]), quantile(result[,ncol(result)-2], c(0.025,0.975)),
           mean(result[,ncol(result)-1]), quantile(result[,ncol(result)-1], c(0.025,0.975)),
           mean(result[,ncol(result)-0]), quantile(result[,ncol(result)-0], c(0.025,0.975)),
           max(data$N1), max(data$N2));
# One population dataset (eg. EAS_self)
otherMAF = "EURMAF";
otherMAF = "EASMAF";
result = c();
for (trial in 0:99) {
  x = data[, paste0("Zsim1_",trial)];
  x = x / pmax(abs(x),1);
  z = !is.na(data[,otherMAF]) & data[,otherMAF] >= 0.05;
  n = 1*10^5;
  p = gwaspowerncp(n, n*x^2/(1-x^2), 5*10^-8);
  y = (!is.na(p) & p >= 0.05);
  w100K = snpstowindow(data$chr[y], data$pos[y], p[y]);
  w100Ksharedvariant = snpstowindow(data$chr[y&z], data$pos[y&z], p[y&z]);
  n = 2*10^5;
  p = gwaspowerncp(n, n*x^2/(1-x^2), 5*10^-8);
  y = (!is.na(p) & p >= 0.05);
  w200K = snpstowindow(data$chr[y], data$pos[y], p[y]);
  w200Ksharedvariant = snpstowindow(data$chr[y&z], data$pos[y&z], p[y&z]);
  n = 5*10^5;
  p = gwaspowerncp(n, n*x^2/(1-x^2), 5*10^-8);
  y = (!is.na(p) & p >= 0.05);
  w500K = snpstowindow(data$chr[y], data$pos[y], p[y]);
  w500Ksharedvariant = snpstowindow(data$chr[y&z], data$pos[y&z], p[y&z]);
  n = max(data$N1);
  p = gwaspowerncp(n, n*x^2/(1-x^2), 5*10^-8);
  y = (!is.na(p) & p >= 0.05);
  wmaxN = snpstowindow(data$chr[y], data$pos[y], p[y]);
  wmaxNsharedvariant = snpstowindow(data$chr[y&z], data$pos[y&z], p[y&z]);
  result = c(result,
             sum(w100K$pow), sum(w200K$pow), sum(w500K$pow),
             sum(w100K$pow) -
               sum(w100Ksharedvariant$pow),
             sum(w200K$pow) -
               sum(w200Ksharedvariant$pow),
             sum(w500K$pow) -
               sum(w500Ksharedvariant$pow),
             sum(wmaxN$pow), sum(wmaxN$pow)-sum(wmaxNsharedvariant$pow))
}
result = matrix(result, ncol=8, byrow=TRUE);
output = c(apply(result[,1:6], 2, mean), # no. of loci does not need scaling for M
           mean(result[,ncol(result)-1]), quantile(result[,ncol(result)-1], c(0.025,0.975)),
           mean(result[,ncol(result)-0]), quantile(result[,ncol(result)-0], c(0.025,0.975)),
           max(data$N1));


n = 1*10^5;
n = 2*10^5;
n = 5*10^5;
b = 0.05 # 0.03 # 0.04
x = seq(-b, b, 0.001)
y = sapply(x, function(x){ gwaspowerncp(n, n*x^2/(1-x^2), 5*10^-8) })
ggplot(data.frame(x,y),
       aes(x=x, y=0, fill=y)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="brown", limits=c(0,1)) +
  labs(fill="Power") +
  scale_x_continuous(
    "Effect-size",
    limits=c(-b,b),
    position="top")
