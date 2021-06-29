Sweave("vorinostat_countsummary.Rnw");
library(tools);

texi2dvi("vorinostat_countsummary.tex",pdf=TRUE);

