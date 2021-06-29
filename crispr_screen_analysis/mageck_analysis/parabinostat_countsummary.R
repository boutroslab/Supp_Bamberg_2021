Sweave("parabinostat_countsummary.Rnw");
library(tools);

texi2dvi("parabinostat_countsummary.tex",pdf=TRUE);

