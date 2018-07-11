install.packages("MCMC.OTU")
library("MCMC.OTU")
# Symbiodinium sp diversity in two coral species at two reefs (banks)
data(green.data)
# removing outliers
goods=purgeOutliers(green.data,
  count.columns=c(4:length(green.data[1,])))
# stacking the data table
gs=otuStack(
  goods,
  count.columns=c(4:length(goods[1,])),
  condition.columns=c(1:3)
)
# fitting the model
mm=mcmc.otu(
  fixed="bank+species+bank:species",
  data=gs
)
# selecting the OTUs that were modeled reliably
acpass=otuByAutocorr(mm,gs)
# calculating effect sizes and p-values:
ss=OTUsummary(mm,gs,summ.plot=FALSE)
# correcting for mutliple comparisons (FDR)
ss=padjustOTU(ss)
# getting significatly changing OTUs (FDR<0.05)
sigs=signifOTU(ss)
# plotting them
ss2=OTUsummary(mm,gs,otus=sigs)
# bar-whiskers graph of relative changes:
# ssr=OTUsummary(mm,gs,otus=signifOTU(ss),relative=TRUE)
# displaying effect sizes and p-values for significant OTUs
ss$otuWise[sigs]
