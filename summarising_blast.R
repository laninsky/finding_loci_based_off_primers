library(readr)
library(dplyr)

primertable <- read_table2((as.matrix(read_table2("params.txt",col_names=FALSE)[2,1])),col_names=FALSE)

blastresults <- list.files(pattern=".blast")

for (i in 1:length(blastresults)) {
  tempblastresults <- read_table2(blastresults[i],col_names=FALSE)
  tempblastresults %>% group_by
