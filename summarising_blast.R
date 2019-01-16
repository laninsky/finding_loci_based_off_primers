library(readr)

primertable <- read_table2((as.matrix(read_table2("params.txt",col_names=FALSE)[2,1])),col_names=FALSE)

blastresults <- list.files(pattern=".blast")
