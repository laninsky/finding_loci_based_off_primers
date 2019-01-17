library(readr)
library(dplyr)

primertable <- read_table2((as.matrix(read_table2("params.txt",col_names=FALSE)[2,1])),col_names=FALSE)
length_slush <- as.numeric(read_table2("params.txt",col_names=FALSE)[3,1])

blastresults <- list.files(pattern="*.blast$")

output <- as_tibble()

for (i in 1:length(blastresults)) {
  tempblastresults <- read_table2(paste(i,".blast",sep=""),col_names=FALSE)
  
  summarised_matches <- tempblastresults %>% 
    group_by(X2) %>% 
    summarise(nomatches=n(),
              max_evalue=max(X6),
              min_bitscore=min(X5),
              startpos=min(X3,X4),
              endpos=max(X3,X4),
              length=(max(X3,X4)-min(X3,X4)+1))
  
  winning_evalue <- min(summarised_matches$max_evalue)
  winning_bitscore <- max(summarised_matches$min_bitscore) 
              

  summarised_matches <- summarised_matches %>% 
    rowwise %>% 
    mutate(matchrank=ifelse(nomatches>=2,1,0) +
             ifelse(max_evalue==winning_evalue,1,0) + 
             ifelse(min_bitscore==winning_bitscore,1,0) +
             ifelse((length > (primertable[i,4]-length_slush)),
                    ifelse((length < (primertable[i,4]+length_slush)),1,0),0))

  winning_rank <- max(summarised_matches$matchrank)
  
  tempoutput <- cbind(primertable[i,1],filter(summarised_matches,matchrank==winning_rank))
  names(tempoutput)[1:2] <- c("locus_name","scaffold_name")
  output <- rbind(output,tempoutput)
}

write.table(output,"best_blast_matches.txt",quote=FALSE,row.names=FALSE)
