#This script will take a set of multiple alignment
#coronavirus genome files and a list of proteins and generate a table
#That divides the genome alignment files into proteins.
#In its current iteration, this script will divide the coronavirus
#genomes into 25 different protein sequences
library(ggplot2)
library(stringr)
library(Biostrings)
nameList=c()
index=1
setwd("C:/Users/tyson/OneDrive/Desktop/Coronavirus Proteins")
#Your alignment files go here; the naming scheme should be the same
filepaths=c("COVID_genome_translation_rf1 alignment.fasta",
            "COVID_genome_translation_rf2 alignment.fasta",
            "COVID_genome_translation_rf3 alignment.fasta")
#Create objects for each genome sequence across the the three alignments
for(x in filepaths){
  splitOne=str_split(x," ")[[1]][1]
  thisRF=str_split(splitOne, "_")[[1]][[4]]
  print(thisRF)
  con = file(x, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    if(str_detect(line, "^>")){
      nameList=append(nameList,str_split(line, " ")[[1]][1])
    }else{
      assign(paste(thisRF, nameList[index], sep=""), line)
      index=index+1
    }
  }
  close(con)
}
tableExists=FALSE
listOfRows=c()
loopCompleted=FALSE
listOfFiles=list.files("C:/Users/tyson/OneDrive/Desktop/Coronavirus Proteins/Reference_Protein_Sequences", full.names = TRUE)
#Iterate through directory of amino acid sequences for 25 proteins
for(x in listOfFiles){
  if(str_detect(x, "fasta")){
    con = file(x, "r")
    while ( TRUE ) {
      line = readLines(con, n = 1)
      if ( length(line) == 0 ) {
        break
      }else if(str_detect(line, "^>")){
        name=line
      }else if (!str_detect(line, "^>") && !(length(line)==0)){
        thisSequence=line
      }
    }
    close(con)
    #See which rf gets the best alignment
    localAlign1 <-
      pairwiseAlignment(AAString(thisSequence),
                        AAString(`rf1>NC_045512_SARSCoV2_translation`),
                        gapOpening=-1, 
                        gapExtension=-0.1,
                        scoreOnly=FALSE,
                        type = "local")
    localAlign2 <-
      pairwiseAlignment(AAString(thisSequence),
                        AAString(`rf2>NC_045512_SARSCoV2_translation`),
                        gapOpening=-1, 
                        gapExtension=-0.1,
                        scoreOnly=FALSE,
                        type = "local")
    localAlign3 <-
      pairwiseAlignment(AAString(thisSequence),
                        AAString(`rf3>NC_045512_SARSCoV2_translation`),
                        gapOpening=-1, 
                        gapExtension=-0.1,
                        scoreOnly=FALSE,
                        type = "local")
    #The best alignment becomes our de facto alignment
    theMax=max(localAlign1@score,localAlign2@score, localAlign3@score)
    if(localAlign1@score==theMax){
      localAlign=localAlign1
    }
    if(localAlign2@score==theMax){
      localAlign=localAlign2
    }
    if(localAlign3@score==theMax){
      localAlign=localAlign3
    }
    listOfGenomes=c()
    for(q in ls()){
      nameOfProtein=str_remove_all(str_split(x, "/")[[1]][length(str_split(x, "/")[[1]])], ".fasta")
      if(str_detect(q, "rf1>") || str_detect(q, "rf2>") ||str_detect(q, "rf3>")){
        print(paste("Adding", paste(x, q, sep=" "), "to report", sep=" "))
        cat(paste(nameOfProtein, q, sep=" "), file="report.txt", append = TRUE)
        cat("\n", file="report.txt", append = TRUE)
        cat(substring(get(q),
                      localAlign@subject@range@start,
                      localAlign@subject@range@start + localAlign@subject@range@width-1), file = "report.txt", append = TRUE)
        cat("\n", file="report.txt", append = TRUE)
        listOfGenomes=append(listOfGenomes, substring(get(q),
                             localAlign@subject@range@start,
                             localAlign@subject@range@start + localAlign@subject@range@width-1))
        if(!loopCompleted){
          listOfRows=append(listOfRows, q)
        }
      }
    }
    loopCompleted=TRUE
    if(tableExists){
      myTable=as.data.frame(cbind(myTable,listOfGenomes))
      colnames(myTable)[length(colnames(myTable))]=nameOfProtein
    }else{
      tableExists=TRUE
      myTable=as.data.frame(listOfGenomes)
      colnames(myTable)[length(colnames(myTable))]=nameOfProtein
    }
  }
}
row.names(myTable)=listOfRows

finalTable=myTable
f=0
#Remove sequences with stop codons
for(w in 1:length(colnames(finalTable))){
  for (y in 1:length(row.names(finalTable))){
    if(str_detect(finalTable[y,w], "\\*")){
      finalTable[y,w]=" "
    }else{
      f=f+1
    }
  }
}
#Reorder the proteins so they are grouped together
myReorder=c()
for(i in 1:18){
  myReorder=append(myReorder, c(i, i+18, i+36))
}

write.table(finalTable[myReorder,], file="SequenceTable.csv", quote=FALSE, sep=",")
