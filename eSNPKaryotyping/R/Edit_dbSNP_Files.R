#' Edit_dbSNP_Files
#'
#' Before LOH analysis, the dbSNP files need to be edited using the Edit_dbSNP_Files function (done only once):
#' @param Directory the directory were the files are, one GTF file per chromosome
#' @param File_Name the files name, without the number of the chromosomes
#' @param Organism "Human" or "Mouse"
#' @export
#' @return None

Edit_dbSNP_Files<- function(Directory,File_Name,Organism){
  
  if(Organism=="Human"){
    for (i in 1:24){
      print(paste("Chromosme:",i))
      p = paste(Directory,File_Name,i,sep="")
      readData = read.delim(p, as.is = T,header = F,)
      chrRegex ="^chr(\\w+)$"
      chrNum = gsub(chrRegex, "\\1", readData$V1)
      chrNum[chrNum=="X"] = "23"
      chrNum[chrNum=="Y"] = "24"
      snpTable = data.frame(chr = chrNum, start = readData$V4, end = readData$V5, snp = array(0,length(readData$V1)))
      snpTable = snpTable[order(as.numeric(snpTable[,1]), as.numeric(snpTable[,2])),]
      d=duplicated(snpTable,MARGIN = 1)
      snpTable=snpTable[d==FALSE,]
      size=snpTable$end-snpTable$start
      snpTable=snpTable[size==0,]
      output=paste(Directory,"Edited_Common_chr",i,sep="")
      write.table(snpTable, output, sep="\t",row.names=F,quote=F)
    }}
  
  if(Organism=="Mouse"){
    for (i in 1:21){
      print(paste("Chromosme:",i))
      p = paste(Directory,File_Name,i,sep="")
      readData = read.delim(p, as.is = T,header = F,)
      chrRegex ="^chr(\\w+)$"
      chrNum = gsub(chrRegex, "\\1", readData$V1)
      chrNum[chrNum=="X"] = "20"
      chrNum[chrNum=="Y"] = "21"
      snpTable = data.frame(chr = chrNum, start = readData$V4, end = readData$V5, snp = array(0,length(readData$V1)))
      snpTable = snpTable[order(as.numeric(snpTable[,1]), as.numeric(snpTable[,2])),]
      d=duplicated(snpTable,MARGIN = 1)
      snpTable=snpTable[d==FALSE,]
      size=snpTable$end-snpTable$start
      snpTable=snpTable[size==0,]
      output=paste(Directory,"Edited_Common_chr",i,sep="")
      write.table(snpTable, output, sep="\t",row.names=F,quote=F)
    }}
  
}