#' EditVCF
#'
#' Edit the VCF (Variant Call Format) file, Creates file with SNPs data at the BAM directory called variantTable.csv
#' @param Directory The Path of to the BAM Directory
#' @param Organism "Human" or "Mouse"
#' @export
#' @return None

EditVCF<-function(Directory,Organism){
  print("Editing VCF File")
  Directory  
  Dir=Directory
  file = "variants_output.vcf"
  path = paste(Dir, file, sep="")
  readData = read.delim(path,as.is=T)
  readData=as.character(readData[-c(1:which(readData=="#CHROM")-1),1])
  
  jump = 10
  startChr = 1+jump
  startPos = 2+jump
  startInfo = 10+jump
  
  len = length(readData)
  chrRegex ="^chr(\\w+)$"
  infoRegex = "^([01])\\/([01]):(\\d+)\\,(\\d+):(\\d+):\\d+:\\d+\\,\\d+\\,\\d+$"
  
  chrVector =  readData[startChr]
  posVector =  readData[startPos]
  infoVector = readData[startInfo]
  
  while (startInfo + jump < len) {
    startChr = startChr + jump
    startPos = startPos + jump
    startInfo = startInfo + jump
    chrVector = append(chrVector, readData[startChr])
    posVector = append(posVector, readData[startPos])
    infoVector = append(infoVector, readData[startInfo])
  }
  
  chrNum = gsub(chrRegex, "\\1", chrVector)
  if (Organism=="Human"){
    chrNum[chrNum=="X"]="23"
    chrNum[chrNum=="Y"]="24"}
  
  if (Organism=="Mouse"){
    chrNum[chrNum=="X"]="20"
    chrNum[chrNum=="Y"]="21"}
  
  chrNum =  as.numeric(chrNum)
  Karyotape = 10*abs(as.numeric(gsub(infoRegex, "\\1", infoVector))-as.numeric(gsub(infoRegex, "\\2", infoVector)))
  AD1 = as.numeric(gsub(infoRegex, "\\3", infoVector))
  AD2 = as.numeric(gsub(infoRegex, "\\4", infoVector))
  DP = as.numeric(gsub(infoRegex, "\\5", infoVector))
  
  posVector = as.numeric(posVector)
  
  table = data.frame("chr" = chrNum, "position" = posVector, "AD1" = AD1, "AD2" = AD2, "DP" = DP, "Karyotape" = Karyotape)
  
  fileName = "variantTable.csv"
  pathToSave = paste(Dir, fileName, sep="")
  table[is.na(table)] = 0
  write.table(table,pathToSave , sep="\t",row.names=F,quote=F)
  return(table)
}