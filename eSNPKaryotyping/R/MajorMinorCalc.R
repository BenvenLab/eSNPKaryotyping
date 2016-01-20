#' MajorMinorCalc
#'
#' Major to Minor Calculation function
#' @param Table The variable containing the table of SNPs
#' @param minDP minimal reading depth required from accepted SNP, usually set to 20
#' @param maxDP maximal reading depth required from accepted SNP, usually set to very high number
#' @param minAF minimal minor allele frequency, usually set to 0.2-0.25
#' @export
#' @return None

MajorMinorCalc<-function(Table,minDP,maxDP,minAF){
  
  Table[is.na(Table)] = 0
  newTable = Table[Table$DP >= minDP,]
  newTable = newTable[newTable$DP <= maxDP,]
  AF1 = newTable$AD1/newTable$DP
  AF2 = newTable$AD2/newTable$DP
  newTable = data.frame("chr" = newTable$chr, "position" = newTable$position, "AF1" = AF1, "AF2" = AF2)
  frequncyTable = newTable[newTable$AF1 >= minAF,]
  frequncyTable = frequncyTable[frequncyTable$AF2 >= minAF,]
  orderedTable = Sort_major_minor(data=frequncyTable, col1=3, col2=4)
  MajorMinor = orderedTable$AF1/orderedTable$AF2
  orderedTable["MajorMinor"] = MajorMinor
  return(orderedTable)
}