#' DeletionTable
#'
#' Intersect the observed SNPs with the common SNPs table from dbSNPs, Creates file with the LOH data called Deletions.txt
#' @param Directory The Path of to the BAM Directory, containing the variantTable.csv file
#' @param Table The variable containing the output of the MajorMinorCalc function
#' @param dbSNP_Data_Directory The path for the directory where the edited dbSNP file are (created previously by the Edit_dbSNP_Files function)
#' @param dbSNP_File_Name The edited dbSNP file names, without the chromosome number
#' @param Genome_Fa_dict the path for the dictionary file created for the whole genome FASTA file.
#' @param Organism "Human" or "Mouse"
#' @export
#' @return None

DeletionTable<-function(Directory,Table,dbSNP_Data_Directory,dbSNP_File_Name,Genome_Fa_dict,Organism){
  print("Reading SNPs table")
  i=1
  
  if(Organism=="Human"){mx=25}
  if(Organism=="Mouse"){mx=22}
  
  while (i <mx){
    print(paste("Chromosome:",i))
    chrTable = read.delim(paste(dbSNP_Data_Directory,dbSNP_File_Name,i,sep=""))
    if (i==1) {snpTable = chrTable}
    if (i>1) {snpTable = rbind(snpTable,chrTable)}
    i=i+1
  }
  
  
  table2=Table
  colnames(table2)[2]="start"
  print("Merging Tables")
  x=merge(snpTable,table2,by = c("chr","start"),all.x=T)
  x=x[order(x$chr,x$start),]
  
  
  setwd(Directory)
  
  print("Indexing BAM File")
  system("samtools index accepted_hits.bam")
  
  dict=read.csv(Genome_Fa_dict,as.is=T)
  dict_type=grep(pattern = "chr",x = dict[1,1])
  if(length(dict_type)==0){dict_type=0}
  
  
  i=1
  while (i<mx){
    print(paste("Chromosome ",i, "| ",Sys.time(),sep=""))
    x1=x[x$chr==i,]
    
    
    if(dict_type==0){
      loc=paste(x1$chr[1],":",x1$start[1],"-",x1$start[dim(x1)[1]],sep="")
      if(Organism=="Human"){
        if (i==23) {loc=paste("X:",x1$start[1],"-",x1$start[dim(x1)[1]],sep="")}
        if (i==24) {loc=paste("Y:",x1$start[1],"-",x1$start[dim(x1)[1]],sep="")}}
      if(Organism=="Mouse"){
        if (i==20) {loc=paste("X:",x1$start[1],"-",x1$start[dim(x1)[1]],sep="")}
        if (i==21) {loc=paste("Y:",x1$start[1],"-",x1$start[dim(x1)[1]],sep="")}}
    }
    
    if(dict_type==1){
      loc=paste("chr",x1$chr[1],":",x1$start[1],"-",x1$start[dim(x1)[1]],sep="")
      if(Organism=="Human"){
        if (i==23) {loc=paste("chrX:",x1$start[1],"-",x1$start[dim(x1)[1]],sep="")}
        if (i==24) {loc=paste("chrY:",x1$start[1],"-",x1$start[dim(x1)[1]],sep="")}}
      if(Organism=="Mouse"){
        if (i==20) {loc=paste("chrX:",x1$start[1],"-",x1$start[dim(x1)[1]],sep="")}
        if (i==21) {loc=paste("chrY:",x1$start[1],"-",x1$start[dim(x1)[1]],sep="")}}
    }
    
    
    
    command=paste("samtools depth -r ",loc ," accepted_hits.bam > reads-per-position.txt",sep="")
    system(command)
    system ("awk -F \" \" '($3 >20){print $0}' reads-per-position.txt >reads-per-position2.txt")
    chr=read.delim("reads-per-position2.txt",header = F)
    colnames(chr)=c("chr","start","Depth")
    x2=merge(x1,chr,by="start")
    x2=cbind(x2,Depth_group=x2$Depth)
    x2$Depth_group[x2$Depth<50] ="20-50"
    x2$Depth_group[x2$Depth>49 & x2$Depth<100] ="50-100"
    x2$Depth_group[x2$Depth>99 & x2$Depth<200] ="100-200"
    x2$Depth_group[x2$Depth>199 & x2$Depth<500] ="200-500"
    x2$Depth_group[x2$Depth>499] =">500"
    x2$Depth_group=as.factor(x2$Depth_group)
    if (i==1) {tbl=x2}
    if (i>1) {tbl=rbind(tbl,x2)}
    i=i+1
  }
  
  
  print("Writing Table")
  tbl[is.na(tbl)]=0
  write.table(tbl,"Deletions.txt" , sep="\t",row.names=F,quote=F)
  return(tbl)
}