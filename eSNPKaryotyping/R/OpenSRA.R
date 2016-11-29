#' OpenSRA
#'
#' Open SRA file into FASTQ file. The script will create new directory called FASTQ with the fastq files.
#' @param File The Path to the SRA File
#' @param Library_Type "Paired" or "Single"
#' @param SRAPath path for the SRAtools BIN directory (example: SRAPath = "~/SRAtoolkit/sratoolkit.2.5.4-ubuntu64/bin/")
#' @export
#' @return None

OpenSRA<-function(File,Library_Type,SRAPath){
  Output_Dir=unlist(strsplit(File,".sra",fixed=T))
  dir.create(Output_Dir)
  Output_Dir=paste(Output_Dir,"/FASTQ/",sep="")
  dir.create(Output_Dir)
  command = "Wrong Library Type"  
  
  if(Library_Type=="Single"){
    command = paste(SRAPath,"fastq-dump ", File, " -O ", Output_Dir,sep="")}
  
  if(Library_Type=="Paired"){
    command = paste(SRAPath,"fastq-dump --split-files ", File, " -O ", Output_Dir,sep="")}
  
  print(command)
  
  if ( command != "Wrong Library Type"  ){system(command)}
}
