#' Tophat
#'
#' Align reads to genome using Tophat2. The Script creates new directory called BAM in the same path were the FASTQ dir is located.
#' @param Directory The Path of the directory were the FASTQ files are
#' @param Library_Type "Paired" or "Single"
#' @param Threads Number of threads used for the alignment
#' @param Transcripts_Annotations path for the transcript annotation GTF file.
#' @param Bowtie_Genome_Index The basename of the genome index to be searched. The basename is the name of any of the index files up to but not including the first period.
#' @export
#' @return None


Tophat<-function(Directory,Library_Type,Threads,Transcripts_Annotations,Bowtie_Genome_Index){
  Output_Dir=strsplit(Directory,"/F",fixed=T)
  Output_Dir=unlist(Output_Dir)
  Output_Dir=Output_Dir[1]
  Output_Dir=paste(Output_Dir,"/BAM",sep="")
  dir.create(Output_Dir)
  command = "Wrong Library Type"  
  
  if(Library_Type=="Single"){
    setwd(Directory)
    files=system("ls",intern=TRUE)
    setwd(Output_Dir)
    read=paste(Directory, files, sep="")
    command=paste("tophat2 -o ./ --no-coverage-search --max-multihits 1 -G ",Transcripts_Annotations," -p ",Threads," ",Bowtie_Genome_Index ," ", read, sep="")
  }
  
  if(Library_Type=="Paired"){
    setwd(Directory)
    first_read= system("ls *_1.*",intern=TRUE)
    second_read = system("ls *_2.*",intern=TRUE)
    setwd(Output_Dir)
    read1=paste(Directory, first_read, sep="")
    read2=paste(Directory, second_read, sep="")
    command = paste("tophat2 -o ./ --no-coverage-search --max-multihits 1 -G ",Transcripts_Annotations," -p ",Threads," ",Bowtie_Genome_Index, " ", read1," ", read2, sep="")
  }
  
  
  
  print(command)
  if ( command != "Wrong Library Type"  ){
    system(command)
    system ("samtools flagstat accepted_hits.bam > flagstat.txt")
  }
  
  
}
