#' SRAtoSNP
#'
#' In case of analyzing SRA file, it is also possible to use the following function that performs step 1-4 in one function
#' @param   File The Path to the SRA File
#' @param   Library_Type "Paired" or "Single"
#' @param   TopHat_Threads Number of threads used for the alignment
#' @param   Transcripts_Annotations path for the transcript annotation GTF file.
#' @param   Bowtie_Genome_Index The basename of the genome index to be searched. The basename is the name of any of the index files up to but not including the first period.
#' @param   Genome_Fa Path for whole genome FASTA file
#' @param   Picard_Path Path to the Picard directory, with all the JAR files
#' @param   GATK_Path Path to the GATK directory, with GenomeAnalysisTK.jar file
#' @param   Organism "Human" or "Mouse"
#' @export
#' @return None

SRAtoSNP<-function(File,Library_Type,TopHat_Threads,Transcript_Anotations,Bowtie_Genome_Index,Genome_FA,Picard_Path,GATK_Path,Organism){
  OpenSRA(File,Library_Type)
  Output_Dir=unlist(strsplit(File,".sra",fixed=T))
  
  dir=paste(Output_Dir,"/FASTQ/",sep="")
  Tophat(Directory=dir,Library_Type=Library_Type,Threads = TopHat_Threads,Transcripts_Annotations = Transcript_Anotations,Bowtie_Genome_Index = Bowtie_Genome_Index)
  dir=paste(Output_Dir,"/BAM/",sep="")
  CreatVCF(Directory=dir,Genome_Fa = Genome_FA,Picard_Path = Picard_Path,GATK_Path = GATK_Path)
  table = EditVCF(Directory=dir,Organism=Organism)
  return(table)
}
