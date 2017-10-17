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
#' @param   dbSNP_Data_Directory The path for the directory where the edited dbSNP file are (created previously by the Edit_dbSNP_Files function)
#' @param   dbSNP_File_Name The edited dbSNP file names, without the chromosome number
#' @param   Genome_Fa_dict the path for the dictionary file created for the whole genome FASTA file.
#' @param   Organism "Human" or "Mouse"
#' @export
#' @return None

SRAtoSNP<-function(File,Library_Type,TopHat_Threads,Transcript_Anotations,Bowtie_Genome_Index,Genome_FA,SRAPath,Picard_Path,GATK_Path,dbSNP_Data_Directory,dbSNP_File_Name,Genome_Fa_dict,Organism){
  OpenSRA(File,Library_Type,SRAPath)
  Output_Dir=unlist(strsplit(File,".sra",fixed=T))
  
  dir=paste(Output_Dir,"/FASTQ/",sep="")
  Tophat(Directory=dir,Library_Type=Library_Type,Threads = TopHat_Threads,Transcripts_Annotations = Transcript_Anotations,Bowtie_Genome_Index = Bowtie_Genome_Index)
  dir=paste(Output_Dir,"/BAM/",sep="")
  CreateVCF(Directory=dir,Genome_Fa = Genome_FA,Picard_Path = Picard_Path,GATK_Path = GATK_Path)
  table = EditVCF(Directory=dir,Organism=Organism)
  table$chr=as.numeric(table$chr)
  table=table[order(table$chr,table$position),]
  table=table[table$chr>0,]
  table2=MajorMinorCalc(Table =table ,minDP =20 ,maxDP = 10000,minAF =0.2 )
  tbl=DeletionTable(Directory =dir ,Table = table2 ,dbSNP_Data_Directory = dbSNP_Data_Directory,dbSNP_File_Name = dbSNP_File_Name,Genome_Fa_dict = Genome_Fa_dict,Organism = Organism)
  
}
