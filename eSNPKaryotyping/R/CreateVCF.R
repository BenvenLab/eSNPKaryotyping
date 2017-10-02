#' CreateVCF
#'
#' Create the VCF (Variant Call Format) file 
#' @param Directory The Path of to the BAM Directory
#' @param Genome_Fa Path for whole genome FASTQ file
#' @param Picard_Path Path to the Picard directory, with all the JAR files
#' @param GATK_Path Path to the GATK directory, where the GenomeAnalysisTK.jar file is located
#' @export
#' @return None


CreateVCF<-function(Directory,Genome_Fa,Picard_Path,GATK_Path){

  setwd(newDir)
  comm = paste("java -jar ",Picard_Path, "picard.jar AddOrReplaceReadGroups  I=",Directory,"accepted_hits.bam O=accepted_hits_rg.bam ID=\"n\" LB=\"lb\" PL=\"illumina\" PU=\"pu\" SM=\"ES\"",sep="")
  print("======================")
  print(comm)
  print("======================")
  system(comm)
  setwd("~/")
  
  comm = paste("java -jar ",Picard_Path, "picard.jar ReorderSam I=",Directory, "accepted_hits_rg.bam O=",newDir, "accepted_hits_rg_sorted.bam R=",Genome_Fa, sep="")
  print("======================")
  print(comm)
  print("======================")
  system(comm)
  
  comm = paste("java -jar ",Picard_Path, "picard.jar BuildBamIndex I=",Directory, "accepted_hits_rg_sorted.bam", sep="")
  print("======================")
  print(comm)
  print("======================")
  system(comm)
  
  comm = paste("java -jar ",GATK_Path, "GenomeAnalysisTK.jar -T SplitNCigarReads -R ",Genome_Fa," -I ",Directory, "accepted_hits_rg_sorted.bam -o ",Directory, "split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS", sep="") 
  altComm = paste(comm,"--fix_misencoded_quality_scores")
  print("======================")
  print(comm)
  print("======================")
  system(comm)
  if (!file.exists( paste(Directory, "split.bam",sep=""))) {
    print("======================")
    print("Alternative Command")
    print(altComm)
    print("======================")
    system(altComm)
    
  }
  
  
  comm = paste("java -jar ",GATK_Path, "GenomeAnalysisTK.jar -T HaplotypeCaller -R ",Genome_Fa," -I ",Directory, "split.bam -dontUseSoftClippedBases -stand_call_conf 10.0 -o ",Directory, "variants_output.vcf", sep="")
  print("======================")
  print(comm)
  print("======================")
  system(comm)
}
