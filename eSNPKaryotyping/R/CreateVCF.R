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
  newDir = substring(Directory,3)
  Genome_Fa=substring(Genome_Fa,3)
  
  setwd(paste( "~/",newDir, sep=""))
  comm = paste("java -jar ",Picard_Path, "picard.jar AddOrReplaceReadGroups.jar  I=accepted_hits.bam O=accepted_hits_rg.bam ID=\"n\" LB=\"lb\" PL=\"illumina\" PU=\"pu\" SM=\"ES\"",sep="")
  print("======================")
  print(comm)
  print("======================")
  system(comm)
  setwd("~/")
  
  comm = paste("java -jar ",Picard_Path, "picard.jar ReorderSam.jar I=",newDir, "accepted_hits_rg.bam O=",newDir, "accepted_hits_rg_sorted.bam R=",Genome_Fa, sep="")
  print("======================")
  print(comm)
  print("======================")
  system(comm)
  
  comm = paste("java -jar ",Picard_Path, "picard.jar BuildBamIndex.jar I=",newDir, "accepted_hits_rg_sorted.bam", sep="")
  print("======================")
  print(comm)
  print("======================")
  system(comm)
  
  comm = paste("java -jar ",GATK_Path, "GenomeAnalysisTK.jar -T SplitNCigarReads -R ",Genome_Fa," -I ",newDir, "accepted_hits_rg_sorted.bam -o ",newDir, "split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS", sep="") 
  altComm = paste(comm,"--fix_misencoded_quality_scores")
  print("======================")
  print(comm)
  print("======================")
  system(comm)
  if (!file.exists( paste(newDir, "split.bam",sep=""))) {
    print("======================")
    print("Alternative Command")
    print(altComm)
    print("======================")
    system(altComm)
    
  }
  
  
  comm = paste("java -jar ",GATK_Path, "GenomeAnalysisTK.jar -T HaplotypeCaller -R ",Genome_Fa," -I ",newDir, "split.bam -recoverDanglingHeads -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o ",newDir, "variants_output.vcf", sep="")
  print("======================")
  print(comm)
  print("======================")
  system(comm)
}
