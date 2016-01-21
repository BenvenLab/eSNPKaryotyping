

**********************
# eSNPKaryotyping
Analysis of Chromosomal Aberrations from RNA-Seq sequencing data.


  Weissbein et al, 2015.  


  Before Use  

# The following programs should be installed on your computer:
# 1. Tophat2 (install to the Linux path. If not, give the full path to the tophat2 function in the Tophat function)
# 2. Bowtie2
# 3. GATK
# 4. Picard
# 5. SAMtools
# 6. SRAtools (for extracting SRA files into FASTQ)


# The following Files Should be on your computer:
# 1. Bowtie2 Genome Index (Available for download at Illumina iGenome Website or build using Bowtie2)
# 2. Transcription annotation file (GTF) (Available at Illumina iGenome Website, USCS or any other source)
# 3. Whole Genome FASTA (Available at Illumina iGenome Website or any other source)
# 4. Picard genome dictionary file, should be at the same directory where the genome FASTA file is, created using Picrad "CreateSequenceDictionary" function.
# 5. directory with dbSNP build 142 common SNP files (GTF Format). A file for each chromosome (Available from UCSC Table browser tool, called snp142common table). Each file name should end with the chromosome number, X and Y should get a name by their number (For example: in human chromosome X is 23 and Y is 24)

# Install and load the following R packages :
install.packages("zoo")
library("zoo")
install.packages("gplots")
library("gplots")

# Before LOH analysis, the dbSNP files need to be edited using the Edit_dbSNP_Files function (done only once):
Edit_dbSNP_Files(Directory = , File_Name = , Organism = )
# Argument: 1. Directory - the directory were the files are, one GTF file per chromosome
#           2. File_Name - the files name, without the number of the chromosomes
#           3. Organism - "Human" or "Mouse"



**********************
***   Genomic integrity analysis workflow   ***
**********************
 
# 1. Open SRA file into FASTQ file. The script will create new directory called FASTQ with the fastq files:
OpenSRA(File = ,Library_Type = ,SRAPath = )
# Argument: 1. File - The Path to the SRA File
#           2. Library_Type - "Paired" or "Single"
#           3. SRAPath - path for the SRAtools BIN directory (example: SRAPath = "~/SRAtoolkit/sratoolkit.2.5.4-ubuntu64/bin/")


# 2. Align reads to genome using Tophat2. The Script creates new directory called BAM in the same path were the FASTQ dir is located
Tophat(Directory = ,Library_Type = ,Threads = ,Transcripts_Annotations = ,Bowtie_Genome_Index = )
# Argument: 1. Directory - The Path of the directory were the FASTQ files are
#           2. Library_Type - "Paired" or "Single"
#           3. Threads - Number of threads used for the alignment
#           4. Transcripts_Annotations - path for the transcript annotation GTF file.
#           5. Bowtie_Genome_Index - The basename of the genome index to be searched. The basename is the name of any of the index files up to but not including the first period.


# 3. Create the VCF (Variant Call Format) file 
CreatVCF(Directory = ,Genome_Fa = ,Picard_Path = ,GATK_Path = )
# Argument: 1. Directory - The Path of to the BAM Directory
#           2. Genome_Fa - Path for whole genome FASTQ file
#           3. Picard_Path - Path to the Picard directory, with all the JAR files
#           4. GATK_Path - Path to the GATK directory, where the GenomeAnalysisTK.jar file is located


# 4. Edit the VCF (Variant Call Format) file, Creates file with SNPs data at the BAM directory called variantTable.csv
table=EditVCF(Directory = ,Organism = )
# Argument: 1. Directory - The Path of to the BAM Directory
#           2. Organism - "Human" or "Mouse"


# To run the sample again, just read the variantTable.csv file into a variable using the read.delim function


# 5. Run the following lines:
table$chr=as.numeric(table$chr)
table=table[order(table$chr,table$position),]
table=table[table$chr>0,]


# 6. Run the MajorMinorCalc function
table2=MajorMinorCalc(Table = ,minDP = ,maxDP = ,minAF = )
# Argument: 1. Table - The variable containing the table of SNPs
#           2. minDP - minimal reading depth required from accepted SNP, usually set to 20
#           3. maxDP - maximal reading depth required from accepted SNP, usually set to very high number
#           4. minAF - minimal minor allele frequency, usually set to 0.2-0.25


# 7. Plot Allelic ratio along the genome for duplication detection:
PlotGenome(orderedTable = ,Window = ,Ylim = ,PValue = ,Organism = )
# Argument: 1. orderedTable - The variable containing the output of the MajorMinorCalc function
#           2. Window - an odd number determining the window of the moving average plot, usually 151
#           3. Ylim - the plot y axes maximal limit, usually 3
#           4. PValue - Determine if to add the P-Value bar to the graph, accepts "TRUE" or "FALSE" values
#           5. Organism - "Human" or "Mouse"


# 8. intersect the observed SNPs with the common SNPs table from dbSNPs, Creates file with the LOH data called Deletions.txt
tbl=DeletionTable(Directory = ,Table = ,dbSNP_Data_Directory = ,dbSNP_File_Name = ,Genome_Fa_dict = ,Organism = )
# Argument: 1. Directory - The Path of to the BAM Directory, containing the variantTable.csv file
#           2. Table - The variable containing the output of the MajorMinorCalc function
#           3. dbSNP_Data_Directory - The path for the directory where the edited dbSNP file are (created previously by the Edit_dbSNP_Files function)
#           4. dbSNP_File_Name - The edited dbSNP file names, without the chromosome number
#           5. Genome_Fa_dict - the path for the dictionary file created for the whole genome FASTA file.
#           6. Organism - "Human" or "Mouse"


# To run the LOH analysis again, just read the Deletions.txt file into a variable using the read.delim function


# 9. Plot each SNP, without any summarization
Plot_Zygosity_Sinle(Table = ,Organism = )
# Argument: 1. Table - The LOH table containing the output of the DeletionTable function
#           2. Organism - "Human" or "Mouse"


# 10. Plot blocks of heterozygous and homozygous SNPs
Plot_Zygosity_Blocks(Table = ,Window = ,Max = ,Max2 = ,Organism = )
# Argument: 1. Table - The deletion table containing the output of the DeletionTable function
#           2. window - the block size in bp, usually 1500000
#           3. Max - How many Heterozygouse SNP need to be in a block to get the full color, usually 6
#           4. Max2 - How many Homozygouse SNP need to be in a block to get the full color, usually 60
#           5. Organism - "Human" or "Mouse"




# In case of analyzing SRA file, it is also possible to use the following function that performs step 1-4 in one function
SRAtoSNP(File = ,Library_Type = ,TopHat_Threads = ,Transcript_Anotations = ,Bowtie_Genome_Index = ,Genome_FA = ,Picard_Path = ,GATK_Path = ,Organism = )
# Argument: 1. File - The Path to the SRA File
#           2. Library_Type - "Paired" or "Single"
#           3. TopHat_Threads - Number of threads used for the alignment
#           4. Transcripts_Annotations - path for the transcript annotation GTF file.
#           5. Bowtie_Genome_Index - The basename of the genome index to be searched. The basename is the name of any of the index files up to but not including the first period.
#           6. Genome_Fa - Path for whole genome FASTA file
#           7. Picard_Path - Path to the Picard directory, with all the JAR files
#           8. GATK_Path - Path to the GATK directory, with GenomeAnalysisTK.jar file
#           9. Organism - "Human" or "Mouse"
#  

