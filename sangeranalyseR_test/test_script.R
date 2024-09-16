

library(sangeranalyseR)

# There is only four F AB1 files in dataDir
# To test if the program can deal with all F AB1 files
# You may have to change dataDir to your own path where you put the test data
dataDir <- 'D:/work/content/sq/article/sanger/R/gittest/t'

my_contigs <- SangerAlignment(ABIF_Directory     = dataDir,
                              REGEX_SuffixForward = "_[0-9]*_F.ab1$",
                              REGEX_SuffixReverse = "_[0-9]*_R.ab1$")
#launchApp(my_contigs)
writeFasta(my_contigs)
generateReport(my_contigs)



#**********************         console output:           **********************#



# INFO [2024-09-13 19:18:43] #################################################
# INFO [2024-09-13 19:18:43] #### Start creating SangerAlignment instance ####
# INFO [2024-09-13 19:18:43] #################################################
# INFO [2024-09-13 19:18:43]   >> You are using Regular Expression Method to group AB1 files!
# INFO [2024-09-13 19:18:43] ========================================================
# INFO [2024-09-13 19:18:43] ================ Creating 'SangerContig' ===============
# INFO [2024-09-13 19:18:43] ========================================================
# INFO [2024-09-13 19:18:43]   >> Contig Name: 'Achl_ACHLO006-09'
# SUCCESS [2024-09-13 19:18:43] --------------------------------------------------------
# SUCCESS [2024-09-13 19:18:43] -------- 'SangerRead' S4 instance is created !! --------
# SUCCESS [2024-09-13 19:18:43] --------------------------------------------------------
# SUCCESS [2024-09-13 19:18:43]    >> 'Achl_ACHLO006-09_1_F.ab1' is created (Forward Read; ABIF).
# INFO [2024-09-13 19:18:43]    >> The number of reads detected: 1
# WARN [2024-09-13 19:18:43] REGEX_MATCH_WARN
# Your 'contigName' and 'REGEX_SuffixReverse' regular expression parameters can not match any reverse reads.
# ERROR [2024-09-13 19:18:43] READ_NUMBER_ERROR
# The number of your total reads is 1.
# Number of total reads has to be equal or more than 2 ('minReadsNum' that you set)
# INFO [2024-09-13 19:18:43] ========================================================
# INFO [2024-09-13 19:18:43] ================ Creating 'SangerContig' ===============
# INFO [2024-09-13 19:18:43] ========================================================
# INFO [2024-09-13 19:18:43]   >> Contig Name: 'Achl_ACHLO007-09'
# SUCCESS [2024-09-13 19:18:43] --------------------------------------------------------
# SUCCESS [2024-09-13 19:18:43] -------- 'SangerRead' S4 instance is created !! --------
# SUCCESS [2024-09-13 19:18:43] --------------------------------------------------------
# SUCCESS [2024-09-13 19:18:43]    >> 'Achl_ACHLO007-09_1_F.ab1' is created (Forward Read; ABIF).
# INFO [2024-09-13 19:18:43]    >> The number of reads detected: 1
# WARN [2024-09-13 19:18:43] REGEX_MATCH_WARN
# Your 'contigName' and 'REGEX_SuffixReverse' regular expression parameters can not match any reverse reads.
# ERROR [2024-09-13 19:18:43] READ_NUMBER_ERROR
# The number of your total reads is 1.
# Number of total reads has to be equal or more than 2 ('minReadsNum' that you set)
# INFO [2024-09-13 19:18:43] ========================================================
# INFO [2024-09-13 19:18:43] ================ Creating 'SangerContig' ===============
# INFO [2024-09-13 19:18:43] ========================================================
# INFO [2024-09-13 19:18:43]   >> Contig Name: 'Achl_ACHLO040-09'
# SUCCESS [2024-09-13 19:18:43] --------------------------------------------------------
# SUCCESS [2024-09-13 19:18:43] -------- 'SangerRead' S4 instance is created !! --------
# SUCCESS [2024-09-13 19:18:43] --------------------------------------------------------
# SUCCESS [2024-09-13 19:18:43]    >> 'Achl_ACHLO040-09_1_F.ab1' is created (Forward Read; ABIF).
# INFO [2024-09-13 19:18:43]    >> The number of reads detected: 1
# WARN [2024-09-13 19:18:43] REGEX_MATCH_WARN
# Your 'contigName' and 'REGEX_SuffixReverse' regular expression parameters can not match any reverse reads.
# ERROR [2024-09-13 19:18:43] READ_NUMBER_ERROR
# The number of your total reads is 1.
# Number of total reads has to be equal or more than 2 ('minReadsNum' that you set)
# INFO [2024-09-13 19:18:43] ========================================================
# INFO [2024-09-13 19:18:43] ================ Creating 'SangerContig' ===============
# INFO [2024-09-13 19:18:43] ========================================================
# INFO [2024-09-13 19:18:43]   >> Contig Name: 'Achl_ACHLO041-09'
# SUCCESS [2024-09-13 19:18:44] --------------------------------------------------------
# SUCCESS [2024-09-13 19:18:44] -------- 'SangerRead' S4 instance is created !! --------
# SUCCESS [2024-09-13 19:18:44] --------------------------------------------------------
# SUCCESS [2024-09-13 19:18:44]    >> 'Achl_ACHLO041-09_1_F.ab1' is created (Forward Read; ABIF).
# INFO [2024-09-13 19:18:44]    >> The number of reads detected: 1
# WARN [2024-09-13 19:18:44] REGEX_MATCH_WARN
# Your 'contigName' and 'REGEX_SuffixReverse' regular expression parameters can not match any reverse reads.
# ERROR [2024-09-13 19:18:44] READ_NUMBER_ERROR
# The number of your total reads is 1.
# Number of total reads has to be equal or more than 2 ('minReadsNum' that you set)
# SUCCESS [2024-09-13 19:18:44] #############################################################
# SUCCESS [2024-09-13 19:18:44] ######## 'SangerAlignment' S4 instance is created !! ########
# SUCCESS [2024-09-13 19:18:44] #############################################################
# ERROR [2024-09-13 19:18:44] CONTIG_NUMBER_ZERO_ERROR
# The number of your total contig is 0.
# Please check your name matching parameters.
# DEBUG [2024-09-13 19:18:44]    >> For more information, please run 'object'.
# DEBUG [2024-09-13 19:18:44]    >> Run 'object@objectResults@readResultTable' to check the results of each Sanger reads
# > #launchApp(my_contigs)
# > writeFasta(my_contigs)
# INFO [2024-09-13 19:18:44] Your input is 'SangerAlignment' S4 instance
# INFO [2024-09-13 19:18:44] >>> outputDir : C:\Users\44139\AppData\Local\Temp\RtmpwRzaIz
# INFO [2024-09-13 19:18:44] Start to write 'SangerAlignment' to FASTA format ...
# INFO [2024-09-13 19:18:44] >> Writing 'alignment' to FASTA ...
# INFO [2024-09-13 19:18:44] >> Writing 'contigs' to FASTA ...
# INFO [2024-09-13 19:18:44] >> Writing all single reads to FASTA ...
# INFO [2024-09-13 19:18:44] Finish writing 'SangerAlignment' to FASTA format
# > generateReport(my_contigs)
# INFO [2024-09-13 19:18:44] Your input is 'SangerAlignment' S4 instance
# INFO [2024-09-13 19:18:44] >>> outputDir : C:\Users\44139\AppData\Local\Temp\RtmpwRzaIz
# INFO [2024-09-13 19:18:44] dir.exists(outputDirSA): TRUE
# 
# 
# processing file: SangerAlignment_Report.Rmd
# |.............................................................                                                     |  53% [unnamed-chunk-4]Error in `BrowseSeqs()`:
# ! No sequence information to display.
# Backtrace:
#   1. DECIPHER::BrowseSeqs(...)
# 
# Quitting from lines 123-127 [unnamed-chunk-4] (SangerAlignment_Report.Rmd)
# > 