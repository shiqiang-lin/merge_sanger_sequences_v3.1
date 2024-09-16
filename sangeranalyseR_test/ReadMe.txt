Citation: Kuan-Hao Chao, Kirston Barton, Sarah Palmer, and Robert Lanfear (2021). "sangeranalyseR: simple and interactive processing of Sanger sequencing data in R" in Genome Biology and Evolution. DOI: doi.org/10.1093/gbe/evab028

Test procedure:

1. Install R and Rstudio
2. Download source code of sangeranalyseR from https://github.com/roblanf/sangeranalyseR
3. Unzip source code and get a directory 'sangeranalyseR-master'
4. Make a new directory 't'
5. Copy the four F (forward direction) ab1 files in 'sangeranalyseR-master/inst/extdata/Allolobophora_chlorotica/ACHLO' to directory 't'.
This step can also be done by download the four F ab1 files from https://github.com/roblanf/sangeranalyseR/tree/master/inst/extdata/Allolobophora_chlorotica/ACHLO.
6. Modify the dataDir in the test_script.R to the path where directory 't' resides
7. Run the script test_script.R in RStudio
8. The output in the console is shown below. The directory 'RtmpwRzaIz' is where the sangeranalyseR puts its output content. 

INFO [2024-09-13 19:18:43] #################################################
INFO [2024-09-13 19:18:43] #### Start creating SangerAlignment instance ####
INFO [2024-09-13 19:18:43] #################################################
INFO [2024-09-13 19:18:43]   >> You are using Regular Expression Method to group AB1 files!
INFO [2024-09-13 19:18:43] ========================================================
INFO [2024-09-13 19:18:43] ================ Creating 'SangerContig' ===============
INFO [2024-09-13 19:18:43] ========================================================
INFO [2024-09-13 19:18:43]   >> Contig Name: 'Achl_ACHLO006-09'
SUCCESS [2024-09-13 19:18:43] --------------------------------------------------------
SUCCESS [2024-09-13 19:18:43] -------- 'SangerRead' S4 instance is created !! --------
SUCCESS [2024-09-13 19:18:43] --------------------------------------------------------
SUCCESS [2024-09-13 19:18:43]    >> 'Achl_ACHLO006-09_1_F.ab1' is created (Forward Read; ABIF).
INFO [2024-09-13 19:18:43]    >> The number of reads detected: 1
WARN [2024-09-13 19:18:43] REGEX_MATCH_WARN
Your 'contigName' and 'REGEX_SuffixReverse' regular expression parameters can not match any reverse reads.
ERROR [2024-09-13 19:18:43] READ_NUMBER_ERROR
The number of your total reads is 1.
Number of total reads has to be equal or more than 2 ('minReadsNum' that you set)
INFO [2024-09-13 19:18:43] ========================================================
INFO [2024-09-13 19:18:43] ================ Creating 'SangerContig' ===============
INFO [2024-09-13 19:18:43] ========================================================
INFO [2024-09-13 19:18:43]   >> Contig Name: 'Achl_ACHLO007-09'
SUCCESS [2024-09-13 19:18:43] --------------------------------------------------------
SUCCESS [2024-09-13 19:18:43] -------- 'SangerRead' S4 instance is created !! --------
SUCCESS [2024-09-13 19:18:43] --------------------------------------------------------
SUCCESS [2024-09-13 19:18:43]    >> 'Achl_ACHLO007-09_1_F.ab1' is created (Forward Read; ABIF).
INFO [2024-09-13 19:18:43]    >> The number of reads detected: 1
WARN [2024-09-13 19:18:43] REGEX_MATCH_WARN
Your 'contigName' and 'REGEX_SuffixReverse' regular expression parameters can not match any reverse reads.
ERROR [2024-09-13 19:18:43] READ_NUMBER_ERROR
The number of your total reads is 1.
Number of total reads has to be equal or more than 2 ('minReadsNum' that you set)
INFO [2024-09-13 19:18:43] ========================================================
INFO [2024-09-13 19:18:43] ================ Creating 'SangerContig' ===============
INFO [2024-09-13 19:18:43] ========================================================
INFO [2024-09-13 19:18:43]   >> Contig Name: 'Achl_ACHLO040-09'
SUCCESS [2024-09-13 19:18:43] --------------------------------------------------------
SUCCESS [2024-09-13 19:18:43] -------- 'SangerRead' S4 instance is created !! --------
SUCCESS [2024-09-13 19:18:43] --------------------------------------------------------
SUCCESS [2024-09-13 19:18:43]    >> 'Achl_ACHLO040-09_1_F.ab1' is created (Forward Read; ABIF).
INFO [2024-09-13 19:18:43]    >> The number of reads detected: 1
WARN [2024-09-13 19:18:43] REGEX_MATCH_WARN
Your 'contigName' and 'REGEX_SuffixReverse' regular expression parameters can not match any reverse reads.
ERROR [2024-09-13 19:18:43] READ_NUMBER_ERROR
The number of your total reads is 1.
Number of total reads has to be equal or more than 2 ('minReadsNum' that you set)
INFO [2024-09-13 19:18:43] ========================================================
INFO [2024-09-13 19:18:43] ================ Creating 'SangerContig' ===============
INFO [2024-09-13 19:18:43] ========================================================
INFO [2024-09-13 19:18:43]   >> Contig Name: 'Achl_ACHLO041-09'
SUCCESS [2024-09-13 19:18:44] --------------------------------------------------------
SUCCESS [2024-09-13 19:18:44] -------- 'SangerRead' S4 instance is created !! --------
SUCCESS [2024-09-13 19:18:44] --------------------------------------------------------
SUCCESS [2024-09-13 19:18:44]    >> 'Achl_ACHLO041-09_1_F.ab1' is created (Forward Read; ABIF).
INFO [2024-09-13 19:18:44]    >> The number of reads detected: 1
WARN [2024-09-13 19:18:44] REGEX_MATCH_WARN
Your 'contigName' and 'REGEX_SuffixReverse' regular expression parameters can not match any reverse reads.
ERROR [2024-09-13 19:18:44] READ_NUMBER_ERROR
The number of your total reads is 1.
Number of total reads has to be equal or more than 2 ('minReadsNum' that you set)
SUCCESS [2024-09-13 19:18:44] #############################################################
SUCCESS [2024-09-13 19:18:44] ######## 'SangerAlignment' S4 instance is created !! ########
SUCCESS [2024-09-13 19:18:44] #############################################################
ERROR [2024-09-13 19:18:44] CONTIG_NUMBER_ZERO_ERROR
The number of your total contig is 0.
Please check your name matching parameters.
DEBUG [2024-09-13 19:18:44]    >> For more information, please run 'object'.
DEBUG [2024-09-13 19:18:44]    >> Run 'object@objectResults@readResultTable' to check the results of each Sanger reads
> #launchApp(my_contigs)
> writeFasta(my_contigs)
INFO [2024-09-13 19:18:44] Your input is 'SangerAlignment' S4 instance
INFO [2024-09-13 19:18:44] >>> outputDir : C:\Users\44139\AppData\Local\Temp\RtmpwRzaIz
INFO [2024-09-13 19:18:44] Start to write 'SangerAlignment' to FASTA format ...
INFO [2024-09-13 19:18:44] >> Writing 'alignment' to FASTA ...
INFO [2024-09-13 19:18:44] >> Writing 'contigs' to FASTA ...
INFO [2024-09-13 19:18:44] >> Writing all single reads to FASTA ...
INFO [2024-09-13 19:18:44] Finish writing 'SangerAlignment' to FASTA format
> generateReport(my_contigs)
INFO [2024-09-13 19:18:44] Your input is 'SangerAlignment' S4 instance
INFO [2024-09-13 19:18:44] >>> outputDir : C:\Users\44139\AppData\Local\Temp\RtmpwRzaIz
INFO [2024-09-13 19:18:44] dir.exists(outputDirSA): TRUE


processing file: SangerAlignment_Report.Rmd
  |.............................................................                                                     |  53% [unnamed-chunk-4]Error in `BrowseSeqs()`:
! No sequence information to display.
Backtrace:
 1. DECIPHER::BrowseSeqs(...)

Quitting from lines 123-127 [unnamed-chunk-4] (SangerAlignment_Report.Rmd)
                                                                                                                                                                          
> 



=========================================================================================================
~~~***~~~***~~~***~~~***~~~***~~~***   End of Console Output   ~~~***~~~***~~~***~~~***~~~***~~~***~~~***  
=========================================================================================================



Positions of Error report of the program

https://github.com/roblanf/sangeranalyseR/blob/master/R/ClassSangerContig.R

Line 689-691:
        if (readNumber >= minReadsNum) {
            msg <- ""
            if (readNumber >= 2) {

From these lines, we know that for each contig, the sum of forwardNumber + reverseNumber should be no less than max(minReadsNum,2).

Line 767: msg <- paste0("The number of your total reads is ", readNumber, ".",
produce Error report we see in the Console output shown above. We pick one of them and shown below.

ERROR [2024-09-13 19:18:43] READ_NUMBER_ERROR
The number of your total reads is 1.


Also, the three .fa files in the output folder RtmpwRzaIz are all empty, indicating that sangeranalyseR cannot deal with this kind of Single-End sequencing data.



Conclusion
It can be seen from above results that the program mainly analyses the Paired-end Sequencing data. It is better suited for the analysis of gene diversity such as SNP (single nucleotide polymorphism), indel (insertion and deletions ), and so on. In the sequencing confirmation of gene sequence during gene cloning, Single-end sequencing is often used. The sequencing data of this kind cannot be processed with sangeranalyseR since it requires Paired-end Sequencing data. Therefore, we do not add the comparison of MergeSangerV3.1 with sangeranalyseR in the revised manuscript.




