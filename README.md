### A Python program to merge Sanger sequences: an update
This program is used to merge multiple DNA sequence files. The DNA sequence file supports seq/txt/fasta formats, and the program can automatically recognize the file format and read the sequence content correctly. This program can obtain overlapping regions of multiple DNA sequences to get the correct merging result. The program can run on macOS/Linux/Windows and the results are the same. Please follow the instructions below to install and run the program.

The program supports two modes: GUI (Graphical User Interface) and CLI (Command Line Interface).


### 1. Installation
#### <font color=red> 1.1 in macOS/Linux </font>

(1) make sure you have installed Python 3.12.

```zsh
% python3.12
Python 3.12.1 (v3.12.1:2305ca5144, Dec  7 2023, 17:23:38) [Clang 13.0.0 (clang-1300.0.29.30)] on darwin
Type "help", "copyright", "credits" or "license" for more information.
```
or
```zsh
% python
Python 3.12.1 (v3.12.1:2305ca5144, Dec  7 2023, 17:23:38) [Clang 13.0.0 (clang-1300.0.29.30)] on darwin
Type "help", "copyright", "credits" or "license" for more information.
```
Please note that some computers may have more than one Python version installed. You need to check the version number of the Python. If you enter "python3.12" in the command prompt and do not find it, you can try to use the "python" command. The following steps are based on "python3.12" commands, and you can decide whether to change them to "python" commands based on your specific environment.

(2) make sure you have installed the dependent package biopython 1.83 and numpy 2.1.1.

```zsh
% python3.12 -m pip list
biopython                                1.83
numpy                                    2.1.1
pip                                      24.2
```

The above output indicates that the biopython 1.83 and numpy 2.1.1 package has been installed.
If the dependency is not installed, install it with the following command:

```zsh
% python3.12 -m pip install biopython==1.83
```
If you have installed the biopython package but the version is inconsistent, run the following command to update the version of biopython:
```zsh
% python3.12 -m pip install --upgrade biopython==1.83
```
Confirm whether the version is correct after the update is complete.

```zsh
% python3.12 -m pip show biopython
Name: biopython
Version: 1.83
Summary: Freely available tools for computational molecular biology.
Home-page: https://biopython.org/
Author: The Biopython Contributors
Author-email: biopython@biopython.org
```

In most cases, pip will automatically install numpy when installing biopython, but if pip does not automatically install numpy in your environment, you need to install numpy manually.


```zsh
% python3.12 -m pip list
biopython                                1.83
pip                                      24.2
```

If numpy is not found, execute the following command.

```zsh
% python3.12 -m pip install numpy==2.1.1
```

(3) Download the program file Merge_Sanger_v3.1.py and example sequence files to your computer from GitHub https://github.com/shiqiang-lin/merge_sanger_sequences_v3.1.

(4) Make a new directory on the desktop.

```zsh
% cd ~              # Go to the user's Home path
% cd Desktop        # Go to the user desktop
% mkdir MerSanSeqs  # folder named MerSanSeqs. You can name it something else
```
(5) Copy the downloaded Merge_Sanger_v3.1.py to MerSanSeq, either by dragging or right-clicking.

#### <font color=red> 1.2 in Windows </font>
(1) make sure you have installed Python 3.12.

```zsh
> python3.12
Python 3.12.1 (v3.12.1:2305ca5144, Dec  7 2023, 17:23:38) [MSC v.1940 64 bit (AMD64)] on win32
Type "help", "copyright", "credits" or "license" for more information.
```
or
```zsh
> python
Python 3.12.1 (v3.12.1:2305ca5144, Dec  7 2023, 17:23:38) [MSC v.1940 64 bit (AMD64)] on win32
Type "help", "copyright", "credits" or "license" for more information.
```

Please note that some computers may have more than one Python version installed. You need to check the version number of the Python. If you enter "python3.12" in the command prompt and do not find it, you can try to use the "python" command. The following steps are based on "python3.12" commands, and you can decide whether to change them to "python" commands based on your specific environment.

(2) make sure you have installed the dependent package biopython 1.83 and numpy 2.1.1.

```zsh
> python3.12 -m pip list
biopython 1.83  
numpy     2.1.1
pip       24.2    
```

The above output indicates that the biopython 1.83 and numpy 2.1.1 package has been installed.
If the dependency is not installed, install it with the following command.

```zsh
> python3.12 -m pip install biopython==1.83
```
If you have installed the biopython package but the version is inconsistent, run the following command to update the version of biopython.

```zsh
> python3.12 -m pip install --upgrade biopython==1.83
```
Confirm whether the version is correct after the update is complete.

```zsh
> python3.12 -m pip show biopython
Name: biopython
Version: 1.83
Summary: Freely available tools for computational molecular biology.
Home-page: https://biopython.org/
Author: The Biopython Contributors
Author-email: biopython@biopython.org
```

In most cases, pip will automatically install numpy when installing biopython, but if pip does not automatically install numpy in your environment, you need to install numpy manually.


```zsh
> python3.12 -m pip list
biopython 1.83
pip       24.2
```

If numpy is not found, execute the following command.

```zsh
> python3.12 -m pip install numpy==2.1.1
```

(3) Download the program file Merge_Sanger_v3.1.py and example sequence files to your computer from GitHub https://github.com/shiqiang-lin/merge_sanger_sequences_v3.1.


(4) Right-click on the desktop to create a new folder: MerSanSeq(can be other name), and copy the downloaded Merge_Sanger_v3.1.py to MerSanSeq, which can be copied and pasted by dragging or right-clicking.

### 2. Usage (GUI mode)
Here are the steps to start the program (GUI mode):

#### <font color=red> 2.1 in macOS/Linux </font>
Open the terminal, go to the Desktop/MerSanSeqs directory, and launch the program.

```zsh
% cd ~
% cd Desktop/MerSanSeqs
% python3.12 Merge_Sanger_v3.1.py
```
or
```zsh
% cd ~
% cd Desktop/MerSanSeqs
% python3.12 Merge_Sanger_v3.1.py -gui
```
#### <font color=red> 2.2 in Windows </font>
Use cmd.exe to open the Command Prompt window. Go to the C:/Users/44139/Desktop/MerSanSeq directory (this directory is the author's environment, please modify it according to the actual path), and start the program.

```zsh
> cd C:/Users/44139/Desktop/MerSanSeq
> python3.12 Merge_Sanger_v3.1.py
```
or
```zsh
> cd C:/Users/44139/Desktop/MerSanSeq
> python3.12 Merge_Sanger_v3.1.py -gui
```

#### <font color=red> 2.3 in X11 forwarding </font>

The GUI mode of our program supports X11 forwarding to connect to a remote Linux server. Here we show how to use a Windows 11 client to connect to a remote Linux server using the SecureCRT client.

(1) Install the Xming program on Windows 11;

(2) Use SecureCRT to log into the remote Linux server with SSH2; 

(3) After logging into the Linux server, set the environment variable DISPLAY. 

The 172.27.55.36 is the IP address of the author’s Windows computer. Please replace it with the IP address of your own computer.


```zsh
% export DISPLAY=172.27.55.36:0.0
% 
```

(4) Run the xclock program to see if the X11 forwarding configuration is normal, if the xclock clock can be displayed, the configuration is correct.

```zsh
% xclock
% 
```
(5) Open the terminal, go to the Desktop/MerSanSeqs directory, and launch the program.

```zsh
% cd ~
% cd Desktop/MerSanSeqs
% python3.12 Merge_Sanger_v3.1.py
```

After launching the program, the program will open the main GUI interface, and here are the steps to use.

#### <font color=red> 2.4 in macOS/Linux/Windows </font>

(1) Click the Add File button to add example sequences one by one.
Select one file at a time. You need to add at least two sequence files. 
Please pay attention to the order in which the files are added. The program will align and merge them in the order in which they are added. 
If you want to delete a sequence file, click the Remove button. If there are less than two files, the program will pop up a prompt when merging.

(2) Select the sequence direction F/R.
Please select the sequence direction F/R correctly. The wrong choice will lead to wrong results. If one or more sequences do not have the F/R attribute selected, the program will pop up a prompt when merging.

(3) Set the number of Consecutive Matches, the default value is 50.
This value can be left unset, the default value is 50. If you need to use different consecutive matches, you can set this value. Set the range to [30, 60], beyond which the program will use the default value of 50.

(4) Click the Merge button to merge the sequence.
The program will align the added sequence files, and display the log in the log area on the right of GUI. During the running process, the program will create a merged_sequence& (& stands for a timestamp) subdirectory in the directory where the first sequence file is located. The log, process files, and merged result file, are saved in the  erged_sequence*, and all of these files can be opened using a text editor.

(5) Click the Save Result button to save the merged sequence.
To save the merged sequence to a specified directory, you can select one of the three file formats: *.seq, *.fasta, *.txt.

(6) Click the Reset button to clear the loaded sequence files and log information.
The program can be used multiple times without the need to close-open again. Clicking the Reset button will erase all the files and results of the previous group, and you can continue to join the next set of file sequences.

### 3. Usage (CLI mode)
Here are the steps to start the program(CLI mode):

#### <font color=red> 3.1 in macOS/Linux </font>
Open the terminal, go to the Desktop/MerSanSeqs directory, and launch the program.


There are two ways to run the CLI mode.
```zsh
% cd ~
% cd Desktop/MerSanSeqs
% python3.12 Merge_Sanger_v3.1.py -cli test.txt
```

The test.txt contains the file names of the sequence files to be merged and their F/R attributes. The content of the test.txt is shown as follows.
```zsh
% more test.txt
F104120_F714961b_CS13_T7.seq F
F103749_F714961c_CS13_CS13.W1F(B11369.seq F
```
or
```zsh
% cd ~
% cd Desktop/MerSanSeqs
% python3.12 Merge_Sanger_v3.1.py -cli F104120_F714961b_CS13_T7.seq F "F103749_F714961c_CS13_CS13.W1F(B11369.seq" F 
```
If the file name contains brackets, please enclose the file name in quotes.

After launching the program, the program will merge the sequences and show logs in terminal.

#### <font color=red> 3.2 in Windows </font>
Use cmd.exe to open the Command Prompt window. Go to the C:/Users/44139/Desktop/MerSanSeq directory (this directory is the author's environment, please modify it according to the actual path), and start the program.

There are two ways to run the CLI (command-line interface) mode.
```zsh
> cd C:/Users/44139/Desktop/MerSanSeq
> python3.12 Merge_Sanger_v3.1.py -cli test.txt
```
The test.txt contains the file names of the sequence files to be merged and their F/R attributes. The content of the test.txt is shown as follows.
```zsh
> more test.txt
F104120_F714961b_CS13_T7.seq F
F103749_F714961c_CS13_CS13.W1F(B11369.seq F
```
or
```zsh
> cd C:/Users/44139/Desktop/MerSanSeq
> python3.12 Merge_Sanger_v3.1.py -cli F104120_F714961b_CS13_T7.seq F "F103749_F714961c_CS13_CS13.W1F(B11369.seq" F
```
If the file name contains brackets, please enclose the file name in quotes.

After launching the program, the program will merge the sequences and show logs in Command Prompt window.



