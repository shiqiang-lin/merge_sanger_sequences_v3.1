#usage: python3.7 Merging_sequencing.py 000F.seq 001F.seq 002R.seq 003R.seq
"""
This script is used to merge multiple sucessive sanger DNA sequencing results.

The file names of the sequences to be merged starts with '000' plus 'F' (forward)
or 'R' (reverse), and are in .seq format.
The files are arranged in tandemly according to their real position corresponding
to the sequencing plasmid or linear DNA.
"""

import os,shutil
import sys
import string
from Bio.Emboss.Applications import NeedleCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio import SeqIO


for i in range(32):
    if i % 2 == 1:
        print(i*"(^_^)")
print("\n\n")


#get current path
current_path = os.getcwd()
print("Current path is %s." % current_path)
print("\n")

#make a new directory
folder_name = "merged_sequence"
dirs = os.listdir(current_path)
if folder_name not in dirs:
    os.mkdir(folder_name)
else:
    shutil.rmtree(folder_name)
    os.mkdir(folder_name)

number_of_files = len(sys.argv) -1
if len(sys.argv) -1 < 2:
    print("Must have 2 or more sequences. Please reinput sequences files to be merged.")
    sys.exit()
else:
    print("There are %s sequences to be joined." % number_of_files)

#copy sequence files to the new directory and change to F if it is in R
for file in sys.argv[1:]:
    if file[3] == 'F':
        print(file)
        file_F = open(file)
        file_sequence = file_F.read()
        file_F.close()
        DNA_sequence_tmp_list = []
        
        for i in file_sequence:
            if i.isalpha():
                i = i.upper()
                DNA_sequence_tmp_list.append(i)
        DNA_sequence_tmp_str =''.join(DNA_sequence_tmp_list)        
        dir_sub = os.path.join(current_path,folder_name)
        os.chdir(dir_sub)
        f = open(file[0:4] + '.seq','w')
        f.write(DNA_sequence_tmp_str)
        f.close()
        os.chdir(current_path)
        
    if file[3] == 'R':
        print(file)
        file_R = open(file)
        file_sequence = file_R.read()
        file_R.close()
        DNA_sequence_tmp_list = []
        
        for i in file_sequence:
            if i.isalpha():
                i = i.upper()
                DNA_sequence_tmp_list.append(i)
        DNA_sequence_tmp_str =''.join(DNA_sequence_tmp_list)
        DNA_sequence_tmp_Seq = Seq(DNA_sequence_tmp_str)
        DNA_sequence_tmp_Seq_F = DNA_sequence_tmp_Seq.reverse_complement()
        DNA_sequence_tmp_Seq_F_str = str(DNA_sequence_tmp_Seq_F)
        
        dir_sub = os.path.join(current_path,folder_name)
        os.chdir(dir_sub)
        f = open(file[0:4] + '.seq','w')
        f.write(DNA_sequence_tmp_Seq_F_str)
        f.close()
        os.chdir(current_path)
    

#function for the boundaries 
def needle_align_to_get_boundaries(f1,f2):
    aligment_a_left = 0
    aligment_b_left = 0
    output_file_name = f1.split('.')[0] + f2.split('.')[0] + ".needle"
    
    needle_cline = NeedleCommandline()
    needle_cline.asequence = f1
    needle_cline.bsequence = f2
    needle_cline.gapopen = 10
    needle_cline.gapextend= 0.5
    needle_cline.outfile = output_file_name
    print(needle_cline)
    stdout, stderr = needle_cline()
    print(stdout + stderr)

    #open the needle alignment output file and get boundaries
    file = open(output_file_name)
    file_lines = file.readlines()
    file.close()

    for line in file_lines:
        print(line, end="")

    aligment_a_squence_positions = []
    aligment_b_squence_positions = []

    file = open(output_file_name)

    new_line1 = file.readline()
    new_line2 = file.readline()

    while len(new_line2):
        line_a = new_line1
        line_b = new_line2
        new_line2 = new_line2.strip()
        
        if (50*'|' in new_line2):

            line_b = file.readline()
      
            aligment_a_squence_line_str = line_a.strip()
            aligment_b_squence_line_str = line_b.strip()
            print("The beginning of excellent alignment is shown below.\n")
            
            aligment_a_squence_line_str_split = aligment_a_squence_line_str.split()
            print(aligment_a_squence_line_str_split[0].ljust(5,' '),\
                  aligment_a_squence_line_str_split[1],\
                  aligment_a_squence_line_str_split[2].rjust(6,' '),\
                  sep="")
            aligment_b_squence_line_str_split = aligment_b_squence_line_str.split()
            print(aligment_b_squence_line_str_split[0].ljust(5,' '),\
                  aligment_b_squence_line_str_split[1],\
                  aligment_b_squence_line_str_split[2].rjust(6,' '),\
                  sep="")
     
            print("\n")

            aligment_a_left = int(aligment_a_squence_line_str.split()[0])
            aligment_b_left = int(aligment_b_squence_line_str.split()[0])
            break
        else:
            new_line1 = new_line2         #notice the skill here, one must go step by step through the lines
            new_line2 = file.readline()
    file.close()

    return aligment_a_left, aligment_b_left

#end of function

#align with needle progressively to get boundaries
dir_sub = os.path.join(current_path,folder_name)
os.chdir(dir_sub)
dir_new = os.listdir()
dir_new.sort()

Es_list = []
Es_list.append(1)

for i in range(len(dir_new) - 1):
    E = needle_align_to_get_boundaries(dir_new[i],dir_new[i+1])
    Es_list.append(E[0])
    Es_list.append(E[1])

# length of the last sequence

f1 = open(dir_new[-1])
file_sequence = f1.read()
f1.close()
DNA_sequence_tmp_list = []
for i in file_sequence:
    if i.isalpha():
        i = i.upper()
        DNA_sequence_tmp_list.append(i)
length_of_last_file_sequence = len(DNA_sequence_tmp_list)

Es_list.append(length_of_last_file_sequence + 1)

#till now, got the boundaries for each sequence

#list to store the final merged sequence
merged_sequence_list = []        
        

for i in range(len(dir_new)):
    f = open(dir_new[i])
    file_sequence = f.read()
    f.close()
    DNA_sequence_tmp_list = []
    for j in file_sequence:
                if j.isalpha():
                    j = j.upper()
                    DNA_sequence_tmp_list.append(j)
    merged_sequence_list += DNA_sequence_tmp_list[Es_list[2*i]-1:Es_list[2*i+1]-1]
"""
The content of Es_list, 2n elements in toal
E0=1,E1,    E2,E3,  E4,E5,  E6,E7   ...     E2n-2,E2n-1
S0          S1      S2      S3              Sn-1
mathematical [X,Y)~~mathematical [X,Y-1]~~list elements from X-1 to Y-2~~list [X-1:Y-1]
"""

merged_sequence_str = ''.join(merged_sequence_list)
print("The merged DNA sequence is shown below.\n")
for j in range(0,len(merged_sequence_str),100):
    DNA_string_100_per_line_str = merged_sequence_str[j:j+100]
    print(DNA_string_100_per_line_str)
print("\n")

#write to the file
merged_file_name = "merged" + ".seq"
file=open(merged_file_name,'w')
file.write(merged_sequence_str)
file.close()
os.system('rm *.needle')
os.system('rm 00*')


