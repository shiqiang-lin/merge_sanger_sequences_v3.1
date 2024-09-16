import sys
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext
import os
import shutil
from Bio.Seq import Seq
from Bio import Align
from datetime import datetime

class MultiFileGUI:
    def __init__(self, master):
        self.master = master
        master.title("MergeSanger V3.1")

        # menu
        self.menu_bar = tk.Menu(self.master)
        self.master.config(menu=self.menu_bar)

        # About
        self.about_menu = tk.Menu(self.menu_bar, tearoff=0)
        self.menu_bar.add_cascade(label="About", menu=self.about_menu)
        self.about_menu.add_command(label="About", command=self.show_about)

        self.master.config(bg="#f0f0f0")
        width = root.winfo_screenwidth()
        height = root.winfo_screenheight()
        root.geometry(f"{width}x{height}+0+0")

        self.buttons_frame = tk.Frame(master, bg="#d0d0d0")
        self.buttons_frame.grid(row=0, column=0, sticky="wn", padx=10, pady=10)

        self.open_button = tk.Button(self.buttons_frame, text="Add File", font=("Arial", 10), command=self.open_files,
                                     bg="#a0a0a0")
        self.open_button.grid(row=0, column=0, padx=10, pady=5)

        self.merge_button = tk.Button(self.buttons_frame, text="Merge", font=("Arial", 10), command=self.merge_sanger,
                                      bg="#a0a0a0")
        self.merge_button.grid(row=0, column=1, padx=10, pady=5)

        self.clear_button = tk.Button(self.buttons_frame, text="Reset", font=("Arial", 10), command=self.reset,
                                      bg="#a0a0a0")
        self.clear_button.grid(row=0, column=2, padx=10, pady=5)

        self.save_button = tk.Button(self.buttons_frame, text="Save Result", font=("Arial", 10), command=self.save_file,
                                     bg="#a0a0a0")
        self.save_button.grid(row=0, column=3, padx=10, pady=5)

        self.label = tk.Label(self.buttons_frame, text="Consecutive Matches:", font=("Arial", 10), bg="#d0d0d0")
        self.label.grid(row=0, column=4, padx=5)

        var = tk.StringVar(value="50")
        self.entry_area = tk.Entry(self.buttons_frame, textvariable=var, font=("Arial", 10), width=10)
        self.entry_area.grid(row=0, column=5, padx=5)

        self.text_area = scrolledtext.ScrolledText(master, wrap=tk.WORD, font=("Courier New", 10), bg="#ffffff")
        self.text_area.grid(row=0, column=1,  sticky="nsew", rowspan=500, padx=10, pady=10)
        master.grid_rowconfigure(0, weight=1)
        master.grid_columnconfigure(1, weight=1)

        self.frames = []
        self.merge_result = ""
        self.dna_sequence = []
        self.logdir = ""

    @staticmethod
    def show_about():
        # version
        messagebox.showinfo("About", "MergeSanger Version: 3.1\nPlease send email to 441394376@qq.com\n"
                                     "if you encounter any bugs\n09.2024")
    def remove_frame(self, frame_tuple):
        """
        remove one file
        """
        # frame_tuple[0].grid_forget()
        frame_tuple[0].destroy()

        # refresh the frames
        m = frame_tuple[3]

        for frame_t in self.frames[m:]:
            frame = frame_t[0]
            newrow = frame_t[3] - 1
            # print(f"move frame from row {frame_t[3]} to row {newrow}")
            frame.grid(row=newrow, column=0, sticky='wn', pady=5)
            frame_t[3] = newrow

        self.frames.remove(frame_tuple)

    def check_file_list_count(self):
        """
        check if the file list contains 2 or more files
        """
        file_count = len(self.frames)
        if file_count >= 2:
            print(f"The file list contains {file_count} files, which is greater than or equal to 2.")
        else:
            print(f"The file list does not meet the criteria of having 2 or more files.")
            messagebox.showerror(
                "Error",
                f"The file list contains {file_count} files, which does not meet the criteria of having 2 or more files."
            )
            return -1

    def check_radiobutton_selection(self):
        """
        check if all of the F or R directions are selected
        """
        unselected_groups = []
        for i, frame_tuple in enumerate(self.frames):
            if frame_tuple[2].get() == "none":
                unselected_groups.append(f"Group {i + 1}")

        if unselected_groups:
            messagebox.showerror(
                "Error",
                f"Please make a selection in the following radiobutton groups: {', '.join(unselected_groups)}."
            )
            return -1
        else:
            print("All radiobutton groups have a selection.")

    def add_file_frame(self, filepath, content):
        """
        add one file, shows the filepath ,F/R direction and remove button
        """
        if len(self.frames) % 2 == 0:
            frame = tk.Frame(root, bg="#d0d0d0")
        else:
            frame = tk.Frame(root)
        frame_row = len(self.frames) + 1

        # set log dir
        if frame_row == 1:
            self.logdir = os.path.dirname(filepath)

        # print(f"row={frame_row}")
        frame.grid(row=frame_row, column=0, sticky='wn', pady=5)
        self.master.grid_rowconfigure(frame_row-1, weight=0)
        self.master.grid_rowconfigure(frame_row, weight=1)
        frame.grid_columnconfigure(2, weight=1)

        new_filepath = '\n'.join([filepath[i:i + 80] for i in range(0, len(filepath), 80)])
        filepath_var = tk.StringVar(value=new_filepath)

        label = tk.Label(frame, textvariable=filepath_var, anchor='w', justify='left', font=("Arial", 10), fg="red")
        label.grid(row=0, column=0, columnspan=3, sticky='wn', padx=(0, 5), pady=(5, 0))
        option_var = tk.StringVar(value="none")
        radio1 = tk.Radiobutton(frame, text="F", variable=option_var, value="F", font=("Arial", 10))
        radio1.grid(row=1, column=0, sticky='w', pady=(0, 5))
        radio2 = tk.Radiobutton(frame, text="R", variable=option_var, value="R", font=("Arial", 10))
        radio2.grid(row=1, column=1, sticky='w', pady=(0, 5))

        frame_tuple = [frame, filepath, option_var, frame_row]
        self.frames.append(frame_tuple)

        remove_btn = tk.Button(frame, text='Remove', font=("Arial", 10), bg="#a0a0a0",
                               command=lambda f=frame_tuple: self.remove_frame(f))
        remove_btn.grid(row=1, column=2, sticky='w')

    def open_files(self):
        """
        open one file
        """
        file_path = filedialog.askopenfilename()
        if file_path:
            with open(file_path, 'r') as file:
                content = file.read()
                self.text_area.insert(tk.END, f'\n### File: {file_path} ###\n')
                self.text_area.insert(tk.END, content + '\n')
                self.add_file_frame(file_path, content)

    def reset(self):
        """
        reset program, remove files, clear text area
        """
        self.text_area.delete('1.0', tk.END)
        for frame in self.frames:
            frame[0].destroy()
        self.frames.clear()
        self.dna_sequence.clear()
        self.merge_result = ""

    @staticmethod
    def read_sequence_from_file(file_path):
        """
        read sequence file, support seq, txt and fasta format
        """
        sequences = []
        current_sequence = ''
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()
                if not line.startswith('>'):
                    current_sequence += line
                else:
                    if current_sequence:
                        sequences.append(current_sequence)
                        current_sequence = ''
            if current_sequence:
                sequences.append(current_sequence)
            sequence_str = ''.join([i.upper() for i in sequences if i.isalpha()])
        return sequence_str

    @staticmethod
    def perform_needleman_wunsch_pairwisealigner(seq1, seq2, gap_open, gap_extend):
        """
        align two sequences using Align.PairwiseAligner() function. Be mindful of these parameter settings
        """
        # Create a PairwiseAligner object
        aligner = Align.PairwiseAligner()

        # Set parameters
        aligner.mode = 'global'
        aligner.match_score = 5
        aligner.mismatch_score = -4
        aligner.open_gap_score = -gap_open
        aligner.extend_gap_score = -gap_extend
        aligner.end_gap_score = 0
        # Perform alignment
        alignments = aligner.align(seq1, seq2)

        best_alignment = alignments[0]
        return best_alignment

    @staticmethod
    def find_n_consecutive(aligned_seq1, aligned_seq2, n):
        """
        find the position of n consecutive matching fragment
        """
        if len(aligned_seq1) != len(aligned_seq2):
            print("length of aligned seq1 must equal to aligned seq2.")
            return -1
        for i in range(len(aligned_seq1) - n + 1):
            if aligned_seq1[i:i + n] == aligned_seq2[i:i + n]:
                return i
        return -1

    @staticmethod
    def format_needle_srspair(seq1, seq2, alignment):
        """
        format the alignment by needle file format
        """
        formatted_alignment = []

        fasta = alignment.format("fasta")
        seq_new = fasta.split('>', 3)
        aligned_seq1 = seq_new[1].strip()
        aligned_seq2 = seq_new[2].strip()

        pos_begin1 = 0
        pos_end1 = 0
        length = len(aligned_seq1)

        pos_begin2 = 0
        pos_end2 = 0

        lineLength = 50
        for i in range(0, length, lineLength):
            line_seq1 = aligned_seq1[i:i + lineLength]
            line_seq2 = aligned_seq2[i:i + lineLength]

            if pos_end1 < len(seq1):
                pos_begin1 = pos_end1 + 1
            else:
                pos_begin1 = pos_end1

            if pos_end2 < len(seq2):
                pos_begin2 = pos_end2 + 1
            else:
                pos_begin2 = pos_end2
            pos_indicator = ""

            for j in range(len(line_seq1)):
                if line_seq1[j] == line_seq2[j]:
                    pos_indicator += "|"
                elif line_seq1[j] == '-' or line_seq2[j] == '-':
                    pos_indicator += " "
                else:
                    pos_indicator += "."

                if line_seq1[j] != '-':
                    pos_end1 += 1
                if line_seq2[j] != '-':
                    pos_end2 += 1

            if pos_end1 == 0:
                pos_begin1 = 0
            if pos_end2 == 0:
                pos_begin2 = 0
            formatted_alignment.append(
                f"{pos_begin1:<4}  {line_seq1}  {pos_end1}\n      {pos_indicator}\n{pos_begin2:<4}  {line_seq2}  {pos_end2}\n\n")

        return ''.join(formatted_alignment)

    @staticmethod
    def format_alignment(seq1, seq2, alignment):
        """
        string to print
        """
        alignment_str = ""
        alignment_str += "###########################################\n"
        alignment_str += "#=========================================#\n"
        alignment_str += "#   Aligned_sequences: 1\n"
        alignment_str += "#=========================================#\n"
        alignment_str += "###########################################\n\n"
        alignment_str += f"Sequence: 1;\n{seq1}\n"
        alignment_str += f"Length: {len(seq1)}  Type: N  Check: 0  ..\n"
        alignment_str += "Start: 1  End: {}\n\n".format(len(seq1))
        alignment_str += "###########################################\n"
        alignment_str += "#=========================================#\n"
        alignment_str += "#   Aligned_sequences: 2\n"
        alignment_str += "#=========================================#\n"
        alignment_str += "###########################################\n\n"
        alignment_str += f"Sequence: 2;\n{seq2}\n"
        alignment_str += f"Length: {len(seq2)}  Type: N  Check: 0  ..\n"
        alignment_str += "Start: 1  End: {}\n\n".format(len(seq2))
        alignment_str += "###########################################\n"
        alignment_str += "#=========================================#\n"
        alignment_str += "#   Pairwise alignments\n"
        alignment_str += "#=========================================#\n"
        alignment_str += "###########################################\n\n"
        alignment_str += f"###########################################\n"
        alignment_str += "#                                         #\n"
        alignment_str += f"#  Alignment score: {alignment.score:.1f}\n"
        alignment_str += "#                                         #\n"
        alignment_str += "###########################################\n\n"

        return alignment_str

    @staticmethod
    def check_and_convert(s):
        """
        make sure value of s between 30 and 60
        """
        try:
            num = int(s)
            if num < 30 or num > 60:
                return 50
            else:
                return num
        except ValueError:
            return 50

    def alignment_two_seq(self, seq1, seq2, save_file):
        """
        align seq1 and seq2, save alignment file
        """
        gap_open = 10.0
        gap_extend = 0.5

        # Perform Needleman-Wunsch alignment
        alignment = self.perform_needleman_wunsch_pairwisealigner(seq1, seq2, gap_open, gap_extend)
        # Format alignment for output
        alignment_output = self.format_alignment(seq1, seq2, alignment)

        # Print or save alignment output
        print(alignment_output)
        self.text_area.insert(tk.END, alignment_output + "\n")

        #get aligned seq1 and seq2
        fasta = alignment.format("fasta")
        # print(fasta)

        seq_new = fasta.split('>', 3)
        aligned_seq1 = seq_new[1].strip()
        aligned_seq2 = seq_new[2].strip()

        entry_value = self.entry_area.get()
        n = 50

        n = self.check_and_convert(entry_value)
        print("n=", n)

        aligned_pos = self.find_n_consecutive(aligned_seq1, aligned_seq2, n)
        if aligned_pos == -1:
            print(f"can not find consecutive {n} equal in two sequences.")
            self.text_area.insert(tk.END, f"can not find consecutive {n} equal in two sequences." + "\n")
            print(self.format_needle_srspair(seq1, seq2, alignment))
            self.text_area.insert(tk.END, self.format_needle_srspair(seq1, seq2, alignment) + "\n")
            f = open(save_file, 'w')
            f.write(alignment_output)
            f.write(self.format_needle_srspair(seq1, seq2, alignment))
            f.close()
            return -1, -1
        seq1_pos = aligned_pos - aligned_seq1[:aligned_pos].count('-') + 1 # start from 1
        seq2_pos = aligned_pos - aligned_seq2[:aligned_pos].count('-') + 1

        print(self.format_needle_srspair(seq1, seq2, alignment))
        self.text_area.insert(tk.END, self.format_needle_srspair(seq1, seq2, alignment) + "\n")
        f = open(save_file, 'w')
        f.write(alignment_output)
        f.write(self.format_needle_srspair(seq1, seq2, alignment))
        f.close()
        return seq1_pos, seq2_pos
    # function for the boundaries

    def align_to_get_boundaries(self, f1, f2, seq1, seq2):
        """
        get the starts of seq1 and seq2 align fragment
        """
        alignment_a_left = 0
        alignment_b_left = 0
        output_file_name = f1.rsplit('.', 1)[0] + '_' + f2.rsplit('.', 1)[0] + ".align"
        print("Saved alignment file name is: ", output_file_name)
        self.text_area.insert(tk.END, "Saved alignment file name is: %s" % output_file_name + "\n")

        alignment_a_left, alignment_b_left = self.alignment_two_seq(seq1, seq2, output_file_name)

        return alignment_a_left, alignment_b_left
    # end of function

    def merge_sanger(self):
        """
        merge DNA sequence files
        """
        self.text_area.delete('1.0', tk.END)
        self.dna_sequence.clear()
        self.merge_result = ""

        # get current path
        current_path = os.getcwd()
        print(f"Script's path is {current_path}")
        self.text_area.insert(tk.END, "Script's path is %s." % current_path + '\n')

        # make a new directory
        print(f"Log's path is {self.logdir}")
        self.text_area.insert(tk.END, "Log's path is %s." % self.logdir + '\n')

        now = datetime.now()
        formatted_time = now.strftime("%Y%m%d%H%M%S")
        folder_name = "merged_sequence" + formatted_time
        os.chdir(self.logdir)
        dirs = os.listdir(self.logdir)
        if folder_name not in dirs:
            os.mkdir(folder_name)
        else:
            shutil.rmtree(folder_name)
            os.mkdir(folder_name)

        logdir_sub = os.path.join(self.logdir, folder_name)
        os.chdir(logdir_sub)

        print("Files to be merged:")
        # print(f"{'\n'.join([tup[1] for tup in self.frames])}")
        for tup in self.frames:
            print(tup[1])
        self.text_area.insert(tk.END, "Files to be merged:" + '\n')
        self.text_area.insert(tk.END, "###########################################" + '\n')
        self.text_area.insert(tk.END, '\n'.join([tup[1] for tup in self.frames]) + '\n')
        self.text_area.insert(tk.END, "###########################################" + '\n')

        if self.check_file_list_count() == -1:
            self.text_area.insert(tk.END,
                                  "Must have 2 or more sequences. Please reinput sequences files to be merged. \n")
            return
        if self.check_radiobutton_selection() == -1:
            self.text_area.insert(tk.END,
                                  "Please choose format of sequences files (F or R) to be merged. \n")
            return

        for tup in self.frames:
            DNA_sequence_tmp_str = self.read_sequence_from_file(tup[1])
            if tup[2].get() == "F":
                self.dna_sequence.append(DNA_sequence_tmp_str)

                f = open(os.path.basename(tup[1]).rsplit('.', 1)[0] + '.newseq', 'w')
                f.write(DNA_sequence_tmp_str)
                f.close()

            if tup[2].get() == "R":
                DNA_sequence_tmp_Seq = Seq(DNA_sequence_tmp_str)
                DNA_sequence_tmp_Seq_F = DNA_sequence_tmp_Seq.reverse_complement()
                DNA_sequence_tmp_Seq_F_str = str(DNA_sequence_tmp_Seq_F)
                self.dna_sequence.append(DNA_sequence_tmp_Seq_F_str)

                f = open(os.path.basename(tup[1]).rsplit('.', 1)[0] + '.newseq', 'w')
                f.write(DNA_sequence_tmp_Seq_F_str)
                f.close()

        # align with PairwiseAligner progressively to get boundaries
        # read from self.dna_sequence list

        Es_list = []
        Es_list.append(1)

        print('Current directory is: ', os.getcwd())
        self.text_area.insert(tk.END, "Current directory is: " + os.getcwd() + "\n")

        for i in range(0, len(self.dna_sequence)-1, 1):
            file1 = self.frames[i][1]
            file2 = self.frames[i+1][1]
            direction1 = self.frames[i][2].get()
            direction2 = self.frames[i+1][2].get()
            newfile1 = os.path.basename(file1).rsplit('.', 1)[0] + '.newseq'
            newfile2 = os.path.basename(file2).rsplit('.', 1)[0] + '.newseq'
            print("Align %s:%s and %s:%s" % (newfile1, direction1, newfile2, direction2))
            self.text_area.insert(tk.END, "Align %s:%s and %s:%s \n" % (newfile1, direction1, newfile2, direction2))
            E = self.align_to_get_boundaries(os.path.basename(file1), os.path.basename(file2), self.dna_sequence[i], self.dna_sequence[i + 1])
            if E[0] == -1 and E[1] == -1:
                print(f"Merges fail found: {file1} and {file2} merges fail.\n")
                self.text_area.insert(tk.END, f"Merges fail found: {file1} and {file2} merges fail. \n\n")
            Es_list.append(E[0])
            Es_list.append(E[1])

        # length of the last sequence
        length_of_last_file_sequence = len(self.dna_sequence[-1])
        Es_list.append(length_of_last_file_sequence + 1)
        print("Boundaries list:" + str(Es_list))
        self.text_area.insert(tk.END, "Boundaries list: " + str(Es_list) + '\n\n')
        print("\n")
        # till now, got the boundaries for each sequence

        # if Es_list contains -1, report error message and stop merging
        for ii in range(len(Es_list)):
            if Es_list[ii] == -1:
                print("Merges fail found, please check pairwise alignments logs above.")
                print("You can merge sequences separately or tune the Consecutive Matches parameter.")
                self.text_area.insert(tk.END, "Merges fail found, please check pairwise alignments logs above. \n")
                self.text_area.insert(tk.END, "You can merge sequences separately or tune the Consecutive Matches parameter. \n")
                # save log file
                merged_log_name = "_merged" + ".log"
                file = open(merged_log_name, 'w')
                log_content = self.text_area.get("1.0", tk.END)
                file.write(log_content)
                file.close()

                os.chdir(current_path)
                return

        # list to store the final merged sequence
        merged_sequence_list = []

        for i in range(len(self.dna_sequence)):
            merged_sequence_list += self.dna_sequence[i][Es_list[2 * i] - 1:Es_list[2 * i + 1] - 1]
        """
        The content of Es_list, 2n elements in total, where n=len(self.dna_sequence) is the number of
        sequences to be merged.

        E0=1,E1;    E2,E3;    E4,E5;    E6,E7;    ...,    E2n-2,E2n-1
        S0          S1        S2        S3        ...,    Sn-1

        The above S0,S1,...,Sn-1 denote the sequences to be merged.
        The E0,E1,E2,...,E2n-1, which are stored in list Es_list, represent the coordinates
        of intervals to join to full-length sequence. The alignment program 
        starts nucleotide numbers from 1.
        However, in Python the list starts from 0. Thus, a mathematical interval [X,Y), which is
        equal to a mathematical interval [X,Y-1], corresponds to the list elements from
        X-1 to Y-2, i.e., list[X-1:Y-1].
        Thus, self.dna_sequence[i][Es_list[2 * i] - 1:Es_list[2 * i + 1] - 1] stands for the nucleotides
        that will be used for merging for each Sanger sequence file that has been changed to
        forward sequencing.
        """
        # print the merged sequence
        merged_sequence_str = ''.join(merged_sequence_list)
        print("The merged DNA sequence is shown below:")
        self.text_area.insert(tk.END,"The merged DNA sequence is shown below: \n")
        for j in range(0, len(merged_sequence_str), 100):
            DNA_string_100_per_line_str = merged_sequence_str[j:j + 100]
            print(DNA_string_100_per_line_str)
            self.text_area.insert(tk.END, DNA_string_100_per_line_str + "\n")
        print("\n")
        self.text_area.insert(tk.END, "\n")
        # write to the file
        merged_file_name = "merged" + ".seq"
        file = open(merged_file_name, 'w')
        file.write(merged_sequence_str)
        file.close()
        self.merge_result = merged_sequence_str

        # save log file
        merged_log_name = "_merged" + ".log"
        file = open(merged_log_name, 'w')
        log_content = self.text_area.get("1.0", tk.END)
        file.write(log_content)
        file.close()

        os.chdir(current_path)
        fullname = os.path.join(logdir_sub, merged_file_name)
        os.path.join(logdir_sub, merged_file_name)
        messagebox.showinfo("merge successfully", f"Merged Sanger file is saved as: {fullname}\n"
                                                  f"Please check the logs in the right field")
    def save_file(self):
        """
        save DNA merged sequence file, 3 formats
        """
        filetypes = (
            ("Seq files", "*.seq"),
            ("Fasta files", "*.fasta"),
            ("Text files", "*.txt")
        )

        save_path = filedialog.asksaveasfilename(
            title="Save As",
            initialfile="merge.seq",
            filetypes=filetypes,
            defaultextension=".seq"
        )

        if save_path:
            with open(save_path, 'w') as file:
                # file.write(self.text_area.get('1.0', tk.END))
                file.write(self.merge_result)
                file.close()
            messagebox.showinfo("save successfully", f"Merge Sanger file is saved as: {save_path}")


class MultiFileCLI:
    def __init__(self, txt_file):
        self.filename = txt_file
        self.results = None
        self.merge_result = ""
        self.dna_sequence = []
        self.logdir = ""

    def read_and_process_file(self):
        results = []
        try:
            with open(self.filename, 'r') as file:
                for line in file:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        seq_filename = parts[0]
                        fr_flag = parts[1]
                        results.append((seq_filename, fr_flag))
            return results
        except FileNotFoundError:
            print(f"Error: File {self.filename} Not found.")
            return None
        except Exception as e:
            print(f"Error when open file: {e}")
            return None

    def getlogdir(self):
        print(f"filename is {self.filename}")
        self.logdir = os.path.dirname(self.filename)
        if self.logdir == "":
            self.logdir = os.getcwd()
        print(f"log dir is {self.logdir}")

    def alignment_two_seq(self, seq1, seq2, save_file):
        """
        align seq1 and seq2, save alignment file
        """
        gap_open = 10.0
        gap_extend = 0.5

        # Perform Needleman-Wunsch alignment
        alignment = MultiFileGUI.perform_needleman_wunsch_pairwisealigner(seq1, seq2, gap_open, gap_extend)
        # Format alignment for output
        alignment_output = MultiFileGUI.format_alignment(seq1, seq2, alignment)

        # Print or save alignment output
        print(alignment_output)

        #get aligned seq1 and seq2
        fasta = alignment.format("fasta")
        # print(fasta)

        seq_new = fasta.split('>', 3)
        aligned_seq1 = seq_new[1].strip()
        aligned_seq2 = seq_new[2].strip()

        n = 50

        print("n=", n)

        aligned_pos = MultiFileGUI.find_n_consecutive(aligned_seq1, aligned_seq2, n)
        if aligned_pos == -1:
            print(f"can not find consecutive {n} equal in two sequences.")
            print(MultiFileGUI.format_needle_srspair(seq1, seq2, alignment))
            f = open(save_file, 'w')
            f.write(alignment_output)
            f.write(MultiFileGUI.format_needle_srspair(seq1, seq2, alignment))
            f.close()
            return -1, -1
        seq1_pos = aligned_pos - aligned_seq1[:aligned_pos].count('-') + 1 # start from 1
        seq2_pos = aligned_pos - aligned_seq2[:aligned_pos].count('-') + 1

        print(MultiFileGUI.format_needle_srspair(seq1, seq2, alignment))

        f = open(save_file, 'w')
        f.write(alignment_output)
        f.write(MultiFileGUI.format_needle_srspair(seq1, seq2, alignment))
        f.close()
        return seq1_pos, seq2_pos

    def align_to_get_boundaries(self, f1, f2, seq1, seq2):
        """
        get the starts of seq1 and seq2 align fragment
        """
        alignment_a_left = 0
        alignment_b_left = 0
        output_file_name = f1.rsplit('.', 1)[0] + '_' + f2.rsplit('.', 1)[0] + ".align"
        print("Saved alignment file name is: ", output_file_name)

        alignment_a_left, alignment_b_left = self.alignment_two_seq(seq1, seq2, output_file_name)

        return alignment_a_left, alignment_b_left

    def merge_sanger(self):
        self.getlogdir()
        results = self.read_and_process_file()
        if results is not None:
            print("Files to be merged:")
            for filename, flag in results:
                print(f"file name: {filename}, F/R flag: {flag}")
        else:
            print("read argv txt file error")
            return

        self.merge_result = ""

        # get current path
        current_path = os.getcwd()
        print(f"Script's path is {current_path}")
        print(f"Log's path is {self.logdir}")

        now = datetime.now()
        formatted_time = now.strftime("%Y%m%d%H%M%S")
        folder_name = "merged_sequence" + formatted_time
        os.chdir(self.logdir)
        dirs = os.listdir(self.logdir)
        if folder_name not in dirs:
            os.mkdir(folder_name)
        else:
            shutil.rmtree(folder_name)
            os.mkdir(folder_name)

        logdir_sub = os.path.join(self.logdir, folder_name)
        os.chdir(logdir_sub)

        if len(results) < 2:
            print("Must have 2 or more sequences. Please reinput sequences files to be merged.")
            return
        for filename, flag in results:
            if flag not in ("F", "R"):
                print(f"Error. There is a Flag {flag} not in F or R.")
                return
            os.chdir(current_path)
            DNA_sequence_tmp_str = MultiFileGUI.read_sequence_from_file(filename)
            os.chdir(logdir_sub)
            if flag == "F":
                self.dna_sequence.append(DNA_sequence_tmp_str)

                f = open(os.path.basename(filename).rsplit('.', 1)[0] + '.newseq', 'w')
                f.write(DNA_sequence_tmp_str)
                f.close()

            if flag == "R":
                DNA_sequence_tmp_Seq = Seq(DNA_sequence_tmp_str)
                DNA_sequence_tmp_Seq_F = DNA_sequence_tmp_Seq.reverse_complement()
                DNA_sequence_tmp_Seq_F_str = str(DNA_sequence_tmp_Seq_F)
                self.dna_sequence.append(DNA_sequence_tmp_Seq_F_str)

                f = open(os.path.basename(filename).rsplit('.', 1)[0] + '.newseq', 'w')
                f.write(DNA_sequence_tmp_Seq_F_str)
                f.close()

        # align with PairwiseAligner progressively to get boundaries
        # read from self.dna_sequence list
        Es_list = []
        Es_list.append(1)

        print('Current directory is: ', os.getcwd())

        for i in range(0, len(self.dna_sequence)-1, 1):
            file1 = results[i][0]
            file2 = results[i+1][0]
            direction1 = results[i][1]
            direction2 = results[i+1][1]
            newfile1 = os.path.basename(file1).rsplit('.', 1)[0] + '.newseq'
            newfile2 = os.path.basename(file2).rsplit('.', 1)[0] + '.newseq'
            print("Align %s:%s and %s:%s" % (newfile1, direction1, newfile2, direction2))
            E = self.align_to_get_boundaries(os.path.basename(file1), os.path.basename(file2), self.dna_sequence[i], self.dna_sequence[i + 1])
            if E[0] == -1 and E[1] == -1:
                print(f"Merges fail found: {file1} and {file2} merges fail.\n")
            Es_list.append(E[0])
            Es_list.append(E[1])

        # length of the last sequence
        length_of_last_file_sequence = len(self.dna_sequence[-1])
        Es_list.append(length_of_last_file_sequence + 1)
        print("Boundaries list:" + str(Es_list))

        print("\n")
        # till now, got the boundaries for each sequence
        # if Es_list contains -1, report error message and stop merging
        for ii in range(len(Es_list)):
            if Es_list[ii] == -1:
                print("Merges fail found, please check pairwise alignments logs above.")
                print("You can merge sequences separately or tune the Consecutive Matches parameter.")

                # save log file
                # merged_log_name = "_merged" + ".log"
                # file = open(merged_log_name, 'w')
                # log_content = self.text_area.get("1.0", tk.END)
                # file.write(log_content)
                # file.close()

                os.chdir(current_path)
                return

        # list to store the final merged sequence
        merged_sequence_list = []

        for i in range(len(self.dna_sequence)):
            merged_sequence_list += self.dna_sequence[i][Es_list[2 * i] - 1:Es_list[2 * i + 1] - 1]

        # print the merged sequence
        merged_sequence_str = ''.join(merged_sequence_list)
        print("The merged DNA sequence is shown below:")

        for j in range(0, len(merged_sequence_str), 100):
            DNA_string_100_per_line_str = merged_sequence_str[j:j + 100]
            print(DNA_string_100_per_line_str)

        print("\n")

        # write to the file
        merged_file_name = "merged" + ".seq"
        file = open(merged_file_name, 'w')
        file.write(merged_sequence_str)
        file.close()
        self.merge_result = merged_sequence_str

        # save log file
        # merged_log_name = "_merged" + ".log"
        # file = open(merged_log_name, 'w')
        # log_content = self.text_area.get("1.0", tk.END)
        # file.write(log_content)
        # file.close()

        os.chdir(current_path)
        fullname = os.path.join(logdir_sub, merged_file_name)

        print(f"merge successfully!\nMerged Sanger file is saved as: {fullname}\nPlease check the logs upside.")


if __name__ == "__main__":
    """
    support command:
    1.python Merge_Sanger_v3.1.py # no other argv, start gui
    2.python Merge_Sanger_v3.1.py -gui
    3.python Merge_Sanger_v3.1.py -cli test.txt 
        test.txt:
            001.seq F
            002.seq F
            003.seq F
            004.seq R
            005.seq R
            006.seq R
    4.python Merge_Sanger_v3.1.py -cli 001.seq F 002.seq F 003.seq F 004.seq R 005.seq R 006.seq R
    If the file name contains brackets, please enclose the file name in quotes, like:
    python Merge_Sanger_v3.1.py -cli "001(.seq" F "0(02).seq" F
    """
    if len(sys.argv) > 1:
        if sys.argv[1] == "-cli":
            print("CLI argv:")
            print(sys.argv[1:])
            if len(sys.argv) < 3:
                print("wrong argv, please use -cli filelist.txt format")
                sys.exit()
            elif len(sys.argv) == 3:
                filename, extension = os.path.splitext(os.path.basename(sys.argv[2]))
                if extension == ".txt":
                    print("txt")
                    cli = MultiFileCLI(sys.argv[2])
                    cli.merge_sanger()
                    sys.exit()
                else:
                    # python Merge_Sanger_v3.1.py -cli test.txt
                    print("wrong argv, please use -cli filelist.txt format")
                    sys.exit()
            elif len(sys.argv) % 2 != 0:
                # python Merge_Sanger_v3.1.py -cli 001.seq F 002.seq F 003.seq F 004.seq R 005.seq R 006.seq R
                print("wrong argv, please use -cli file1.seq F file2.seq F file3.seq F file4.seq R file5.seq R file6.seq R format")
                sys.exit()
            else:
                pairs = []
                for i in range(2, len(sys.argv), 2):
                    filename = sys.argv[i]
                    mode = sys.argv[i + 1]
                    if mode not in ['F', 'R']:
                        print(f"Error. There is a Flag {mode} not in F or R.")
                        sys.exit()
                    pairs.append((filename, mode))

                seqfilelist = "seqfilelist.txt"
                with open(seqfilelist, 'w') as file:
                    for pair in pairs:
                        file.write(f"{pair[0]} {pair[1]}\n")

                cli = MultiFileCLI(seqfilelist)
                cli.merge_sanger()
                sys.exit()

        elif sys.argv[1] == "-gui":
                print(sys.argv[1:])
                root = tk.Tk()
                app = MultiFileGUI(root)
                root.mainloop()
        else:
            print("wrong argv, please use -gui or -cli filelist.txt")
            sys.exit()
    else:
        print("no other argv")
        root = tk.Tk()
        app = MultiFileGUI(root)
        root.mainloop()
