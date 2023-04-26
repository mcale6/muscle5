#!/usr/bin/env python
import argparse
import os
import random
import string
import subprocess
import pandas as pd
from pandarallel import pandarallel
pandarallel.initialize(use_memory_fs=True, progress_bar=False, nb_workers=20)


def hhfilter_general(ifile, ofile, filter_type, value):
    hhfilter = "/home/adaddi/data/hh-suite/build/src/hhfilter"
    hhfilter_cmd = f"{hhfilter} -i {ifile} -o {ofile} -{filter_type} {value}"
    subprocess.run(hhfilter_cmd, shell=True, check=True)
    return ofile

def count_sequences_in_msa_file(file):
    count = 0
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count/2

def find_closest_cov(input_file, name, output_dir)
    cov_values = [30, 60, 90]
    closest_file = None
    closest_diff = float('inf')

    for cov_value in cov_values:
        ofile = f'{output_dir}/{name}_hhfilter_cov{str(cov_value)}_start.a3m'
        cov_output = hfilter_general(input_file, ofile, "cov", cov_value)
        seq_count = ount_sequences_in_msa_file(cov_output)
        diff = abs(seq_count - 2000)

        if diff < closest_diff:
            closest_diff = diff
            closest_file = cov_output
    return closest_file
    
def random_char(y=3):
    chain_letters = list(string.ascii_uppercase)
    return ''.join(random.choice(chain_letters) for x in range(y))

def parse_fasta(filename):
    """function to parse fasta file"""
    header, sequence = [], []
    with open(filename, "r") as lines:
        for line in lines:
            line = line.rstrip()
            # Ignore lines starting with '#'
            if line.startswith("#"):
                continue
            if line[0] == ">":
                header.append(line[1:])
            else:
                sequence.append(line)

    return header, sequence

def process_dataframe_muscle(header, seqs):
    df = pd.DataFrame({"SequenceName": header, "sequence": seqs})
    sequence_length, n_sequnce_former = len(df.sequence.iloc[0]), len(df)
    #####
    def to_fasta(s):
        return "".join(x for x in s if x.isupper())

    #####
    df["sequence"] = df.parallel_apply(lambda row: to_fasta(row["sequence"]), axis=1)
    df = pd.concat([df.drop_duplicates(subset='SequenceName'), df.drop_duplicates(subset='sequence')]).drop_duplicates()
    return df.SequenceName.to_list(), df.sequence.to_list()

def write_fasta_muscle(names, seqs, outfile):
    with open(outfile, "a") as f:
        for nm, seq in list(zip(names, seqs)):
            if nm == "101":
                nm = "QUERYSEQUENCE"
            f.write(">%s\n%s\n" % (nm, seq))
    return

def resample_muscle(output_dir, file, mode):
    save_folder = f'{output_dir}/muscleMSA_{random_char()}_{mode}'
    subprocess.run(["mkdir", save_folder], capture_output=False, text=True)
    if file != ".":
        a3m_paths = [file] 
    in_ = f'{save_folder}/muscle_unAln.fasta'
    out_ = f'{save_folder}/muscle_Aln.efa'
    out_re = f'{save_folder}/re_aln.@.afa'
    for a3m in a3m_paths:
        header, seqs = parse_fasta(a3m)
        header, seqs = process_dataframe_muscle(header, seqs)
        write_fasta_muscle(header, seqs, in_)
    
    ### Align
    print("Start MUSCLE5 MSA")
    if mode == "stratified":
        command = ["/home/adaddi/scratch/muscle_resampling/muscle5.1.linux_intel64", "-align", in_ , "-output", out_ , "-stratified"]
    if mode == "diversified":
        command = ["/home/adaddi/scratch/muscle_resampling/muscle5.1.linux_intel64", "-align", in_ , "-output", out_ , "-diversified"]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"MUSCLE encountered an error:\n{result.stderr}")
    return
    
if __name__ == "__main__":
    import argparse
    import os
    
    parser = argparse.ArgumentParser(description=os.path.basename(__file__))
    parser.add_argument('--name', type=str, help='Name')
    parser.add_argument('--output_dir', type=str, help='Name')
    parser.add_argument('--a3m_file', type=str, default="mmseqs2_crbn_uniref_env.a3m", help='A3M file name')
    parser.add_argument('--mode', type=str, default="stratified", help='Mode of sampling')
    #parser.add_argument('--value', type=int, default=600, help='Value for hhfilter')
    args = parser.parse_args()
    
    all_a3m_file = os.path.join(os.getcwd(), args.a3m_file)

    # Call the functions
    hhfilter_ofile = find_closest_cov(args.a3m_file, args.name, args.output_dir)
    muscle_foldr = resample_muscle(args.output_dir, hhfilter_ofile, args.mode)