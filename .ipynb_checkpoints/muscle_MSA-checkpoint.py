#!/usr/bin/env python
import argparse
import os
import random
import string
import subprocess
import pandas as pd
from pandarallel import pandarallel
import sys
sys.path.append("/home/adaddi/data/muscle5/")
from MSAFilter import *
pandarallel.initialize(use_memory_fs=True, progress_bar=False, nb_workers=20)
   
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

def MSA_sample_muscle(file, mode):
    file_name = os.path.basename(file)
    save_folder = f'{OUTPUT_DIR}/muscleMSA__{random_char()}-{mode}_{file_name.split(".")[0]}'
    subprocess.run(["mkdir", save_folder], capture_output=False, text=True)
    if file != ".":
        a3m_paths = [file] 
    in_ = f'{save_folder}/unAln-{mode}_{file_name.split(".")[0]}.fasta'
    out_ = f'{save_folder}/Aln-{mode}_{file_name.split(".")[0]}.efa'
    #out_re = f'{save_folder}/re_aln.@.afa'
    for a3m in a3m_paths:
        header, seqs = parse_fasta(a3m)
        header, seqs = process_dataframe_muscle(header, seqs)
        write_fasta_muscle(header, seqs, in_)
    
    ### Align
    print("Start MUSCLE5 MSA")
    if mode == "stratified":
        command = [MUSCLE, "-align", in_ , "-output", out_ , "-stratified"]
    if mode == "diversified":
        command = [MUSCLE, "-align", in_ , "-output", out_ , "-diversified"]
    if mode == "super5":
        out_ = out_[:-4] + ".@.afa"
        for i in [1, 3, 10, 20]:
            command = [MUSCLE, "-super5", in_ , "-output", out_ , "-perturb", str(i), "-perm", "all"]
            result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"MUSCLE encountered an error:\n{result.stderr}")
    return
    
if __name__ == "__main__":
    # python /home/adaddi/data/muscle5/muscle_MSA.py /home/adaddi/scratch/muscle_resampling/c_crbn_mmseqs2_uniref_env_org.a3m c_crbn_uniref_env diversified
    parser = argparse.ArgumentParser(description='New MSA from homologous protein sequences using MUSCLE')
    parser.add_argument('--a3m_file', type=str, help='A3M file name')
    parser.add_argument('--mode', type=str, default="stratified", help='Mode of sampling')
    args = parser.parse_args()
    ####
    MUSCLE = "/home/adaddi/data/muscle5/src/Linux/muscle"
    ###
    main_dir = "/home/adaddi/scratch/muscle_resampling"
    NAME = os.path.basename(args.a3m_file)
    OUTPUT_DIR = f'{main_dir}/Sampling_{args.mode}_{NAME.split(".")[0]}'
    create_filter_dir(OUTPUT_DIR)
    ### Filtering
    hhfilter_ofile = execute_filter(args.a3m_file, OUTPUT_DIR, NAME)
    ### Sampling MSA
    muscle_foldr = MSA_sample_muscle(hhfilter_ofile, args.mode)