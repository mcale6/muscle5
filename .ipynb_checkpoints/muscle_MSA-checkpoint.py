#!/usr/bin/env python
import argparse
import os
import random
import string
import subprocess
import pandas as pd
from pandarallel import pandarallel
pandarallel.initialize(use_memory_fs=True, progress_bar=False, nb_workers=20)


def hhfilter_qsc(ifile, ofile, qsc):
    #qsc increase ==  for more distant homologs
    hhfilter = "/home/adaddi/data/hh-suite/build/src/hhfilter"
    hhfilter_cmd = f"{hhfilter} -i {ifile} -o {ofile} -qsc {qsc/100}"
    subprocess.run(hhfilter_cmd, shell=True, check=True)
    return ofile

def hhfilter_general(ifile, ofile, mode, value):
    #qsc increase ==  for more distant homologs
    hhfilter = "/home/adaddi/data/hh-suite/build/src/hhfilter"
    hhfilter_cmd = f"{hhfilter} -i {ifile} -o {ofile} -{mode} {value}"
    subprocess.run(hhfilter_cmd, shell=True, check=True)
    return ofile
    
def random_char(y=3):
    chain_letters = list(string.ascii_uppercase)
    return ''.join(random.choice(chain_letters) for x in range(y))

def parse_fasta(filename, cardinality):
    """function to parse fasta file"""
    header, sequence = [], []
    with open(filename, "r") as lines:
        for line in lines:
            line = line.rstrip()
            # Ignore lines starting with '#'
            if line.startswith("#") or line.startswith(cardinality):
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

def resample_muscle(file, mode, value, cardinality):
    wrk_dir = os.getcwd()
    save_folder = f'{wrk_dir}/muscle_{random_char()}_{mode}{value}_res'
    subprocess.run(["mkdir", save_folder], capture_output=False, text=True)
    if file != ".":
        a3m_paths = [file] 
    in_ = f'{save_folder}/muscle_unAln.fasta'
    out_ = f'{save_folder}/muscle_Aln.efa'
    out_re = f'{save_folder}/re_aln.@.afa'
    for a3m in a3m_paths:
        header, seqs = parse_fasta(a3m, cardinality)
        header, seqs = process_dataframe_muscle(header, seqs)
        header, seqs = process_dataframe_muscle(header, seqs) ### deduplicate and unaligned sequnces
        write_fasta_muscle(header, seqs, in_)
    ### Align
    print("Start MUSCLE5 MSA")
    command = ["/home/adaddi/scratch/muscle_resampling/muscle5.1.linux_intel64", "-align", in_ , "-output", out_ , "-stratified"]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"MUSCLE encountered an error:\n{result.stderr}")
    return
    

if __name__ == "__main__":
    import argparse
    import os
    
    parser = argparse.ArgumentParser(description=os.path.basename(__file__))
    parser.add_argument('--a3m_file', type=str, default="mmseqs2_crbn_uniref_env.a3m", help='A3M file name')
    parser.add_argument('--cardinality', type=str, default="442", help='Cardinality value')
    parser.add_argument('--mode', type=str, default="diff", help='Mode of sampling')
    parser.add_argument('--value', type=int, default=600, help='Value for hhfilter')
    args = parser.parse_args()
    
    all_a3m_file = os.path.join(os.getcwd(), args.a3m_file)
    hhfilter_ofile = f'{all_a3m_file.split(".")[0]}_hhfilter_{args.mode}{str(args.value)}.a3m'

    # Call the functions
    hhfilter_general(all_a3m_file, hhfilter_ofile, args.mode, args.value)
    muscle_foldr = resample_muscle(hhfilter_ofile, args.mode, args.value, args.cardinality)