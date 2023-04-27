#!/usr/bin/env python
import os
import subprocess
import argparse

def reformat_musclefasta_a3m(muscle_foldr):
    files = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(muscle_foldr)) for f in fn if f[-3:] == "afa"]
    for idx, in_ in enumerate(files):
        out_ = f'{muscle_foldr}/B{idx+1}-{in_.split("-")[-1].split(".")[0]}.a3m'
        re_file = f'{out_.split(".")[0]}_re.fasta'
        ####
        print("Reorder")
        reorder_afa_file(in_, re_file)
        print("MSA FASTA TO A3M")
        reformat_cmd = f'{HHREFROMAT} fas a3m {re_file} {out_}'
        subprocess.run(reformat_cmd, shell=True, check=True)
        #print("Unblock")
        #unblock_a3m(out_)
    return

def unblock_a3m(un_file):
    with open(un_file, "r") as file:
        lines = file.readlines()

    header_lines = lines[:1]
    output_lines = []

    i = 1
    while i < len(lines):
        line = lines[i]
        if line.startswith(">"):
            output_lines.append(line)
            i += 1
            sequence = []
            while i < len(lines) and not lines[i].startswith(">"):
                sequence.append(lines[i].strip())
                i += 1
            output_lines.append("".join(sequence) + "\n")
        else:
            i += 1

    all_lines = header_lines + output_lines
    
    with open(un_file, "w") as file:
        file.writelines(all_lines)
        
    return

def unblock_efa(input_file):
    print("Unblock efa")
    ub_in = input_file.replace("Aln-", "UBAln-")
    # Open input file for reading
    with open(input_file, 'r') as f:
        file_contents = f.read()

    # Reformat the contents of the input file
    lines = file_contents.split('\n')
    reformatted_lines = []

    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith('>') or line.startswith('<'):
            if i > 0:
                reformatted_lines.append('\n' + line)
            else:
                reformatted_lines.append(line)
            i += 1
        else:
            temp = []
            while i < len(lines) and not lines[i].startswith('>') and not lines[i].startswith('<'):
                temp.append(lines[i].strip())
                i += 1
            reformatted_lines.append('\n' + "".join(temp))

    reformatted_text = ''.join(reformatted_lines)

    # Write the reformatted contents to the output file
    with open(ub_in, 'w') as f:
        f.write(reformatted_text)

    return ub_in

def reorder_afa_file(input_file, re_file):
    print("REORDER QUERY ON TOP")
    with open(input_file, "r") as file:
        lines = file.readlines()

    query_lines = []
    other_lines = []
    query_sequence = ""

    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith(">QUERYSEQUENCE"):
            print("Has Querysequence")
            query_sequence = lines[i+1]
            #query_lines.append(f">101\n{query_sequence}")  # when bootstrapping the query sequence (columns) will change
            i += 2  # skip the aln and the next line
        else:
            other_lines.append(line)
            i += 1
            sequence = []
            while i < len(lines) and not lines[i].startswith(">"):
                sequence.append(lines[i].strip())
                i += 1
            other_lines.append("".join(sequence) + "\n")

    reordered_lines = [f'#{str(len(query_sequence))}\t1\n'] + [f">101\n{query_sequence}\n"] + other_lines
    
    with open(re_file, "w") as file:
        file.writelines(reordered_lines)

    return

def resample_muscle(efa):
    in_ = unblock_efa(efa)
    out_re =  f'{in_.split(".")[0]}.@.afa'
    command = [MUSCLE, "-resample", in_ , "-output", out_re, "-replicates", "5", "-minconf",  "0.5", "-max_gap_fract", "0.5"]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"MUSCLE encountered an error:\n{result.stderr}")
    return

def execute_resampling(foldr):
    efa = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(foldr)) for f in fn if f[-3:] == "efa"]
    resample_muscle(efa[0])
    reformat_musclefasta_a3m(foldr)
    res = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(foldr)) for f in fn if f[-3:] == "a3m"]
    return res

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Resample MSA')
    parser.add_argument('--foldr', type=str, help='A3M file name')
    args = parser.parse_args()
    MUSCLE = "/home/adaddi/data/muscle5/src/Linux/muscle"
    HHREFROMAT = "/home/adaddi/data/hh-suite/scripts/reformat.pl"
    #   
    res = execute_resampling(args.foldr)


#os.makedirs(output_foldr, exist_ok=True)
#for file in to_mv:
#    shutil.move(file, os.path.join(output_foldr, os.path.basename(file)))
#return [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(output_foldr)) for f in fn if f[-3:] == "a3m"]