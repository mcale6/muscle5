import os
import subprocess

MUSCLE = "/home/adaddi/data/muscle5/src/Linux/muscle"
HHREFROMAT = "/home/adaddi/data/hh-suite/scripts/reformat.pl"

def reformat_musclefasta_a3m(muscle_foldr, name):
    files = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(muscle_foldr)) for f in fn if f[-3:] == "afa"]
    for idx, in_ in enumerate(files):
        out_ = f'{foldr}/B{idx+1}-{in_.split("-")[-1]}
        re_file = f'{out_.split(".")[0]}_re.fasta'
        ####
        print("Reorder")
        reorder_a3m_file(in_, re_file)
        print("MSA FASTA TO A3M")
        reformat_cmd = f'{HHREFROMAT} fas a3m {re_file} {out_}'
        subprocess.run(reformat_cmd, shell=True, check=True)
        print("Unblock")
        unblock_a3m(out_)
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

def reorder_a3m_file(input_file, re_file):
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
            query_sequence = line[i+1]
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

def resample_muscle(in_):
    out_re =  f'{in_.split(".")[0]}___.@.afa'
    command = [MUSCLE, "-resample", in_ , "-output", out_re, "-replicates", "5", "-minconf",  "1", "-max_gap_fract", "1"]  
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"MUSCLE encountered an error:\n{result.stderr}")
    return

def resampled_MSA(foldr):
    in_ = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(foldr)) for f in fn if f[-3:] == "efa"]
    resample_muscle(in_[0])
    reformat_musclefasta_a3m(foldr)
    res = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(foldr)) for f in fn if f[-3:] == "a3m"]
    return res

#os.makedirs(output_foldr, exist_ok=True)
#for file in to_mv:
#    shutil.move(file, os.path.join(output_foldr, os.path.basename(file)))
#return [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(output_foldr)) for f in fn if f[-3:] == "a3m"]