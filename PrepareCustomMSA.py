#!/usr/bin/env python
import argparse
import os
import random
import string
import pandas as pd
from Bio import SeqIO
from pandarallel import pandarallel
pandarallel.initialize(use_memory_fs=False, progress_bar=False, nb_workers=8)

def random_char(y=3):
    return "".join(random.choice(string.ascii_letters) for x in range(y))

def load_fasta(file):
    seqs, IDs = [], []
    with open(file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq = "".join([x for x in record.seq])
            IDs.append(record.id)
            seqs.append(seq)
    return IDs, seqs


def process_input(file):
    IDs, seqs = load_fasta(file)
    seqs = [  "".join([x for x in s if x.isupper() or x == "-"]) for s in seqs ]
    df = pd.DataFrame({"SequenceName": IDs, "sequence": seqs})
    return df


def process_dataframe(df, gap_cutoff, gap_positions=[], resample=True):
    query_ = df.iloc[:1]
    query_seq_len = len(query_.sequence.values[0])
    df = df.iloc[1:]
    
    if resample:
        df = df.sample(frac=1)

    sequence_length = len(df.sequence.iloc[0])
    df["frac_gaps"] = [x.count("-") / sequence_length for x in df["sequence"]]
    df = df.loc[df.frac_gaps < gap_cutoff]
    if gap_positions != []:
        df['sequence'] = df['sequence'].parallel_apply(lambda row: add_gaps_at_positions(row, gap_positions))
    return query_, df


def write_fasta(names, seqs, outfile, xmer):
    with open(outfile, "w") as f:
        if xmer == "Monomer":
            f.write("#%s\t%s\n" % (str(len(seqs[0])), "1"))
        if xmer == "Homodimer":
            f.write("#%s\t%s\n" % (str(len(seqs[0])), "2"))
        for nm, seq in list(zip(names, seqs)):
            f.write(">%s\n%s\n" % (nm, seq))
    return

def store_pairing_info(filename, start=3):
    # if start =3 the query sequence is not included if 1 yes
    with open(filename, "r") as file:
        lines = file.readlines()

    start_line = None
    end_line = None
    for idx, line in enumerate(lines):
        if line.startswith(">101") and idx != 1:
            end_line = idx
            break
            
    extracted_lines = lines[start:end_line]
    
    seq_names = []
    sequences = []
    for idx, line in enumerate(extracted_lines):
        if line.startswith(">"):
            seq_names.append(line.strip())
            sequences.append(extracted_lines[idx+1].strip())

    return pd.DataFrame({"SequenceName": seq_names, "sequence": sequences })

def write_fasta_unpaired(names, seqs, complex_ids, outfile, cardinality, first_line_seqs):
    with open(outfile, "w") as f:
        f.write(cardinality)
        f.write('>' + "\t".join(complex_ids) + '\n')
        f.write(first_line_seqs)
        for nm, seq in list(zip(names, seqs)):
            f.write(">%s\n%s\n" % (nm, seq))
    return

def write_fasta_unpairedpaired(unpaired_info, paired_info, outfile, cardinality, first_line_seqs, complex_ids):
    with open(outfile, "w") as f:
        f.write(cardinality)
        f.write('>' + "\t".join(complex_ids) + '\n')
        f.write(first_line_seqs)
        for nm, seq in list(zip(paired_info.SequenceName.tolist(), paired_info.sequence.tolist())):
            f.write("%s\n%s\n" % (nm, seq))
        for nm, seq in list(zip(unpaired_info.SequenceName.tolist(), unpaired_info.sequence.tolist())):
            f.write(">%s\n%s\n" % (nm, seq))
    return

def sampling(query_, df, output, n_controls, fname, xmer, msa_mode, sample="NO"):
    tempsU10 = []
    for i in range(n_controls):
        if sample == "uniform":
            tmp = df.sample(n=100)
            tmp = pd.concat([query_, tmp], axis=0)
            write_fasta(tmp.SequenceName.tolist(),tmp.sequence.tolist(),
                outfile=f'{output}/{fname.split("__mmseqs2_uniref_env.a3m")[0]}_SP{str(i)}U100__{msa_mode}.a3m', xmer=xmer)
        #resample
        tmp = df.sample(n=100, replace=True, random_state=i)
        tmp = pd.concat([query_, tmp], axis=0)
        if xmer == "Heterodimer":
            tempsU10.append(tmp)
        else:
            write_fasta(tmp.SequenceName.tolist(),tmp.sequence.tolist(),
                outfile=f'{output}/{fname.split("__mmseqs2_uniref_env")[0]}_SP{str(i)}R100__{msa_mode}.a3m', xmer=xmer)
            
    tempsU100 = []   
    if len(df) > 700:
        for i in range(n_controls):
            if sample == "uniform":
                tmp = df.sample(n=700)
                tmp = pd.concat([query_, tmp], axis=0)
                write_fasta( tmp.SequenceName.tolist(), tmp.sequence.tolist(),
                    outfile=f'{output}/{fname.split("__mmseqs2_uniref_env")[0]}_SP{str(i)}U700__{msa_mode}.a3m', xmer=xmer)
            #resample
            tmp = df.sample(n=700, replace=True, random_state=i)
            tmp = pd.concat([query_, tmp], axis=0)
            if xmer == "Heterodimer":
                tempsU100.append(tmp)
            else:
                write_fasta(tmp.SequenceName.tolist(),tmp.sequence.tolist(),
                    outfile=f'{output}/{fname.split("__mmseqs2_uniref_env")[0]}_SP{str(i)}R700__{msa_mode}.a3m', xmer=xmer)
    
    tempsALL = []
    for i in range(n_controls):
        tmp = pd.concat([query_, df], axis=0)
        if xmer == "Heterodimer":
            tempsALL.append(tmp)
        else:
            write_fasta(tmp.SequenceName.tolist(),tmp.sequence.tolist(),
                outfile=f'{output}/{fname.split("__mmseqs2_uniref_env")[0]}_SP{str(i)}ALL__mmseqs2_uniref_env.a3m', xmer=xmer)

    return tempsU10, tempsU100, tempsALL

def add_gaps1(row, seq_len101, seq_len102):
    if row['id_'] == "101":
        return row['sequence'] + seq_len102 * "-"
    elif row['id_'] == "102":
        return  seq_len101 * "-" + row['sequence']
    else:
        return row['sequence']

def add_gaps_at_positions(sequence, gap_positions):
    return [''.join(['-' if i in gap_positions else c for i, c in enumerate(sequence)])][0]

def add_gap_block(row, complex_dict):
    keys = list(complex_dict.keys())
    idx = keys.index(row['id_'])
    seq_lens = [complex_dict[keys[i]][1] for i in range(len(keys)) if i != idx]
    prefix_gap = sum(seq_lens[:idx]) * "-"
    suffix_gap = sum(seq_lens[idx:]) * "-"
    return prefix_gap + row['sequence'] + suffix_gap

def generate_combined_filename(files, gaps_positions=[]):
    subnames = []
    for i, file in enumerate(files):
        gap_pos = "g"
        if  gaps_positions[i] !=[]:
            gap_pos = 'g'.join(map(str, gaps_positions[i]))
        protein = os.path.basename(file).split('__mmseqs2_uniref_env')[0].split("-")[1]
        #filter_condition = os.path.basename(file).split("-")[0]
        filter_condition = os.path.basename(file).split("-")[1]

        subnames.append(f"{protein}_{gap_pos}_{filter_condition}_{random_char(2)}")

    return "_CW_".join(subnames) + '__mmseqs2_uniref_env'

def execute_sampling(parameters):
    output = f"{parameters['save_location']}" 
    os.makedirs(output, exist_ok=True)

    if parameters['pairing'] != ".":
        paired_info = store_pairing_info(parameters['pairing'])

    if parameters['xmer'] == "Monomer" or "Homodimer":
        df = process_input(parameters['file'][0])
        query_, df = process_dataframe(df, parameters['gap_cutoff'], parameters['gaps_position'], parameters['resample'])
        sampling(query_, df, output, parameters['n_controls'], parameters['fname'], 
            parameters['xmer'], parameters['sample'])
    if parameters['xmer'] == "Heterodimer":
        tempsU10_, tempsU100_, tempsALL_ = [], [], []
        complex_dict = {i: (None, None) for i in parameters['complex_ids']}
        for idx, i in enumerate(parameters['complex_ids']):
            df = process_input(parameters['file'][idx])
            query_, df = process_dataframe(df, parameters['gap_cutoff'], parameters['gaps_position'][idx], parameters['resample'])
            #
            query_.SequenceName.iloc[0] = i
            query_["id_"] = i
            complex_dict[i] = (query_.sequence.iloc[0], len(query_.sequence.iloc[0]))
            df["id_"] = [i] * len(df.index)
            #
            tempsU10, tempsU100, tempsALL = sampling(query_, df, output, parameters['n_controls'], parameters['fname'],
            parameters['xmer'], parameters['sample'])
            tempsU10_.append(tempsU10), tempsU100_.append(tempsU100), tempsALL_.append(tempsALL)

        seq_lens = [v[1] for k, v in complex_dict.items()]
        cardinality = "#" + ",".join(map(str, seq_lens)) + "\t" + ",".join(["1"] * len(complex_dict)) + "\n"
        first_line_seqs = "".join([v[0] for k, v in complex_dict.items()]) + "\n"
        #print(seq_lens, cardinality, first_line_seqs)
        for u_idx, t in enumerate([tempsU10_, tempsU100_, tempsALL_]):
            if u_idx == 0:
                u_idx =10
            if u_idx == 1:
                u_idx = 100
            if u_idx == 2:
                u_idx = "ALL"
            for i in range(parameters["n_controls"]):
                tmp = 0
                try:
                    tmp = pd.concat([t[j][i] for j in range(len(t))], axis=0, ignore_index=True)
                    #print("ok")
                except:
                    print("Error")
                    print(parameters['fname'], u_idx)
                    print([len(t) for t in tempsU10_])
                    continue
                tmp["sequence"] = tmp.parallel_apply(lambda row: add_gap_block(row, complex_dict), axis=1)
                if parameters['pairing'] == ".":
                    write_fasta_unpaired(tmp.SequenceName.tolist(), tmp.sequence.tolist(), parameters['complex_ids'],
                    f"{output}/{parameters['fname'].split('__mmseqs2_uniref_env')[0]}_SP{str(i)}R{u_idx}__{parameters['msa_mode']}.a3m",
                    cardinality, first_line_seqs)
                if parameters['pairing'] != ".":
                    write_fasta_unpairedpaired(tmp, paired_info,
                    f"{output}/{parameters['fname'].split('__mmseqs2_uniref_env')[0]}_SP{str(i)}R{u_idx}__{parameters['msa_mode']}.a3m",
                    cardinality, first_line_seqs,  parameters['complex_ids'])

    return

if __name__ == "__main__":
    #parser = argparse.Argument#Parser(description='New MSA from homologous protein sequences using MUSCLE')
    #parser.add_argument('--out_dir1', type=str)
    #parser.add_argument('--out_dir2', type=str)
    #parser.add_argument('--out_dir', type=str)
    #parser.add_argument('--xmer', type=str)
    #args = #parser.parse_args()
    #files1 = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(args.out_dir1)) for f in fn if f[-3:] == "a3m" and 'mmseqs2' in os.path.join(dp, f)]
    #files2 = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(args.out_dir2)) for f in fn if f[-3:] == "a3m" and 'mmseqs2' in os.path.join(dp, f)]
    a3m_dir1 = "/home/adaddi/scratch/muscle_resampling/Sampling_super5_tcfMETTL3__mmseqs2_uniref_env_UIC/muscleMSA__IVZ-super5_startcov_15-tcfMETTL3__mmseqs2_uniref_env"
    a3m_dir2 = "/home/adaddi/scratch/muscle_resampling/Sampling_super5_tcfMETTL14__mmseqs2_uniref_env_KOP/muscleMSA__GZY-super5_startcov_15-tcfMETTL14__mmseqs2_uniref_env"
    a3m_dir3 = "/home/adaddi/scratch/muscle_resampling/Sampling_super5_tcfCRBN__mmseqs2_uniref_env_GKD/muscleMSA__DRY-super5_startcov_15-tcfCRBN__mmseqs2_uniref_env"
    out_dir =  "/home/adaddi/scratch/muscle_resampling/CRN_METTL3_METTL14"
    pairing_info = "/home/adaddi/scratch/muscle_resampling/tcf_m3_m14_crbn_unpaired_paired.a3m"

    xmer = "Multimer"

    files1 = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(a3m_dir1)) for f in fn if f[-3:] == "a3m" and 'mmseqs2' in os.path.join(dp, f)][:5]
    files2 = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(a3m_dir2)) for f in fn if f[-3:] == "a3m" and 'mmseqs2' in os.path.join(dp, f)][:5]
    files3 = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(a3m_dir3)) for f in fn if f[-3:] == "a3m" and 'mmseqs2' in os.path.join(dp, f)][:5]

    for file1, file2, file3 in zip(files1, files2, files3):
            params = {'file': [file1 , file2, file3], #, file1], ########
            'fname': "f_crbn_m3_m14",
            'gap_cutoff': 0.25,
            'resample': True,
            'n_controls': 3,
            'xmer': xmer,
            'complex_ids': ['101', '102', '103'], #######
            'sample': "resample",
            'pairing': pairing_info,
            'msa_mode': "mmseqs2_uniref_env_unpaired_paired", ####
            'save_location': out_dir,
            'gaps_position': [[],[], []]} ,#,[list(range(10,18)), list(range(276,287))]} # [[],[]]} ######

            params["fname"] = generate_combined_filename(params["file"], params["gaps_position"])
            print("x", params["fname"] )
            output = execute_sampling(params)