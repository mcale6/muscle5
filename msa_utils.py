import pandas as pd
import numpy as np
from pandarallel import pandarallel
import string
import random
from Bio import SeqIO
from Bio import AlignIO
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
pandarallel.initialize(use_memory_fs=False, progress_bar=False, nb_workers=12)
import itertools
from concurrent.futures import ProcessPoolExecutor

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

def store_pairing_info(filename):
    with open(filename, "r") as file:
        lines = file.readlines()

    start_line = None
    end_line = None
    for line in extracted_lines:
        if line.startswith(">"):
            seq_name = line.strip()
            sequence = ""
        else:
            sequence += line.strip()
            data.append({"SequenceName": seq_name, "sequence": sequence})

    df = pd.DataFrame(data)
    
    return df

def process_dataframe_plot(header, seqs, gap_cutoff=-1000):
    alphabet_dict = {
        char: idx for idx, char in enumerate("ACDEFGHIKLMNPQRSTVWYX-")
    }  # to do everyhing not defined is X
    df = pd.DataFrame({"SequenceName": header, "sequence": seqs})
    sequence_length, n_sequnce_former = len(df.sequence.iloc[0]), len(df)
    #####
    def reformat(s):
        return "".join(x for x in s if x.isupper() or x == "-")

    def str2int(s):
        return np.asarray([alphabet_dict[char] for char in s])

    def frac_gaps(s):
        return (np.count_nonzero(np.asarray(s) == 21) / sequence_length) < gap_cutoff

    #####
    df["sequence"] = df.parallel_apply(lambda row: reformat(row["sequence"]), axis=1)
    msa = df.parallel_apply(lambda row: str2int(row["sequence"]), axis=1)
    idxtofilter = msa.parallel_apply(lambda row: frac_gaps(row))
    return pd.DataFrame.from_records(pd.DataFrame(msa.loc[idxtofilter])[0].to_list())


def get_start_end_indices(msa, query_seq_id="QUERYSEQUENCE"):
    query_seq = None
    
    for record in msa:
        if record.id == query_seq_id:
            query_seq = record
            print(query_seq)
            break

    if query_seq is None:
        raise ValueError(f"Query sequence with ID '{query_seq_id}' not found in the MSA")

    start = next((i for i, c in enumerate(query_seq.seq) if c != '-'), None)
    end = next((i for i, c in enumerate(reversed(query_seq.seq)) if c != '-'), None)
    if end is not None:
        end = len(query_seq.seq) - end
    return start, end

def bio_msa_to_numpy(msa):
    return np.array([list(rec.seq) for rec in msa], np.dtype("U1"))

def calculate_letter_confidence(msas, query_seq_id):
    num_msas = len(msas)
    reference_msa = msas[0]  # Assume the first MSA is the reference
    ref_query_seq = [seq_record for seq_record in reference_msa if seq_record.id == query_seq_id][0]
    letter_confidence = []

    for pos in range(len(ref_query_seq.seq)):
        x = ref_query_seq.seq[pos]
        if x == "-":  # Skip gaps
            letter_confidence.append(0)
            continue

        well_aligned_count = 0
        for msa in msas[1:]:  # Exclude the reference MSA
            query_seq = [seq_record for seq_record in msa if seq_record.id == query_seq_id][0]
            column_idx = query_seq.seq.index(x)

            column = tuple(record.seq[column_idx] for record in msa)
            majority_letter = max(set(column), key=column.count)

            if majority_letter == x:
                well_aligned_count += 1

            if num_msas > 1:
                letter_confidence.append(well_aligned_count / (num_msas - 1))
            else:
                print(tuple(record.seq[column_idx] for record in msa))
                letter_confidence.append(0)

    return letter_confidence

def calculate_column_confidence_gaps(msas):
    column_confidence = {}
    num_msas = len(msas)

    for msa in msas:
        for col_idx in range(len(msa[0].seq)):
            column = tuple(record.seq[col_idx] for record in msa if col_idx < len(record.seq))
            if col_idx not in column_confidence:
                column_confidence[col_idx] = 0
            column_confidence[col_idx] += (column.count("-") / len(column)) ** 2

    # Divide the summed values by the total number of MSAs
    for col_idx in column_confidence:
        column_confidence[col_idx] /= num_msas

    return column_confidence

def hamming_similarity(s1, s2):
    return sum(ch1 == ch2 for ch1, ch2 in zip(s1, s2)) / len(s1)

def calc_avg_similarity(column):
    pairwise_similarities = [hamming_similarity(s1, s2) for s1, s2 in itertools.combinations(column, 2)]
    return sum(pairwise_similarities) / len(pairwise_similarities)

def calculate_column_similarity(msas):
    num_columns = len(msas[0][0].seq)
    column_similarity = np.zeros(num_columns)

    # Precompute column representations for each MSA
    precomputed_columns = [
        [np.array([record.seq[i] for record in msa if i < len(record.seq)]) for i in range(num_columns)]
        for msa in msas
    ]

    for msa_columns in precomputed_columns:
        for i, column_i in enumerate(msa_columns):
            if column_i.size > 0:
                for column_j in precomputed_columns:
                    similarity = np.sum(column_i != column_j) / column_i.size
                    column_similarity[i] += similarity

    column_similarity /= len(msas)

    return column_similarity

def calculate_column_confidence(msas):
    column_similarity = calculate_column_similarity(msas)
    column_confidence = {}

    for col_idx in range(column_similarity.shape[0]):
        column_similarity_sum = np.sum(column_similarity[col_idx])
        column_confidence[col_idx] = column_similarity_sum / (len(msas) - 1)

    return column_confidence

def get_query_sequence(msa, query_seq_id="QUERYSEQUENCE"):
    query_seq = None
    for record in msa:
        if record.id == query_seq_id:
            query_seq = str(record.seq)
            break

    return "".join([c for c in query_seq if c.isupper()])

def column_confidence_stats(input_fastas, query_seq_id, threshold):
    msas = [AlignIO.read(input_fasta, "fasta") for input_fasta in input_fastas]
    msas_as_lists = [list(SeqIO.parse(input_fasta, "fasta")) for input_fasta in input_fastas]
    column_confidence = [calculate_column_confidence([msa]) for msa in msas_as_lists]
    #letter_confidence = [calculate_letter_confidence([msa], query_seq_id) for msa in msas_as_lists]

    data = []
    for msa_index, (msa, cc, lc) in enumerate(zip(msas, column_confidence, letter_confidence)):
        below_threshold_cc_count = sum(1 for x in cc if x < threshold)
        average_cc = sum(cc) / len(cc)
        below_threshold_lc_count = sum(1 for x in lc if x < threshold)
        average_lc = sum(lc) / len(lc)
        #nf = calculate_nf(msa)
        #gap_frac = gap_fraction(msa)
        
        data.append({
            'MSA_Index': input_fastas,
            'Below_Threshold_CC_Count': below_threshold_cc_count,
            'Average_CC': average_cc,
            #'Below_Threshold_LC_Count': below_threshold_lc_count,
            #'Average_LC': average_lc,
            #'Nf': nf,
            #'Gap_Fraction': gap_frac
        })

    return pd.DataFrame(data)

def calculate_nf(msa, eff_cutoff=0.8):
    if msa.ndim == 3:
        msa = msa.argmax(-1)
    # pairwise identity
    msa_sm = 1.0 - squareform(pdist(msa, "hamming"))
    msa_sm = (msa_sm >= eff_cutoff).astype(np.float64)
    # Calculate the outer sum
    nf_sum = np.sum(1 / (1 + np.sum(msa_sm, -1)))
    # Calculate Nf
    nf = nf_sum / np.sqrt(msa_sm.shape[1])
    return nf

def gap_fraction(msa):
    gap_count = np.sum(msa == '-')
    total_positions = msa.shape[0] * msa.shape[1]
    return gap_count / total_positions

def plot_msa_v2(feature_dict, sort_lines=True, dpi=100):
    seq = feature_dict["msa"][0]
    if "asym_id" in feature_dict:
        Ls = [0]
        k = feature_dict["asym_id"][0]
        for i in feature_dict["asym_id"]:
            if i == k:
                Ls[-1] += 1
            else:
                Ls.append(1)
            k = i
    else:
        Ls = [len(seq)]
    Ln = np.cumsum([0] + Ls)

    try:
        N = feature_dict["num_alignments"][0]
    except:
        N = feature_dict["num_alignments"]

    msa = feature_dict["msa"][:N]
    gap = msa != 21
    qid = msa == seq
    gapid = np.stack([gap[:, Ln[i] : Ln[i + 1]].max(-1) for i in range(len(Ls))], -1)
    lines = []
    Nn = []
    for g in np.unique(gapid, axis=0):
        i = np.where((gapid == g).all(axis=-1))
        qid_ = qid[i]
        gap_ = gap[i]
        seqid = np.stack(
            [qid_[:, Ln[i] : Ln[i + 1]].mean(-1) for i in range(len(Ls))], -1
        ).sum(-1) / (g.sum(-1) + 1e-8)
        non_gaps = gap_.astype(float)
        non_gaps[non_gaps == 0] = np.nan
        if sort_lines:
            lines_ = non_gaps[seqid.argsort()] * seqid[seqid.argsort(), None]
        else:
            lines_ = non_gaps[::-1] * seqid[::-1, None]
        Nn.append(len(lines_))
        lines.append(lines_)

    Nn = np.cumsum(np.append(0, Nn))
    lines = np.concatenate(lines, 0)
    plt.figure(figsize=(8, 5), dpi=dpi)
    plt.title(f"Sequence coverage, Nf={feature_dict['Nf']}")
    plt.imshow(
        lines,
        interpolation="nearest",
        aspect="auto",
        cmap="rainbow_r",
        vmin=0,
        vmax=1,
        origin="lower",
        extent=(0, lines.shape[1], 0, lines.shape[0]),
    )
    for i in Ln[1:-1]:
        plt.plot([i, i], [0, lines.shape[0]], color="black")
    for j in Nn[1:-1]:
        plt.plot([0, lines.shape[1]], [j, j], color="black")

    plt.plot((np.isnan(lines) == False).sum(0), color="black")
    plt.xlim(0, lines.shape[1])
    plt.ylim(0, lines.shape[0])
    plt.colorbar(label="Sequence identity to query")
    plt.xlabel("Positions")
    plt.ylabel("Sequences")
    plt.show()
    return

def do_featdic(m):
    return {"msa": m,"num_alignments":  m.shape[0], "Nf": round(calculate_nf(m),2)}

def plot_msa(msa, n=50, N=10, replace=False, msa_w=None, sample="False"):
    feature_dict = do_featdic(msa.to_numpy())
    plot = plot_msa_v2(feature_dict)
    res = []
    if sample == "True":
        for i in range(N):
            tmp = msa.sample(n=n, replace=replace, random_state=i, weights=msa_w)
            m = msa.to_numpy()[tmp.index.values]
            feature_dict = do_featdic(m)
            plot = plot_msa_v2(feature_dict)
            #res.append(feature_dict)
        return res

def insepct_msa(file, sample="False"):
    pandarallel.initialize(use_memory_fs=False, progress_bar=False, nb_workers=12)
    header, seqs = parse_fasta(file)
    msa = process_dataframe_plot(header, seqs, 1)
    nf = calculate_nf(msa.to_numpy())
    plot_msa(msa)
