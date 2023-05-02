from Bio import SeqIO
import random
import sys, os
from Bio.Seq import Seq
from ResamplingMSA import reorder_afa_file, msa_fasta_to_a3m, unblock_a3m
from msa_utils import*

def sample_columns(input_fastas, output_fasta, query_seq_id, option='uniform', min_confidence=0.5):
    # Read all MSAs into memory
    msas = [list(SeqIO.parse(input_fasta, "fasta")) for input_fasta in input_fastas]
    
    # Get query sequence
    query_sequence = get_query_sequence(msas[0], query_seq_id)
    
    # Get the start and end indices for each MSA
    indices = [(0, len(msa[0].seq)) for msa in msas]
    
    # Compute the maximum number of columns to sample
    num_columns = min(end - start for start, end in indices)

    # Calculate column confidence (CC) or letter confidence (LC) based on the specified option
    if option == 'cc':
        confidence = calculate_column_confidence_gaps(msas)
    elif option == 'lc':
        confidence = calculate_letter_confidence(msas, query_seq_id)
    elif option == 'uniform':
        confidence = None
    else:
        raise ValueError("Invalid option, must be 'cc', 'lc', or 'uniform'")

    # Initialize the resampled MSA
    resampled_msa = []

    # Iterate through each column
    for col_idx in range(num_columns):
        # If CC or LC is specified, find the MSA with the highest confidence value for the current column
        if confidence is not None:
            max_confidence = -1
            chosen_msa_index = -1

            for msa_index, (msa, (start, end)) in enumerate(zip(msas, indices)):
                if col_idx >= start and col_idx < end:
                    conf = confidence[col_idx]
                    if conf > max_confidence:
                        max_confidence = conf
                        chosen_msa_index = msa_index

            # If the maximum confidence value is below the minimum threshold,
            # select the highest confidence column for each MSA
            if max_confidence < min_confidence:
                chosen_msa_index = random.randint(0, len(msas) - 1)
        else:
            # If no option specified, randomly choose an MSA
            chosen_msa_index = random.randint(0, len(msas) - 1)

        # Select the chosen MSA
        msa = msas[chosen_msa_index]

        # Add the selected column to the resampled MSA
        for record_index, record in enumerate(msa):
            if len(resampled_msa) <= record_index:
                resampled_msa.append(SeqIO.SeqRecord(Seq(''), id=record.id, description=''))

            resampled_msa[record_index].seq += record.seq[col_idx]

    # Write the resampled MSA to the output file
    SeqIO.write(resampled_msa, output_fasta, "fasta")

    # Reorder afa_file
    reorder_sampled_afa_file(output_fasta, query_sequence, query_seq_id)
    return

def reorder_sampled_afa_file(input_file, query_sequence, query_seq_id=">QUERYSEQUENCE"):
    print("REORDER QUERY ON TOP")
    with open(input_file, "r") as file:
        lines = file.readlines()

    query_lines = []
    other_lines = []
    sampled_query_sequence = []

    i = 0
    while i < len(lines):
        line = lines[i]
        if line.strip() == (">"+ query_seq_id):
            print("Has Querysequence")
            i += 1  # skip the line with ">QUERYSEQUENCE"
            while i < len(lines) and not lines[i].startswith(">"):
                sampled_query_sequence.append(lines[i].strip())
                i += 1
        else:
            other_lines.append(line)
            i += 1
            sequence = []
            while i < len(lines) and not lines[i].startswith(">"):
                sequence.append(lines[i].strip())
                i += 1
            other_lines.append("".join(sequence) + "\n")
    #+[f">101\n{''.join(query_sequence)}\n"]+
    reordered_lines = [f'#{str(len(query_sequence))}\t1\n'] + [f">101\n{''.join(sampled_query_sequence)}\n"] + other_lines

    with open(input_file, "w") as file:
        file.writelines(reordered_lines)

    return

if __name__ == "__main__":
    foldr = sys.argv[1]
    input_fastas = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(foldr)) for f in fn if f[-3:] == "afa"]

    output_fasta = sys.argv[-1]

    sample_columns(input_fastas, output_fasta, option='uniform', min_confidence=0.5, query_seq_id="QUERYSEQUENCE")

