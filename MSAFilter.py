import os
import subprocess

def generate_output_file_name(input_file, filter_name, filter_value, output_dir, name):
    return f"{output_dir }/{filter_name}_{str(filter_value)}-{name}"

def create_filter_dir(output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    return

def hhfilter_general(ifile, ofile, mode, value):
    hhfilter = "/home/adaddi/data/hh-suite/build/src/hhfilter"
    hhfilter_cmd = f"{hhfilter} -i {ifile} -o {ofile} -{mode} {value}"
    subprocess.run(hhfilter_cmd, shell=True, check=True)
    return ofile

def count_sequences_in_msa_file(file):
    count = 0
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count/2

def find_closest_cov(input_file,  output_dir, name):
    print("HHFilter find_closest_cov")
    cov_values = [15, 20, 25, 30]
    closest_file = None
    closest_diff = float('inf')

    for cov_value in cov_values:
        ofile = generate_output_file_name(input_file,  "startcov", cov_value,  output_dir, name)
        hhfilter_general(input_file, ofile, "cov", cov_value)
        seq_count = count_sequences_in_msa_file(ofile)
        diff = abs(seq_count - 2000)

        if diff < closest_diff:
            closest_diff = diff
            closest_file = ofile
    print("HHFilter main: ", closest_file, seq_count)
    return closest_file

def filter_msa_A(closest_cov_file,  output_dir, name):
    print("HHFilter A")
    qsc_values = [0.25, 0.5, 0.7]
    qsc_files = []
    
    for qsc_value in qsc_values:
        qsc_output = generate_output_file_name(closest_cov_file, "filterA", qsc_value,  output_dir, name)
        hhfilter_general(closest_cov_file, qsc_output, "qsc", qsc_value)
        seq_count = count_sequences_in_msa_file(qsc_output)
        print("Diff: ", qsc_value, " / ",seq_count)
        qsc_files.append(qsc_output)

    return qsc_files

def filter_msa_B(closest_cov_file, output_dir, name):
    print("HHFilter B")
    diff_value = 400
    diff_output = generate_output_file_name(closest_cov_file, "filterB", diff_value,  output_dir, name)
    hhfilter_general(closest_cov_file , diff_output, "diff", diff_value)

    seq_count = count_sequences_in_msa_file(diff_output)
    print("Diff: ", diff_value, " / ",seq_count)

    return diff_output

def filter_msa_C(closest_cov_file,  output_dir, name ):
    print("HHFilter C")
    qsc_values = [0.25, 0.5, 0.7]
    diff_value = 400
    qsc_diff_files = []

    for qsc_value in qsc_values:
        qsc_output = generate_output_file_name(closest_cov_file, "filterC", qsc_value,  output_dir, name)
        hhfilter_general(closest_cov_file, qsc_output, "qsc", qsc_value)

        qsc_diff_output = generate_output_file_name(qsc_output, "filterC", diff_value,  output_dir, name)
        hhfilter_general(qsc_output, qsc_diff_output, "diff", diff_value)

        seq_count = count_sequences_in_msa_file(qsc_diff_output)
        print("QSC Diff: ", qsc_value, " / ",seq_count)
        qsc_diff_files.append(qsc_diff_output)

    return qsc_diff_files

def filter_msa_D(closest_cov_file,  output_dir, name):
    print("HHFilter D")
    qsc_values = [0.25, 0.5, 0.7]
    diff_value = 400
    diff_output = generate_output_file_name(closest_cov_file, "filterD", diff_value,  output_dir, name)
    hhfilter_general(closest_cov_file, diff_output, "diff", diff_value)
    diff_qsc_files = []

    for qsc_value in qsc_values:
        diff_qsc_output = generate_output_file_name(diff_output, "filterD", qsc_value,  output_dir, name)
        hhfilter_general(diff_output, diff_qsc_output, "qsc", qsc_value)

        seq_count = count_sequences_in_msa_file(diff_qsc_output)
        diff = abs(seq_count/2 - 2000)
        print("Diff QSC: ", qsc_value, " / ",seq_count)
        diff_qsc_files.append(diff_qsc_output)

    return diff_qsc_files

def execute_filter( a3m_file,  output_dir, name):
    hhfilter_ofile = find_closest_cov(a3m_file,  output_dir, name)
    filter_a_files = filter_msa_A(hhfilter_ofile,  output_dir, name)
    filter_b_file =  filter_msa_B(hhfilter_ofile,  output_dir, name)
    filter_c_files = filter_msa_C(hhfilter_ofile,  output_dir, name)
    filter_d_files = filter_msa_D(hhfilter_ofile,  output_dir, name)
    return hhfilter_ofile #, filter_a_files, filter_b_file, filter_c_files, filter_d_files

