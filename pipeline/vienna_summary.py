import os
import sys
import re

def read_lines(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file]

def get_mutant_vs_wild(lines):
    wild_lines = []
    mutant_lines = []

    seq_start_indices = [i for i, line in enumerate(lines) if line.startswith('>')]

    if len(seq_start_indices) < 2:
        raise ValueError("FASTA must contain at least two sequences and the wild sequence should be first")
    wild_lines = lines[seq_start_indices[0]:seq_start_indices[1]]
    mutant_lines = lines[seq_start_indices[1]:]
    return wild_lines, mutant_lines

def get_fold_values(lines):
    MFE = float(lines[2][lines[2].index(" "):].strip()[1:-1].strip()) # value at the end of line 3 in parentheses
    EFE = float(lines[3][lines[3].index(" "):].strip()[1:-1].strip()) # value at the end of line 4 in square brackets
    line5_items = lines[4][lines[4].index(" "):].strip()[1:-1].strip().split()
    CFE = float(line5_items[0]) # value at the end of line 5 in braces before "d="
    CD = float(line5_items[1][2:]) # value at the end of line 5 in braces after "d="
    MEAFE = float(lines[5][lines[5].index(" "):].strip()[1:-1].strip().split()[1][4:]) # value at the end of line 6 in braces after "MEA="
    END = float(lines[6].split()[-1]) # value at the end of line 7
    return MFE, EFE, CFE, CD, MEAFE, END

def make_dict(mfe, efe, cfe, cd, meafe, end):
    d = {}
    d['MFE'] = mfe
    d['EFE'] = efe
    d['CFE'] = cfe
    d['CD'] = cd
    d['MEAFE'] = meafe
    d['END'] = end
    return d

def read_fold_file(vienna_fold_path):
    if not os.path.isfile(vienna_fold_path):
        print(f"WARNING: '{vienna_fold_path}' does not exist")
        return None, None, None
    lines = read_lines(vienna_fold_path)
    wild_lines, mutant_lines = get_mutant_vs_wild(lines)
    wild_MFE, wild_EFE, wild_CFE, wild_CD, wild_MEAFE, wild_END = get_fold_values(wild_lines)
    mutant_MFE, mutant_EFE, mutant_CFE, mutant_CD, mutant_MEAFE, mutant_END = get_fold_values(mutant_lines)
    deltaMFE = round(wild_MFE - mutant_MFE, 2)
    deltaEFE = round(wild_EFE - mutant_EFE, 2)
    deltaCFE = round(wild_CFE - mutant_CFE, 2)
    deltaCD = round(wild_CD - mutant_CD, 2)
    deltaMEAFE = round(wild_MEAFE - mutant_MEAFE, 2)
    deltaEND = round(wild_END - mutant_END, 2)
    wild_dict = make_dict(wild_MFE, wild_EFE, wild_CFE, wild_CD, wild_MEAFE, wild_END)
    mutant_dict = make_dict(mutant_MFE, mutant_EFE, mutant_CFE, mutant_CD, mutant_MEAFE, mutant_END)
    delta_dict = make_dict(deltaMFE, deltaEFE, deltaCFE, deltaCD, deltaMEAFE, deltaEND)
    return wild_dict, mutant_dict, delta_dict

def read_pdist_file(vienna_pdist_path):
    if not os.path.isfile(vienna_pdist_path):
        print(f"WARNING: '{vienna_pdist_path}' does not exist")
        return None
    lines = read_lines(vienna_pdist_path)
    return float(lines[2])

def read_distance_file(vienna_distance_path):
    if not os.path.isfile(vienna_distance_path):
        print(f"WARNING: '{vienna_distance_path}' does not exist")
        return None, None, None
    lines = read_lines(vienna_distance_path)
    MFEED = float(lines[1].replace("f:", ""))
    CFEED = float(lines[3].replace("f:", ""))
    MEAED = float(lines[4].replace("f:", ""))
    return MFEED, CFEED, MEAED

def write_summary(MFEED, CFEED, MEAED, EFEED, wild_dict, mutant_dict, delta_dict, input_folder_path, out_path):
    # MFE EFE CFE CD  MEAFE   END MFEED   CFEED   MEAED EFEED
    with open(out_path, 'w') as file:
        file.write('type\tMFE\tEFE\tCFE\tCD\tMEAFE\tEND\tMFEED\tCFEED\tMEAED\tEFEED\n')
        file.write(f"REF\t{wild_dict['MFE']}\t{wild_dict['EFE']}\t{wild_dict['CFE']}\t{wild_dict['CD']}\t{wild_dict['MEAFE']}\t{wild_dict['END']}\tsee mutant\tsee mutant\tsee mutant\tsee mutant\n")
        file.write(f"ALT\t{mutant_dict['MFE']}\t{mutant_dict['EFE']}\t{mutant_dict['CFE']}\t{mutant_dict['CD']}\t{mutant_dict['MEAFE']}\t{mutant_dict['END']}\t{MFEED}\t{CFEED}\t{MEAED}\t{EFEED}\n")
        file.write(f"wild-mutant\t{delta_dict['MFE']}\t{delta_dict['EFE']}\t{delta_dict['CFE']}\t{delta_dict['CD']}\t{delta_dict['MEAFE']}\t{delta_dict['END']}\n")

if __name__ == "__main__":
    input_folder_path = sys.argv[1]
    out_path = os.path.join(input_folder_path, 'PLX_vienna_summary.tsv')

    vienna_fold_path = os.path.join(input_folder_path, 'vienna_fold.fold')
    vienna_pdist_path = os.path.join(input_folder_path, 'vienna_pdist.txt')
    vienna_distance_path = os.path.join(input_folder_path, 'vienna_distance.txt')

    wild_dict, mutant_dict, delta_dict = read_fold_file(vienna_fold_path)
    EFEED = read_pdist_file(vienna_pdist_path)
    MFEED, CFEED, MEAED = read_distance_file(vienna_distance_path)

    write_summary(MFEED, CFEED, MEAED, EFEED, wild_dict, mutant_dict, delta_dict, input_folder_path, out_path)
    print(f"Wrote summary to '{out_path}'")
