import sys
import os
import json

def read_in_fasta(fasta_file):
    fasta_dict = {}
    with open(fasta_file, 'r') as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    fasta_dict[header] = "".join(seq_lines)
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if header:
            fasta_dict[header] = "".join(seq_lines)
    return fasta_dict

def get_codons(fasta_dict):
    codon_dict = {}
    for seq_id, seq in fasta_dict.items():
        seq = seq.upper().replace("\n", "").replace(" ", "")
        codons = []
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if len(codon) == 3:
                codons.append(codon)
        codon_dict[seq_id] = codons
    return codon_dict
    
def find_difference(codon_dict, seq1_id, seq2_id):
    seq1 = codon_dict[seq1_id]
    seq2 = codon_dict[seq2_id]
    
    diff1 = []
    diff2 = []
    
    max_len = max(len(seq1), len(seq2))
    
    for i in range(max_len):
        codon1 = seq1[i] if i < len(seq1) else None
        codon2 = seq2[i] if i < len(seq2) else None
        
        if codon1 != codon2:
            diff1.append(codon1)
            diff2.append(codon2)
    return {seq1_id: diff1, seq2_id: diff2}

def get_rel_val(codon_efficiency_file):
    rel_dict = {}
    lines = ''
    with open(codon_efficiency_file, 'r') as csvfile:
        lines = csvfile.readlines()
    for line in lines:
        line = line.split(',')
        rel_dict[line[0]] = line[1]
    return rel_dict

def write_difference_file(output, diff_dict):
    print(f'writing differences to {output_file}')
    with open(output, "w") as f:
        for seq_id, diff_codons in diff_dict.items():
            f.write(f"{seq_id}\n")
            for diff_codon in diff_codons:
                f.write(f'{diff_codon}\n')

def get_rna_dict(dna_dict):
    rna_dict = {}
    for seq_id, diff_codons in dna_dict.items():
        rna_codons = []
        for diff_codon in diff_codons:
            new_codon = ''
            for nucleotide in diff_codon:
                # if nucleotide.upper() == "T":
                #     nucleotide = "U"
                #     new_codon += nucleotide
                # else:
                new_codon += nucleotide.upper()
            rna_codons.append(new_codon)
        rna_dict[seq_id] = rna_codons
    return rna_dict

def write_tsv(diff_dict, rvd, output_tsv, seq1_id, seq2_id):
    print(f'writing tsv to {output_tsv}')
    with open(output_tsv, 'w') as file:
        file.write("codon_wild\tcodon_mutant\tefficency_wild\tefficency_mutant\tdifference_between_effmutant_effwild\n")
        for wild_codon, mutant_codon in zip(diff_dict[seq1_id], diff_dict[seq2_id]):
            eff_wild = float(rvd[wild_codon])
            eff_mut = float(rvd[mutant_codon])
            diff = eff_mut - eff_wild
            file.write(f"{wild_codon}\t{mutant_codon}\t{eff_wild}\t{eff_mut}\t{diff}\n")

if __name__ == "__main__":
    fasta_file = sys.argv[1]
    output_dir = sys.argv[2]
    codon_efficiency_file = sys.argv[3]
    fasta_dict = read_in_fasta(fasta_file)
    seq_ids = list(fasta_dict.keys())
    seq1_id, seq2_id = seq_ids[0], seq_ids[1]
    if len(fasta_dict[seq1_id]) != len(fasta_dict[seq2_id]):
        print('You likely have an insertion or deletion mutation')
        sys.exit(1)
    codon_dict = get_codons(fasta_dict)
    if len(seq_ids) < 2:
        sys.exit(1)
    differences = find_difference(codon_dict, seq1_id, seq2_id)
    if len(differences[seq1_id][0]) != len(differences[seq2_id][0]):
        print(differences[seq1_id][0])
        print(differences[seq2_id][0])
        print(f'there is likely an insertion or deletion mutation')
        sys.exit(1)
    variant_name = os.path.basename(os.path.abspath(output_dir))
    output_file = f"{output_dir}/{variant_name}_codon_diff.fasta"
    output_tsv = f"{output_dir}/PLX_{variant_name}_codon_eff_diff.tsv"
    write_difference_file(output_file, differences)
    rel_val_dict = get_rel_val(codon_efficiency_file)
    rna_dict = get_rna_dict(differences)
    write_tsv(rna_dict, rel_val_dict, output_tsv, seq1_id, seq2_id)