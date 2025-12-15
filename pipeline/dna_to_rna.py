import sys
# python dna_to_rna.py -i <input.fasta> -o <output.fasta>

def error_check(arg_lst):
    if len(arg_lst) != 5:
        print('inproper arguments\npython dna_to_rna.py -i <input.fasta> -o <output.fasta>')
        sys.exit()
    else:
        return True

def read_in_fasta(fasta_file):
        lst = []
        with open(fasta_file, 'r') as fasta_file:
            lst = fasta_file.readlines()
            fasta_file.close()
        fasta_dict = {}
        get_fasta_dict(lst, fasta_dict)
        return fasta_dict

def get_sequence(lst):
        if len(lst) == 0:
            return ''
        elif lst[0][0] == '>':
            return ''
        else:
            line = lst[0]
            line = line.strip()
            return line + get_sequence(lst[1:])

def get_fasta_dict(lst, fasta_dict):
        counter = 0
        for line in lst:
            line = line.strip()
            counter += 1
            if line[0] == '>':
                fasta_dict[line] = get_sequence(lst[counter:])

def dna_to_rna(dna_dict):
    rna_dict = {}
    for seq_id, dna_seq in dna_dict.items():
        rna_seq = dna_seq.upper().replace("T", "U")
        rna_dict[seq_id] = rna_seq
    return rna_dict

def write_rna_file(rna_dict, output_file):
    with open(output_file, 'w') as outfile:
        for seq_id, seq in rna_dict.items():
            outfile.write(f"{seq_id}\n")
            outfile.write(f"{seq}\n")

if __name__ == "__main__":
    arg_lst = sys.argv[:]
    error_check(arg_lst)
    input_file = sys.argv[2]
    output_file = sys.argv[4]
    dna_dict = read_in_fasta(input_file)
    rna_dict = dna_to_rna(dna_dict)
    print(rna_dict)
    write_rna_file(rna_dict, output_file)
