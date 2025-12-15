# This script takes a mRNA FASTA file with a REF and ALT sequence and a number and 
# truncates the mRNA sequences to the specified length surounding the single mutation in the ALT.
# If the start or end indices are out of bounds, it adjusts them to fit within the sequence,
# which may result in truncated sequences shorter than the specified length.

import sys
import os

def read_fasta(in_path):
    ref_header = ""
    alt_header = ""
    ref_sequence = ""
    alt_sequence = ""
    with open(in_path, "r") as in_file:
        lines = in_file.readlines()
        # remove any empty lines
        lines = [line for line in lines if line.strip()]
        if len(lines) != 4:
            print(f"Error: FASTA file '{in_path}' does not contain two sequences.")
            sys.exit(1)
        ref_header = lines[0].strip()
        ref_sequence = lines[1].strip()
        alt_header = lines[2].strip()
        alt_sequence = lines[3].strip()
    return ref_header, ref_sequence, alt_header, alt_sequence

def get_diff_index(ref_sequence, alt_sequence):
    for i in range(min(len(ref_sequence), len(alt_sequence))):
        if ref_sequence[i] != alt_sequence[i]:
            return i
    return -1

if __name__ == "__main__":
    in_path = sys.argv[1]
    trunc_length = int(sys.argv[2])
    # optional output directory, defaults to input file directory
    out_dir = sys.argv[3] if len(sys.argv) > 3 else os.path.dirname(in_path)
    out_path = os.path.join(out_dir, os.path.basename(in_path).replace(".fa", f"_trunc{trunc_length}.fa"))

    if not os.path.isfile(in_path):
        print(f"Error: input file '{in_path}' not found")
        sys.exit(1)

    # read in the two sequences
    ref_header, ref_sequence, alt_header, alt_sequence = read_fasta(in_path)
    print(f"Sequences are {len(ref_sequence)} (REF) and {len(alt_sequence)} (ALT) bases long")

    # find the position of the difference between the two sequences
    diff_index = get_diff_index(ref_sequence, alt_sequence)
    if diff_index == -1:
        print("Error: No difference found between REF and ALT sequences")
        sys.exit(1)        
    print(f"Found difference at index {diff_index}")

    # truncate the sequences
    start_index = max(0, diff_index - trunc_length // 2)
    end_index = start_index + trunc_length
    if end_index > len(ref_sequence):
        end_index = len(ref_sequence)
    if start_index < 0:
        start_index = 0
    ref_sequence_trunc = ref_sequence[start_index:end_index]
    alt_sequence_trunc = alt_sequence[start_index:end_index]

    print(f"Truncated sequences to {len(ref_sequence_trunc)} bases")

    # write out the truncated sequences
    with open(out_path, "w") as out_file:
        out_file.write(f"{ref_header}\n{ref_sequence_trunc}\n")
        out_file.write(f"{alt_header}\n{alt_sequence_trunc}\n")
    
    print(f"Wrote truncated sequences to '{out_path}'")