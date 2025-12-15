# This script takes in a VCF file and an output directory and writes a single variant
# VCF file to each subdirectory of the output directory containing the variant chromosome
# and position (e.g. 1-12345)
# This script is meant to be run directly after vcf2fasta.py

import os
import sys

import helper

def read_vcf(vcf_file):
    headers = []
    variants = {}
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                if "contig" in line:
                    # convert contig to chr format if needed
                    chrom_start_index = line.find("ID=") + 3
                    chrom_end_index = line.find(",", chrom_start_index)
                    chrom = line[chrom_start_index:chrom_end_index]
                    new_chrom = helper.convert_chrom_format(chrom, helper.get_chrom_format(chrom), "chr")
                    line = line.replace(f"ID={chrom}", f"ID={new_chrom}")
                headers.append(line)
                continue
            items = line.strip().split('\t')
            chrom = items[0]
            # convert chrom to chr format if needed
            items[0] = helper.convert_chrom_format(chrom, helper.get_chrom_format(chrom), "chr")
            pos = items[1]
            variants[(chrom, pos)] = '\t'.join(items) + '\n'
    return headers, variants

if __name__ == "__main__":
    input_vcf = sys.argv[1]
    output_dir = sys.argv[2]
    # input_vcf = "./inputs/syn_variants.vcf"
    # output_dir = "./outputs/inputs"

    inputs_valid = True
    if not os.path.isfile(input_vcf):
        print(f"Input VCF file {input_vcf} does not exist.")
        inputs_valid = False
    if not os.path.exists(output_dir): # TODO should we just create this directory?
        print(f"Output directory {output_dir} does not exist.")
        inputs_valid = False
    if not inputs_valid:
        print("Exiting.")
        sys.exit(1)
    
    # read in the VCF file
    headers, variants = read_vcf(input_vcf)

    num_vcfs_written = 0
    sub_dirs = sorted([d for d in os.listdir(output_dir) if os.path.isdir(os.path.join(output_dir, d))])
    for sub_dir in sub_dirs:
        sub_dir_path = os.path.join(output_dir, sub_dir)
        if os.path.isdir(sub_dir_path):
            # get the chrom and pos from the subdirectory name
            sub_dir_parts = sub_dir.split('-')
            if len(sub_dir_parts) < 3:
                print(f"Skipping invalid subdirectory name: {sub_dir}")
                continue

            chrom = sub_dir_parts[-2]
            pos = sub_dir_parts[-1]

            # check that chrom and pos are valid
            if not helper.is_chrom_valid(chrom) or not pos.isdigit():
                continue
            
            # check that the variant exists in the VCF
            if not (chrom, pos) in variants:
                print(f"Skipping missing variant: {chrom}-{pos}")
                continue

            # if the headers do not contain contig lines, add them here
            if not any("##contig" in h for h in headers):
                headers.insert(-1, f"##contig=<ID={chrom}>\n")

            # write out the single variant VCF file
            variant_vcf_path = os.path.join(sub_dir_path, f"{sub_dir}.vcf")
            with open(variant_vcf_path, 'w') as outf:
                outf.writelines(headers)
                outf.write(variants[(chrom, pos)])
            # print(f"Wrote single variant VCF: {variant_vcf_path}")
            num_vcfs_written += 1
    
    print(f"Wrote {num_vcfs_written} single variant VCF files.")