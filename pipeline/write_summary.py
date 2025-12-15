import sys
import re
import os
import pysam

# Get variant information from VCF file
def get_variant_info(vcf_file):
    with open(vcf_file, 'r') as vcf:
        line = vcf.readline()
        while line.startswith('#'):
            line = vcf.readline()
    items = line.strip().split('\t')
    chrom = items[0]
    pos = int(items[1])
    rsID = items[2]
    ref = items[3]
    alt = items[4]
    info = items[7]
    return chrom, pos, rsID, ref, alt, info

# Get ClinVar significance for a given variant. 
# Options are "Not in ClinVar", "No sig in ClinVar" (ie, present in Clinvar, but no significance listed), 
# or comma-separated string of the significance values found
def get_clinvar_significance(clinvar_path, chrom, pos, rsID, ref, alt):
    # TODO alt should only be a single allele here, but may be more until I fix vcf2fasta
    sigs = set()

    if chrom.startswith("chr"):
        chrom = chrom[3:]

    clinvar_vcf = pysam.VariantFile(clinvar_path)
    clinvar_records = clinvar_vcf.fetch(chrom, pos - 1, pos)
    clinvar_records = [record for record in clinvar_records]
    found = False
    for record in clinvar_records:
        clinvar_ref = record.ref
        clinvar_alts = list(record.alts)
        if ref == clinvar_ref and alt in clinvar_alts:
            found = True
            sig = record.info.get("CLNSIG")
            if sig:
                if isinstance(sig, (list, tuple)):
                    sigs.update(sig)
                else:
                    sigs.add(sig)
    if not found:
        return "Not in ClinVar"
    sigs = ",".join(sigs) if sigs else "No sig in ClinVar"
    return sigs

# miRanda and PITA have the same interpretation file format, only the scores are formatted differently
# Return "Different miRNAs" if the second line states there are differing miRNAs
# Return "Same miRNAs, different scores" if the second line states there are no differing miRNAs, but there are
# miRNAs with different scores. Otherwise return "No"
def summarize_miRanda_pita(file_path):
    with open(file_path, 'r') as file:
        header = file.readline()
        data_line = file.readline()
        if not data_line.strip():
            return "May have failed, check manually"
        num_differing_miRNAs = int(data_line.strip().split('\t')[3])
        if num_differing_miRNAs > 0:
            return "Different miRNAs"
        other_lines = file.readlines()
    if len(other_lines) < 5:
        return "May have failed, check manually"
    if other_lines[4].strip() == "All scores are the same.":
        return "No"
    else:
        return "Same miRNAs, different scores"

# TargetScan doesn't have scores, just the second line in the interpretation that says
# if there are differing miRNAs. If the second line is empty, neither the REF nor the ALT
# have any predicted miRNA binding sites
def summarize_targetscan(file_path):
    with open(file_path, 'r') as file:
        header = file.readline()
        data_line = file.readline()
        if not data_line.strip():
            return "No"
    num_differing_miRNAs = int(data_line.strip().split('\t')[3])
    if num_differing_miRNAs > 0:
        return "Different miRNAs"
    else:
        return "No"

def summarize_miRNA_agreement(file_path):
    miRNA_two_tool_agreement = set()
    miRNA_three_tool_agreement = set()
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            num_tools = 0
            shared_miRNAs_str = line.split(":")[1]
            if not shared_miRNAs_str:
                continue
            shared_miRNAs = shared_miRNAs_str.split("|")
            if "TargetScan" in line:
                num_tools += 1
            if "PITA" in line:
                num_tools += 1
            if "miRanda" in line:
                num_tools += 1
            if num_tools == 2:
                miRNA_two_tool_agreement.update(shared_miRNAs)
            elif num_tools == 3:
                miRNA_three_tool_agreement.update(shared_miRNAs)
    # remove miRNAs that are in the three tool agreement from the two tool agreement
    miRNA_two_tool_agreement = miRNA_two_tool_agreement - miRNA_three_tool_agreement
    return len(miRNA_two_tool_agreement), len(miRNA_three_tool_agreement)

def int_sei(file_path):
    #return 'yes' if score > 0.18 return 'no' if score < 0.18
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    lines = lines[1:]
    for line in lines:
        line = line.split('\t')
        if float(line[0]) > 0.18:
            return float(line[0])
    return "No"

def int_spliceai(file_path):
    with open(file_path, 'r') as inf:
        for line in inf:
            if line.startswith('#'):
                continue
            items = line.strip().split('\t')
            info_field = items[7]
            info_items = info_field.split(';')
            for item in info_items:
                if item.startswith("SpliceAI="):
                    spliceai_items = item.split('|')
                    acceptor_gain = float(spliceai_items[2])
                    acceptor_loss = float(spliceai_items[3])
                    donor_gain = float(spliceai_items[4])
                    donor_loss = float(spliceai_items[5])
                    max_change = max(acceptor_gain, acceptor_loss, donor_gain, donor_loss)
                    # get which change is the max
                    if max_change == acceptor_gain:
                        return f"a:{max_change}"
                    elif max_change == acceptor_loss:
                        return f"a:{max_change * -1}"
                    elif max_change == donor_gain:
                        return f"d:{max_change}"
                    elif max_change == donor_loss:
                        return f"d:{max_change * -1}"
    return "NA"

def int_vienna(file_path):
    #currently do not have thresholds for vienna, just return bar separated values
    with open(file_path, 'r') as inf:
        header = inf.readline()
        ref_line = inf.readline()
        alt_line = inf.readline()
        alt_line_items = alt_line.strip().split('\t')
        diff_line = inf.readline()
        diff_line_items = diff_line.strip().split('\t')
        delta_MFE = diff_line_items[1]
        delta_EFE = diff_line_items[2]
        delta_CFE = diff_line_items[3]
        delta_CDE = diff_line_items[4]
        delta_MEAFE = diff_line_items[5]
        delta_END = diff_line_items[6]
        MFEED = alt_line_items[7]
        CFEED = alt_line_items[8]
        MEAD = alt_line_items[9]
        EFEED = alt_line_items[10]
        return f"ΔMFE={delta_MFE}|ΔEFE={delta_EFE}|ΔCFE={delta_CFE}|ΔCDE={delta_CDE}|ΔMEAFE={delta_MEAFE}|ΔEND={delta_END}|MFEED={MFEED}|CFEED={CFEED}|MEAD={MEAD}|EFEED={EFEED}"
        
    return 'Check manually'

def int_codon_diff(file_path):
    # currently do not have thresholds for codon efficiency differences, just return difference
    with open(file_path, 'r') as inf:
        header = inf.readline()
        data_line = inf.readline()
        line_items = data_line.strip().split('\t')
        return line_items[-1]

def read_in_fasta(fasta_file):
    fasta_dict = {}
    with open(fasta_file, 'r') as f:
        header = ""
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

def int_extramp(file_path):
    fasta_dict = read_in_fasta(file_path)
    if not fasta_dict:
        return "No"
    elif len(fasta_dict) == 1:
        header = list(fasta_dict.keys())[0]
        if header.startswith(">REF"):
            return "ramp loss"
        if header.startswith(">ALT"):
            return "ramp gain"
        return "Yes"
    else:
        headers = list(fasta_dict.keys())
        seq1 = fasta_dict[headers[0]]
        seq2 = fasta_dict[headers[1]]
        if len(seq1) != len(seq2):
            return "ramp length change"
        else:
            return "No"

#this code looks into a given folder path and finds all the files with PLX at the beginning and are a .txt
def find_files(folder_path):
    pattern = re.compile(r"^PLX.*\.(txt|tsv|fasta|fa|vcf)$")

    if not os.path.exists(folder_path) or not os.path.isdir(folder_path):
        print(f"Invalid folder path: {folder_path}")
        return None

    lst_of_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f)) and pattern.match(f)]
    return lst_of_files

if __name__ == '__main__':
    # find interpreted files in the given directory
    folder_input = sys.argv[1]
    clinvar_path = sys.argv[2] if len(sys.argv) > 2 else "" # TODO make sure this is added to pipeline!

    lst_of_files = find_files(folder_input)
    print(lst_of_files)
    if not lst_of_files:
        print("No files found, exiting")
        sys.exit(1)
    out_path = os.path.join(folder_input, "summary.tsv")

    variant_id = folder_input.rstrip('/').split('/')[-1]
    # get variant info from single variant VCF
    vcf_path = os.path.join(folder_input, f"{variant_id}.vcf")
    chrom, pos, rsID, ref, alt, info = get_variant_info(vcf_path)

    # get significance from clinvar
    clinvar_significance = get_clinvar_significance(clinvar_path, chrom, pos, rsID, ref, alt) if clinvar_path else "NA"

    targetscan_summary = 'NA'
    pita_summary = 'NA'
    miRanda_summary = 'NA'
    miRNA_two_tool_agreement = 'NA'
    miRNA_three_tool_agreement = 'NA'
    sei_summary = 'NA'
    spliceai_summary = 'NA'
    vienna_summary = 'NA'
    codon_eff_summary = 'NA'
    extramp_summary = 'NA'
    for filename in lst_of_files:
        if re.search(r"targetscan.*\.txt$", filename):
            targetscan_summary = summarize_targetscan(filename)
        elif re.search(r"pita.*\.txt$", filename):
            pita_summary = summarize_miRanda_pita(filename)
        elif re.search(r"miRanda.*\.txt$", filename):
            miRanda_summary = summarize_miRanda_pita(filename)
        elif re.search(r"sorted\.(.+)\.sequence_class_scores\.tsv$", filename):
            sei_summary = int_sei(filename)
        elif re.search(r"spliceai.*\.vcf$", filename):
            spliceai_summary = int_spliceai(filename)
        elif re.search(r"vienna.*\.(txt|tsv)$", filename):
            vienna_summary = int_vienna(filename)
        elif re.search(r"codon_eff_diff.*\.(txt|tsv)$", filename):
            codon_eff_summary = int_codon_diff(filename)
        elif re.search(r"extramp.*\.(fasta|fa)$", filename):
            extramp_summary = int_extramp(filename)
        elif re.search(r"PLX_miRNA_prediction_overlaps\.txt$", filename):
            miRNA_two_tool_agreement, miRNA_three_tool_agreement = summarize_miRNA_agreement(filename)
        else:
            print(f'{filename} could not be interpreted correctly')
    
    if codon_eff_summary == "NA" or extramp_summary == "NA":
        codon_eff_log = os.path.join(folder_input, "log_codon_efficiency.log")
        if os.path.isfile(codon_eff_log):
            with open(codon_eff_log, 'r') as logf:
                log_contents = logf.read()
            # if the last line contains "No CDS file found.", then set codon_eff_summary to "No CDS file"
            if "No CDS file found." in log_contents:
                codon_eff_summary = "No CDS file"
                extramp_summary = "No CDS file"
            elif "No sequences passed the initial filter." in log_contents:
                extramp_summary = "CDS too short or not multiple of 3"

    with open(out_path, 'w') as outf:
        outf.write(f'variant\tchrom\tpos\trsID\tref\talt\tinfo\tTargetScan\tPITA\tmiRanda\tmiRNA two tool agreement\tmiRNA three tool agreement\tsei\tspliceAI\tVienna\tCodon Efficiency\tExtRamp\tClinVar\n')
        outf.write(f'{variant_id}\t{chrom}\t{pos}\t{rsID}\t{ref}\t{alt}\t{info}\t{targetscan_summary}\t{pita_summary}\t{miRanda_summary}\t{miRNA_two_tool_agreement}\t{miRNA_three_tool_agreement}\t{sei_summary}\t{spliceai_summary}\t{vienna_summary}\t{codon_eff_summary}\t{extramp_summary}\t{clinvar_significance}')

    print('summary written')

