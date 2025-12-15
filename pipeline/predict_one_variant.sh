#!/bin/bash
#bash predict_one_variant.sh /Users/zkjz/Desktop/Coding/Bioinformatics/pipeline_docker_prac/inputs_test /Users/zkjz/Desktop/Coding/Bioinformatics/pipeline_docker_prac/outputs_test ref/miR_Family_Info_Human.fa 500
#dir of variant
INPUT_DIR=$1
#dir where you want it go
OUTPUT_DIR=$2
#ref/miR_Family_Info_Human.fa
MOTIF_FILE=$3
#500
truncate_mRNA_num=$4
#can leave blank
CLINVAR_VCF=$5

# remove trailing slashes if they exist
INPUT_DIR="${INPUT_DIR%/}"
OUTPUT_DIR="${OUTPUT_DIR%/}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p $OUTPUT_DIR

#inputfile creates the three kinds of files we want
echo "Running pipeline for: $INPUT_DIR"
dirname=${INPUT_DIR##*/}

single_var_vcf_file="${INPUT_DIR}/${dirname}.vcf"
mRNA_seq_file=$(find "$INPUT_DIR" -maxdepth 1 -type f -name "mRNA-*.fa" ! -name "*trunc*")
CDS_seq_file=$(find "$INPUT_DIR" -maxdepth 1 -type f -name "CDS-*.fa")

# truncate mRNA to the 1000 bases surrounding the variant
# TODO make this part optional
echo ""
echo "Truncating mRNA to $truncate_mRNA_num bases surrounding the variant"
python truncate_mRNA.py "$mRNA_seq_file" "$truncate_mRNA_num" "$OUTPUT_DIR"
# get the name of the truncated mRNA file
mRNA_seq_file_name=$(basename "$mRNA_seq_file" .fa)
mRNA_seq_file="$OUTPUT_DIR/${mRNA_seq_file_name}_trunc${truncate_mRNA_num}.fa"
# check if the file exists
if [ ! -f "$mRNA_seq_file" ]; then
    echo "Error: truncated mRNA file not found: $mRNA_seq_file"
    exit 1
fi

echo -e "\nsingle variant vcf file path: '$single_var_vcf_file'"
echo "mRNA file path: '$mRNA_seq_file'"
echo "CDS file path: '$CDS_seq_file'"
echo "Script dir: '$SCRIPT_DIR'"

# single_var_vcf_file scripts
# TODO test
echo -e "\nStarting splice site prediction with SpliceAI"
(
    start_spliceai=$(date +%s)
    spliceai_log="${OUTPUT_DIR}/log_SpliceAI.log"
    cd $SCRIPT_DIR/spliceai
    # made all paths absolute inside the container
    bash spliceai_snp.sh "$single_var_vcf_file" "$OUTPUT_DIR/PLX_spliceai_output.vcf" > "$spliceai_log" 2>&1
    if [ $? -ne 0 ]; then
        echo "Error: SpliceAI prediction failed. See the log file: $spliceai_log" | tee -a "$spliceai_log"
        exit 1
    else
        end_spliceai=$(date +%s)
        echo "SpliceAI prediction finished in $((end_spliceai - start_spliceai)) seconds" | tee -a "$spliceai_log"
        cd "$SCRIPT_DIR"
    fi
) &

# sei (TFBS)
# TODO test
echo -e "\nStarting TFBS prediction with Sei"
(
    start_tfbs=$(date +%s)
    tfbs_log="${OUTPUT_DIR}/log_TFBS.log"
    cd $SCRIPT_DIR/sei/sei-framework-main
    # I switched everything to an absolute path inside the container
    echo "Moved to sei dir: $(pwd)" > "$tfbs_log" 2>&1
    bash run-variantEffectScores.sh "$single_var_vcf_file" "hg38" "$OUTPUT_DIR" >> "$tfbs_log" 2>&1
    if [ $? -ne 0 ]; then
        echo "Error: TFBS prediction failed. See the log file: $tfbs_log" | tee -a "$tfbs_log"
        exit 1
    else
        cd "$SCRIPT_DIR"
        echo "Moved back to main dir: $(pwd)" >> "$tfbs_log" 2>&1
        python rename_plx.py "$OUTPUT_DIR"
        end_tfbs=$(date +%s)
        echo "TFBS prediction finished in $((end_tfbs - start_tfbs)) seconds" | tee -a "$tfbs_log"
    fi
) &

# mRNA_seq_file
# miRNA
# TODO test
echo -e "\nStarting miRNA prediction with miRanda, PITA, and TargetScan"
(
    start_mirna=$(date +%s)
    # TODO make the convertMRNAFASTATargetScan.py and convertMiRNAFASTATargetScan.py scripts the same file
    mirna_log="${OUTPUT_DIR}/log_miRNA.log"
    bash run_miRNA_prediction.sh "$MOTIF_FILE" "$mRNA_seq_file" "$OUTPUT_DIR" > "$mirna_log" 2>&1
    if [ $? -ne 0 ]; then
        echo "Error: miRNA prediction failed. See the log file: $mirna_log" | tee -a "$mirna_log"
        exit 1
    else
        end_mirna=$(date +%s)
        echo "miRNA prediction finished in $((end_mirna - start_mirna)) seconds" | tee -a "$mirna_log"
    fi
) &

# mRNA secondary structure with vienna
echo -e "\nStarting mRNA secondary structure prediction with ViennaRNA"
(
    start_mrna_structure=$(date +%s)
    mrna_structure_log="${OUTPUT_DIR}/log_mrna_structure.log"
    bash run_vienna.sh "$mRNA_seq_file" "$OUTPUT_DIR" >> "$mrna_structure_log" 2>&1
    if [ $? -ne 0 ]; then
        echo "Error: mRNA secondary structure prediction failed. See the log file: $mrna_structure_log"  | tee -a "$mrna_structure_log"
        exit 1
    else
        end_mrna_structure=$(date +%s)
        echo "mRNA secondary structure prediction finished in $((end_mrna_structure - start_mrna_structure)) seconds"  | tee -a "$mrna_structure_log"
    fi
) &

# CDS_seq_file
# codon_efficiency
# TODO test
echo -e "\nStarting codon efficiency predictions"
(
    codon_efficiency_file="ref/GRCh38_longest_isoforms-wij.csv"
    start_codon_efficiency=$(date +%s)
    codon_efficiency_log="${OUTPUT_DIR}/log_codon_efficiency.log"
    echo "Predicting codon efficiency" > "$codon_efficiency_log" 2>&1
    # check if the CDS file exists
    if [ ! -f "$CDS_seq_file" ]; then
        echo "No CDS file found. Skipping." | tee -a "$codon_efficiency_log"
        exit 1
    fi

    python find_codon_efficiency.py "$CDS_seq_file" "$OUTPUT_DIR" "$codon_efficiency_file" >> "$codon_efficiency_log" 2>&1
    if [ $? -ne 0 ]; then
        echo "Error: Codon efficiency prediction failed. See the log file: $codon_efficiency_log" | tee -a "$codon_efficiency_log"
    fi

    echo "" >> "$codon_efficiency_log" 2>&1
    python ExtRamp.py -v -i "$CDS_seq_file" -o "$OUTPUT_DIR/PLX_extramp.fasta" -a "$codon_efficiency_file" >> "$codon_efficiency_log" 2>&1
    if [ $? -ne 0 ]; then
        echo "Error: ExtRamp prediction failed. See the log file: $codon_efficiency_log" | tee -a "$codon_efficiency_log"
        exit 1
    else
        end_codon_efficiency=$(date +%s)
        echo "Codon efficiency prediction finished in $((end_codon_efficiency - start_codon_efficiency)) seconds" | tee -a "$codon_efficiency_log"
    fi
) &

# wait for all background processes to finish
wait

# summary
python write_summary.py "$OUTPUT_DIR" "$CLINVAR_VCF"

exit 0