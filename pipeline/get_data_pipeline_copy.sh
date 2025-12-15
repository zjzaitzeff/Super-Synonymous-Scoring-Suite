#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=2G
#SBATCH -J "synonymousEffectsPipeline"
#SBATCH --output=slurm-%x.%j.out
#SBATCH --mail-user=cloward9@byu.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# VCF -> SpliceAI
# gene FASTA -> TFBS
# mRNA FASTA -> miRNA, mRNA ss
# CDS -> codon efficiency, ExtRamp

# INPUT_FILE=inputs/syn_variants.vcf
# INPUTS_FOLDER=outputs/inputs
# MOTIF_FILE=ref/miR_Family_Info_Human.fa
# truncate_mRNA_num=500
# SLURM=True
# TIME_PER_VARIANT="00:30:00"

INPUT_FILE=${1:-"inputs/syn_variants.vcf"}
INPUTS_FOLDER=${2:-"outputs/inputs"}
OUTPUTS_FOLDER=${3:-"outputs"}
MOTIF_FILE=${4:-"ref/miR_Family_Info_Human.fa"}
truncate_mRNA_num=${5:-500}
CLINVAR_VCF=${6:-""}
SLURM=${7:-True}
TIME_PER_VARIANT=${8:-"00:30:00"}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# loop through each directory in inputs
for INPUT_DIR in $INPUTS_FOLDER/*; do
    if [[ -d "$INPUT_DIR" ]]; then
        dirname=${INPUT_DIR##*/}
        OUTPUT_DIR="$OUTPUTS_FOLDER/$dirname"
        mkdir -p $OUTPUT_DIR
        mkdir -p $OUTPUT_DIR/slurm
        
        # if slurm
        if [ "$SLURM" = True ]; then
            slurm_file_path="$OUTPUT_DIR/slurm_predict_one_variant.sh"
            echo "#!/bin/bash" > "$slurm_file_path"
            echo "#SBATCH --time=$TIME_PER_VARIANT" >> "$slurm_file_path"
            echo "#SBATCH --ntasks=4" >> "$slurm_file_path"
            echo "#SBATCH --nodes=1" >> "$slurm_file_path"
            echo "#SBATCH --mem-per-cpu=8G" >> "$slurm_file_path"
            echo "#SBATCH -J \"predict-$dirname\"" >> "$slurm_file_path"
            echo "#SBATCH --output=$OUTPUT_DIR/slurm/slurm-%x.%j.out" >> "$slurm_file_path"
            echo "" >> "$slurm_file_path"
            echo "bash predict_one_variant.sh \"$INPUT_DIR\" \"$OUTPUT_DIR\" \"$MOTIF_FILE\" \"$truncate_mRNA_num\" \"$CLINVAR_VCF\"" >> "$slurm_file_path"
            sbatch "$slurm_file_path"
            echo -e "\t$slurm_file_path"
        else
            bash predict_one_variant.sh "$INPUT_DIR" "$OUTPUT_DIR" "$MOTIF_FILE" "$truncate_mRNA_num" "$CLINVAR_VCF"
        fi
    fi
done

# TODO make a master summary file at the end
echo -e "\nAll done!"