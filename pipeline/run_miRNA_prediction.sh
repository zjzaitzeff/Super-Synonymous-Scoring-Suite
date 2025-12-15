# usage: bash run_miRNA_prediction.sh <miRNA sequence file> <mRNA sequence file> <output directory>
# examples: 
# bash run_miRNA_prediction.sh inputs/miRNATests/BCL2L12/miR-671-5p.fasta inputs/miRNATests/BCL2L12/BCL2L12-seq-isoform6.fasta outputs/miRNATests/BCL2L12/
# bash run_miRNA_prediction.sh inputs/miRNATests/IRGM/miR-196.txt inputs/miRNATests/IRGM/IRGM.fa outputs/miRNATests/IRGM/

# get the inputs from the command line
miRNA=$1
mRNA=$2
output=$3

echo "Running miRNA prediction with inputs:"
echo -e "\tmiRNA file: $miRNA"
echo -e "\tmRNA file: $mRNA"
echo -e "\toutput dir: $output"

# check that the correct number of inputs were given
if [ $# -ne 3 ]; then
    echo "Error: 3 inputs are required"
    exit 1
fi

# check that the inputs are valid
if [ ! -f $miRNA ]; then
    echo "Error: miRNA file not found"
    exit 1
fi
if [ ! -f $mRNA ]; then
    echo "Error: mRNA file not found"
    exit 1
fi

# add a / to the end of the output directory if it doesn't already have one
if [ "${output: -1}" != "/" ]; then
    output="${output}/"
fi
# make the output directory
mkdir -p $output

############################################
# make sure miRNA and mRNA files exist in both TargetScan and FASTA formats
# if the miRNA sequence file is in the TargetScan format, convert it to FASTA format
miRNAExtension="${miRNA##*.}"
inputFolder=$(dirname $miRNA)
miRNAName=$(basename $miRNA ".$miRNAExtension")
if [ $miRNAExtension == "txt" ]; then
    miRNATargetScan=$miRNA
    miRNAFasta=$inputFolder/${miRNAName}.fa
    if [ ! -f $miRNAFasta ]; then
        echo "Creating FASTA miRNA file"
        python3 convertMiRNAFASTATargetScan.py $miRNATargetScan $miRNAFasta
        # exit if the conversion fails
        if [ $? -ne 0 ]; then
            echo "Error: failed to convert miRNA file to FASTA format"
            exit 1
        fi
    fi
# if the miRNA sequence file is in FASTA format (fa, fasta, FASTA), convert it to TargetScan format
elif [ $miRNAExtension == "fa" ] || [ $miRNAExtension == "fasta" ] || [ $miRNAExtension == "FASTA" ] || [ $miRNAExtension == "fna" ]; then
    miRNAFasta=$miRNA
    miRNATargetScan=$inputFolder/${miRNAName}.txt
    if [ ! -f $miRNATargetScan ]; then
        echo "Creating TargetScan miRNA file"
        python3 convertMiRNAFASTATargetScan.py $miRNAFasta $miRNATargetScan
        # exit if the conversion fails
        if [ $? -ne 0 ]; then
            echo "Error: failed to convert miRNA file to TargetScan format"
            exit 1
        fi
    fi
else
    echo "Error: miRNA file must be in FASTA or TargetScan format"
    exit 1
fi

# if the mRNA sequence file is in the TargetScan format, convert it to FASTA format
mRNAExtension="${mRNA##*.}"
inputFolder2=$(dirname $mRNA)
mRNAName=$(basename $mRNA ".$mRNAExtension")
if [ $mRNAExtension == "txt" ]; then
    mRNATargetScan=$mRNA
    mRNAFasta=$inputFolder2/${mRNAName}.fa
    if [ ! -f $mRNAFasta ]; then
        echo "Creating FASTA mRNA file"
        python3 convertMRNAFASTATargetScan.py $mRNATargetScan $mRNAFasta
        # exit if the conversion fails
        if [ $? -ne 0 ]; then
            echo "Error: failed to convert mRNA file to FASTA format"
            exit 1
        fi
    fi
# if the mRNA sequence file is in FASTA format (fa, fasta, FASTA, fna), convert it to TargetScan format
elif [ $mRNAExtension == "fa" ] || [ $mRNAExtension == "fasta" ] || [ $mRNAExtension == "FASTA" ] || [ $mRNAExtension == "fna" ]; then
    mRNAFasta=$mRNA
    mRNATargetScan=$inputFolder2/${mRNAName}.txt
    if [ ! -f $mRNATargetScan ]; then
        echo "Creating TargetScan mRNA file"
        python3 convertMRNAFASTATargetScan.py $mRNAFasta $mRNATargetScan
        # exit if the conversion fails
        if [ $? -ne 0 ]; then
            echo "Error: failed to convert mRNA file to TargetScan format"
            exit 1
        fi
    fi
else
    echo "Error: mRNA file must be in FASTA or TargetScan format"
    exit 1
fi
############################################

echo ""
echo "miRNA FASTA file: $miRNAFasta"
echo "miRNA TargetScan file: $miRNATargetScan"
echo "mRNA FASTA file: $mRNAFasta"
echo "mRNA TargetScan file: $mRNATargetScan"
echo "Output directory: $output"

# run predictions
echo -e "\nStarting miRanda prediction"
(
    start_miranda=$(date +%s)
    miranda $miRNAFasta $mRNAFasta -out "${output}miRanda_output.txt" -sc 50 > "${output}miRanda.log" 2>&1

    if [ $? -ne 0 ]; then
        echo "Error: miRanda prediction failed. See the log file: ${output}miRanda.log"
        exit 1
    else
        end_miranda=$(date +%s)
        echo "miRanda finished in $((end_miranda - start_miranda)) seconds"
    fi
) &

echo -e "\nStarting PITA prediction"
(
    start_pita=$(date +%s)
    # cd to output directory to run PITA
    currentDir=$(pwd)
    cd $output
    "/pipeline/pita/pita_prediction.pl" -utr "/pipeline/$mRNAFasta" -mir "/pipeline/$miRNAFasta" -prefix "pipeline/$output" > "${output}/pita.log" 2>&1
    if [ $? -ne 0 ]; then
        echo "Error: PITA prediction failed. See the log file: ${output}pita.log"
        exit 1
    else
        cd $currentDir
        end_pita=$(date +%s)
        echo "PITA finished in $((end_pita - start_pita)) seconds"
    fi
) &

echo -e "\nStarting TargetScan prediction"
(
    start_targetscan=$(date +%s)
    # if the output file already exists, delete it
    if [ -f "${output}targetscan_output.tsv" ]; then
        rm "${output}targetscan_output.tsv"
    fi
    /pipeline/TargetScan/targetscan_70.pl $miRNATargetScan $mRNATargetScan "${output}targetscan_output.tsv" > "${output}targetscan.log" 2>&1
    if [ $? -ne 0 ]; then
        echo "Error: TargetScan prediction failed. See the log file: ${output}targetscan.log"
        exit 1
    else
        end_targetscan=$(date +%s)
        echo "TargetScan finished in $((end_targetscan - start_targetscan)) seconds"
    fi
) &

# TODO need to auto-interpret the output file
# TODO paccmit-cds expects CDS fasta, not mRNA fasta
# echo -e "\nRunning PACCMIT-CDS prediction"
# (
#     start_paccmit=$(date +%s)
#     currentDir=$(pwd)
#     relativePathToCurrentDir=$(realpath --relative-to="$output" "$currentDir")
#     cd $output
#     "${relativePathToCurrentDir}/paccmit-cds/paccmit-cds" -g "${relativePathToCurrentDir}/$mRNAFasta" -m "${relativePathToCurrentDir}/$miRNAFasta" -i 8 > "${relativePathToCurrentDir}/${output}paccmit.log" 2>&1
#     cd $currentDir
#     end_paccmit=$(date +%s)
#     echo "PACCMIT-CDS finished in $((end_paccmit - start_paccmit)) seconds"
# ) &

# wait for all background processes to finish
wait

# run output interpreter
echo -e "\nRunning output interpreter:"
python interpretOutputs.py $mRNAFasta $output