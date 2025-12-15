input=$1
outdir=$2
RNAFOLD="$outdir/vienna_fold.fold"
RNAPDIST="$outdir/vienna_pdist.txt"
RNADISTANCE="$outdir/vienna_distance.txt"

# TODO what does this do?
add_suffix_before_ext () {
    local filename="$1"
    local suffix="$2"
    local base="${filename%.*}"
    local ext="${filename##*.}"
    if [[ "$base" == "$filename" ]]; then
        echo "${filename}${suffix}"
    else
        echo "${base}${suffix}.${ext}"
    fi
}

if [ -z "$input" ] || [ -z "$RNAFOLD" ] || [ -z "$RNAPDIST" ] || [ -z "$RNADISTANCE" ]; then
    echo "Error: One or more input files are missing."
    exit 1
fi

python vienna_input_verify.py $RNAFOLD $RNAPDIST $RNADISTANCE
RNAfold -p --MEA -i $input > $RNAFOLD

# if RNAFOLD is empty
if [ ! -s "$RNAFOLD" ]; then
    echo "Error: RNAfold output is empty."
    exit 1
fi

RNApdist < $input > $RNAPDIST
RNAdistance < $RNAFOLD > $RNADISTANCE

echo "Summarizing ViennaRNA results..."
python vienna_summary.py $outdir
