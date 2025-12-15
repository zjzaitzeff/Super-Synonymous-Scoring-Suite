# This script finds all the output directories and runs the summary script on each one

OUTPUT_DIR=$1
if [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: $0 <output_directory>"
    exit 1
fi

# remove the trailing slash if it exists
OUTPUT_DIR="${OUTPUT_DIR%/}"
CLINVAR_VCF=${2:-""}

num_dirs=0
for DIR in $OUTPUT_DIR/*; do
    if [ -d "$DIR" ]; then
        # make sure the dir has at least two dashes in it
        if [[ $(basename "$DIR") == *-*-* ]]; then
            python write_summary.py "$DIR" "$CLINVAR_VCF" > /dev/null
            if [ $? -ne 0 ]; then
                echo "Error: Writing summary failed for $DIR"
                exit 1
            fi
        fi
        num_dirs=$((num_dirs + 1))
        if (( num_dirs % 10 == 0 )); then
            echo "$num_dirs directories processed..."
        fi
    fi
done
echo "All done! Processed $num_dirs directories."