import os
import sys
import shutil
import re
results_dir = sys.argv[1]

def main(results_dir):
    pattern = r"^sorted\.(.+)\.sequence_class_scores\.tsv$"

    found = False
    for filename in os.listdir(results_dir):
        match = re.match(pattern, filename)
        if match:
            program_name = match.group(1)
            src = os.path.join(results_dir, filename)
            dst = os.path.join(results_dir, f"PLX_sorted.{program_name}.sequence_class_scores.tsv")
            shutil.copy2(src, dst)
            found = True
            break

    if not found:
        print("No file matching pattern 'sorted.<program>.sequence_class_scores.tsv' found.")

if __name__ == "__main__":
    main(results_dir)