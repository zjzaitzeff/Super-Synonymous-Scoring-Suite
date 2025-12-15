# This script takes in the output directory of the pipeline and checks if each subdirectory
# has a summary.tsv file.

import os
import sys


if __name__ == "__main__":
    results_dir = sys.argv[1]

    missing_summary = 0
    has_summary = 0
    for folder in os.listdir(results_dir):
        folder_path = os.path.join(results_dir, folder)
        if os.path.isdir(folder_path):
            # only check the folder if the nanem contains at least 2 hyphens
            if folder.count('-') >= 2:
                # check if the directory contains a summary.tsv file
                summary_path = os.path.join(folder_path, "summary.tsv")
                if not os.path.isfile(summary_path): 
                    missing_summary += 1
                else:
                    has_summary += 1

    print(f"Missing summary files: {missing_summary}")
    print(f"Found summary files: {has_summary}")
