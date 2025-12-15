import sys
from pathlib import Path

if __name__ == "__main__":
    # TODO make clinvar optional, especially when the species is not human
    # input_folder = "./outputs"
    # output = "./outputs/master_summary_v3.tsv"

    # input_folder = sys.argv[1].rstrip("/")
    # output = sys.argv[2]

    # TODO remove
    input_folder = "./outputs_jag1"
    output = "./outputs_jag1/master_summary.tsv"

    summary_files = list(Path(input_folder).glob("*/summary.tsv"))
    summary_paths = [str(f) for f in summary_files]
    num_files_processed = 0

    header_written = False
    with open(output, "w") as outf:
        for file in summary_paths:
            with open(file, "r") as f:
                header = f.readline()
                if not header_written:
                    outf.write(header)
                    header_written = True
                outf.write(f.readline().strip() + "\n")
            num_files_processed += 1
            if num_files_processed % 100 == 0:
                print(f"{num_files_processed} files processed...")
        print(f"{num_files_processed} files processed.")
    print(f"Wrote master summary to '{output}'")