import sys

FASTA_FILE_EXTENSIONS = ["fa", "fasta", "fna", "FASTA"]

# TODO test this function
def convertTargetScanToFASTA(inPath, outPath):
    with open(inPath) as inf, open(outPath, "w") as outf:
        for line in inf:
            lineItems = line.strip().split("\t")
            header = lineItems[0]
            seq = lineItems[2]
            outf.write(f">{header}\n{seq}\n")

# TODO test this function
def convertFASTAToTargetScan(inPath, outPath):
    with open(inPath) as inf, open(outPath, "w") as outf:
        for line in inf:
            if line[0] == ">":
                header = line[1:].strip().replace("\t", "    ")
                seq = inf.readline()
                outf.write(f"{header}\t{9606}\t{seq}")

if __name__ == "__main__":
    inPath = sys.argv[1]
    outPath = sys.argv[2]

    # get the file extension for the input file
    fileExtension = inPath.split(".")[-1]

    if fileExtension == "txt":
        convertTargetScanToFASTA(inPath, outPath)
    elif fileExtension in FASTA_FILE_EXTENSIONS:
        convertFASTAToTargetScan(inPath, outPath)
    else:
        print(f"Error: file extension {fileExtension} not supported.")
        sys.exit(1)
        
    print("done.")