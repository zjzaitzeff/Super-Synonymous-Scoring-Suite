import sys

FASTA_FILE_EXTENSIONS = ["fa", "fasta", "fna", "FASTA"]

def convertTargetScanToFASTA(inPath, outPath):
    with open(inPath) as inf, open(outPath, "w") as outf:
        for line in inf:
            # determine if there is a header line
            if line.startswith("miR\tfamily"):
                continue
            fullSequence = line.strip().split("\t")[4]
            outf.write(f">{line.replace('\t', '|')}")
            outf.write(f"{fullSequence}\n")

def convertFASTAToTargetScan(inPath, outPath):
    with open(inPath) as inf, open(outPath, "w") as outf:
        for line in inf:
            if line[0] == ">":
                line = line.replace("|", "\t")[1:]
                outf.write(line)

if __name__ == "__main__":
    # inPath = "../TargetScan/miR_Family_Info.txt"
    # outPath = "miR.fa"
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