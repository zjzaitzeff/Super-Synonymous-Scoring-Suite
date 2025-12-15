import sys
from os import path
import numpy as np

def parsePITAOutput(path):
    mRNAsToMiRNAs = {}
    mRNAsToMiRNAsToScores = {}
    with open(path) as inf:
        header = inf.readline()
        for line in inf:
            line = line.strip()
            lineItems = line.split("\t")
            mRNA = lineItems[0].split(" ")[0]
            # miRNA = lineItems[1].split("|")[0]
            # mRNA = lineItems[0] # TODO
            miRNA = lineItems[1] # TODO
            startPos = int(lineItems[2])
            endPos = int(lineItems[3])
            scores = lineItems[6:]
            scores = [float(score) for score in scores]
            miRNA = f"{miRNA}_{startPos}_{endPos}"
            if mRNA not in mRNAsToMiRNAs:
                mRNAsToMiRNAs[mRNA] = set()
            mRNAsToMiRNAs[mRNA].add(miRNA)
            if mRNA not in mRNAsToMiRNAsToScores:
                mRNAsToMiRNAsToScores[mRNA] = {}
            mRNAsToMiRNAsToScores[mRNA][miRNA] = scores
    return mRNAsToMiRNAs, mRNAsToMiRNAsToScores

def parseTargetScanOutput(path):
    mRNAsToMiRNAs = {}
    with open(path) as inf:
        header = inf.readline()
        for line in inf:
            line = line.strip()
            lineItems = line.split("\t")
            mRNA = lineItems[0].split(" ")[0]
            # mRNA = lineItems[0]
            miRNA = lineItems[1] # TODO this only gets the family- how to get full miRNA name?
            if mRNA not in mRNAsToMiRNAs:
                mRNAsToMiRNAs[mRNA] = set()
            mRNAsToMiRNAs[mRNA].add(miRNA)
    return mRNAsToMiRNAs

def parseSingleMiRandaScan(mRNAsToMiRNAs, mRNAsToMiRNAsToScores, inf, line):
    # get the two sequences compared
    if line.startswith("Performing Scan:"):
        miRNA, mRNA = line.strip().replace("Performing Scan: ", "").split(" vs ")
    elif line.startswith("Scores for this hit:"):
        line = inf.readline().strip()
        miRNA, mRNA = line.split("\t")[0:2]
    else:
        print(f"unexpected line in miRanda output: '{line}'")
        exit()
    # miRNA = miRNA.split("|")[0].lstrip(">") # TODO
    miRNA = miRNA.lstrip(">")
    if mRNA not in mRNAsToMiRNAs:
        mRNAsToMiRNAs[mRNA] = set()
    if mRNA not in mRNAsToMiRNAsToScores:
        mRNAsToMiRNAsToScores[mRNA] = {}
    
    line = inf.readline() # skip the =-=-=-= line (Performing Scan) or empty line (Scores for this hit)
    line = inf.readline().strip()
    if line == "": # this indicates a result was found: extract it!
        line = inf.readline().strip()
        if line.startswith("Forward:"):
            lineItems = line[line.index("R:")+2:].split(" ")
            startPos = int(lineItems[0])
            endPos = int(lineItems[2])
            score = float(line.split(" ")[1])
            miRNA = f"{miRNA}_{startPos}_{endPos}"
            mRNAsToMiRNAs[mRNA].add(miRNA)
            mRNAsToMiRNAsToScores[mRNA][miRNA] = score
        else:
            print(f"unexpected line in miRanda output: '{line}'")
            exit()

# TODO update this to use logic found here: https://www.biostars.org/p/251946/
# ie: grep -A 1 "Scores for this hit:" miranda_out.txt | sort | grep '>'
# which gives ex: >mirna_name      transcript_target        143.00  -22.87  2 18    339 357 16      75.00%  87.50%
# interpreted as: mirna Target  Score Energy-Kcal/Mol Query-Aln(start-end) Subjetct-Al(Start-End) Al-Len Subject-Identity Query-Identity
def parseMiRandaOutput(path):
    mRNAsToMiRNAs = {}
    mRNAsToMiRNAsToScores = {}
    with open(miRandaOutputPath) as inf:
        line = inf.readline()
        while line:
            if line.startswith("Performing Scan:") or line.startswith("Scores for this hit:"):
                parseSingleMiRandaScan(mRNAsToMiRNAs, mRNAsToMiRNAsToScores, inf, line)
            line = inf.readline()
    return mRNAsToMiRNAs, mRNAsToMiRNAsToScores

# def paccmitGetResultsFileName(path):
#     """
#     Gets the name of the most complete output file if there are results, otherwise returns None.
#     """
#     with open(path) as inf:
#         lines = inf.readlines()
#         if "no active genes left - terminating" in lines[-1]:
#             return None
#         else:
#             for line in reversed(lines):
#                 if "writing results to" in line:
#                     fileName = line.split(" ")[-1].strip().replace("'", "")
#                     return fileName

# def parsePACCMITOutput(path):
#     mRNAsToMiRNAs = {}
#     mRNAsToMiRNAsToScores = {}
#     with open(path, "r") as inf:
#         for line in inf:
#             line = line.strip()
#             lineItems = line.split("\t")
#             mRNA = lineItems[0]
#             miRNA = lineItems[-1]
#             if mRNA not in mRNAsToMiRNAs:
#                 mRNAsToMiRNAs[mRNA] = set()
#             mRNAsToMiRNAs[mRNA].add(miRNA)
#             if mRNA not in mRNAsToMiRNAsToScores:
#                 mRNAsToMiRNAsToScores[mRNA] = {}
#             score = float(lineItems[5])
#             mRNAsToMiRNAsToScores[mRNA][miRNA] = score
#     return mRNAsToMiRNAs, mRNAsToMiRNAsToScores

def interpretResults(inputFileOrderedmRNAs, mRNAsToMiRNAs, mRNAsToMiRNAsToScores, interpretedResultsPath):
    comparisonResults = {}
    uniqueMRNAs = list(mRNAsToMiRNAs.keys())
    # order the unique mRNAs based on the order of inputFileOrderedmRNAs
    uniqueMRNAs.sort(key=lambda x: inputFileOrderedmRNAs.index(x))
    # compare results between mRNAs
    comparisonsCompleted = set()
    allGainedTargets = set()
    allLostTargets = set()
    with open(interpretedResultsPath, "w") as outf:
        outf.write("comparison\ttotal unique miRNAs\tshared miRNAs\tdiffering miRNAs\tnum unique to mRNA1\tnum unique to mRNA2\tunique to mRNA1\tunique to mRNA2\n")
        if len(uniqueMRNAs) == 1:
            mRNA = uniqueMRNAs[0]
            outf.write(f"{mRNA} vs NONE\t{len(mRNAsToMiRNAs[mRNA])}\t0\t0\t{len(mRNAsToMiRNAs[mRNA])}\t0\t{sorted(mRNAsToMiRNAs[mRNA])}\tNone\n")
        else:
            for mRNA in uniqueMRNAs:
                otherMRNAs = set(uniqueMRNAs) - {mRNA}
                for othermRNA in otherMRNAs:
                    if (mRNA, othermRNA) in comparisonsCompleted or (othermRNA, mRNA) in comparisonsCompleted:
                        continue
                    comparisonsCompleted.add((mRNA, othermRNA))
                    miRNAsTargetingThis = mRNAsToMiRNAs[mRNA]
                    miRNAsTargetingOther = mRNAsToMiRNAs[othermRNA]
                    totalUniqueMiRNAs = miRNAsTargetingThis.union(miRNAsTargetingOther)
                    sharedMiRNAs = miRNAsTargetingThis.intersection(miRNAsTargetingOther)
                    totalDiffMiRNAs = miRNAsTargetingThis.symmetric_difference(miRNAsTargetingOther)
                    diffMiRNAs1 = miRNAsTargetingThis - miRNAsTargetingOther
                    diffMiRNAs2 = miRNAsTargetingOther - miRNAsTargetingThis
                    allGainedTargets.update(diffMiRNAs1)
                    allLostTargets.update(diffMiRNAs2)
                    comparisonResults[(mRNA, othermRNA)] = [diffMiRNAs1, diffMiRNAs2]
                    outf.write(f"{mRNA} vs {othermRNA}\t{len(totalUniqueMiRNAs)}\t{len(sharedMiRNAs)}\t{len(totalDiffMiRNAs)}\t{len(diffMiRNAs1)}\t{len(diffMiRNAs2)}\t{sorted(diffMiRNAs1)}\t{sorted(diffMiRNAs2)}\n")

            # which mRNAs have the same results?
            outf.write("\n")
            sameResults = {}
            for mRNA in uniqueMRNAs:
                otherMRNAs = set(uniqueMRNAs) - {mRNA}
                for othermRNA in sorted(otherMRNAs):
                    if mRNAsToMiRNAs[mRNA] == mRNAsToMiRNAs[othermRNA]:
                        if mRNA not in sameResults:
                            sameResults[mRNA] = set()
                        sameResults[mRNA].add(othermRNA)
            if sameResults:
                outf.write("mRNAs with the same miRNA interactions:\n")
                alreadyPrinted = set()
                for mRNA in sorted(sameResults):
                    if mRNA in alreadyPrinted:
                        continue
                    # add the mRNA and all the other mRNAs that have the same results
                    alreadyPrinted.add(mRNA)
                    for sameResult in sameResults[mRNA]:
                        alreadyPrinted.add(sameResult)
                    results = sorted([mRNA] + list(sameResults[mRNA]))
                    output = f"{results}"
                    outf.write(f"{output}\n")
            else:
                outf.write("No mRNAs have the same results.\n")

            # compare the scores between mRNAs that have overlapping miRNAs
            if mRNAsToMiRNAsToScores:
                outf.write("\n")
                for comparisonGroup in comparisonResults:
                    mRNA1, mRNA2 = comparisonGroup
                    if mRNA1 not in mRNAsToMiRNAsToScores or mRNA2 not in mRNAsToMiRNAsToScores:
                        continue
                    miRNAs1 = mRNAsToMiRNAs[mRNA1]
                    miRNAs2 = mRNAsToMiRNAs[mRNA2]
                    sharedMiRNAs = miRNAs1.intersection(miRNAs2)
                    if not sharedMiRNAs:
                        continue
                    outf.write(f"Compare different scores between shared miRNAs in {mRNA1} and {mRNA2}:\n")
                    scoresAllSame = True
                    for miRNA in sharedMiRNAs:
                        scores1 = mRNAsToMiRNAsToScores[mRNA1][miRNA]
                        scores2 = mRNAsToMiRNAsToScores[mRNA2][miRNA]
                        if scores1 == scores2:
                            continue
                        outf.write(f"\t{miRNA}: {scores1} - {scores2} = {np.array(scores1) - np.array(scores2)}\n")
                        scoresAllSame = False
                    if scoresAllSame:
                        outf.write("\tAll scores are the same.\n")

            # compare the comparison groups (ie: are the differences between wild and mutant1 the same as wild and mutant2?)
            if len(comparisonResults) > 1:
                comparisonGroupsThatAreTheSame = set()
                comparisonComparisonsCompleted = set()
                for comparisonGroup in sorted(comparisonResults):
                    diffsList = comparisonResults[comparisonGroup]
                    otherComparisonGroups = set(comparisonResults.keys()) - {comparisonGroup}
                    for otherComparisonGroup in sorted(otherComparisonGroups):
                        if (comparisonGroup, otherComparisonGroup) in comparisonComparisonsCompleted or (otherComparisonGroup, comparisonGroup) in comparisonComparisonsCompleted:
                            continue
                        comparisonComparisonsCompleted.add((comparisonGroup, otherComparisonGroup))
                        otherDiffsList = comparisonResults[otherComparisonGroup]
                        if diffsList[0] == otherDiffsList[0] and diffsList[1] == otherDiffsList[1]:
                            comparisonGroupsThatAreTheSame.add((comparisonGroup, otherComparisonGroup))
                outf.write("\n")
                if comparisonGroupsThatAreTheSame:
                    outf.write("Comparison groups that are the same:\n")
                    for comparisonGroupPair in sorted(comparisonGroupsThatAreTheSame):
                        output = f"'{comparisonGroupPair[0][0]} vs {comparisonGroupPair[0][1]}' == '{comparisonGroupPair[1][0]} vs {comparisonGroupPair[1][1]}'"
                        outf.write(f"{output}\n")
                else:
                    outf.write("No comparison groups are the same.\n")
    return allGainedTargets, allLostTargets, comparisonResults

if __name__ == "__main__":
    # TODO write this check better
    if len(sys.argv) != 3:
        print(f"ERROR: Invalid number of arguments.")
        print(f"Usage: python {sys.argv[0]} <mRNAFastaPath> <outputFolder>")
        exit(1)
    mRNAFastaPath = sys.argv[1]
    if not path.exists(mRNAFastaPath):
        print(f"ERROR: mRNA fasta file does not exist at: {mRNAFastaPath}")
        exit(1)
    outputFolder = sys.argv[2]

    # get the order of the mRNAs as they are in the input file
    inputFileOrderedmRNAs = []
    with open(mRNAFastaPath) as inf:
        for line in inf:
            line = line.strip()
            if line.startswith(">"):
                lineItems = line[1:].split("\t")
                mRNA = lineItems[0].split(" ")[0]
                inputFileOrderedmRNAs.append(mRNA)

    miRandaOutputPath = path.join(outputFolder, "miRanda_output.txt")
    if path.exists(miRandaOutputPath):
        miRandaInterpretedResultsPath = path.join(outputFolder, "PLX_miRanda_output_interpreted.txt")
        # TODO update this to include scores
        miRandamRNAsToMiRNAs, miRandaMRNAsToMiRNAsToScores = parseMiRandaOutput(miRandaOutputPath)
        miRandaAllGainedTargets, miRandaAllLostTargets, miRandaComparisonResults = interpretResults(inputFileOrderedmRNAs, miRandamRNAsToMiRNAs, miRandaMRNAsToMiRNAsToScores, miRandaInterpretedResultsPath)
        print("finished interpreting miRanda output")
    else:
        print(f"WARNING: skipping miRanda interpretation. No results at: {miRandaOutputPath}")

    pitaOutputPath = path.join(outputFolder, "_pita_results.tab")
    if path.exists(pitaOutputPath):
        pitaInterpretedResultsPath = path.join(outputFolder, "PLX_pita_output_interpreted.txt")
        pitamRNAsToMiRNAs, pitaMRNAsToMiRNAsToScores = parsePITAOutput(pitaOutputPath)
        pitaAllGainedTargets, pitaAllLostTargets, ptaComparisonResults = interpretResults(inputFileOrderedmRNAs, pitamRNAsToMiRNAs, pitaMRNAsToMiRNAsToScores, pitaInterpretedResultsPath)
        print("finished interpreting PITA output")
    else:
        print(f"WARNING: skipping PITA interpretation. No results at: {pitaOutputPath}")

    targetScanOutputPath = path.join(outputFolder, "targetscan_output.tsv")
    if path.exists(targetScanOutputPath):
        targetScanInterpretedResultsPath = path.join(outputFolder, "PLX_targetscan_output_interpreted.txt")
        targetScanmRNAsToMiRNAs = parseTargetScanOutput(targetScanOutputPath)
        targetScanAllGainedTargets, targetScanAllLostTargets, targetScanComparisonResults = interpretResults(inputFileOrderedmRNAs, targetScanmRNAsToMiRNAs, None, targetScanInterpretedResultsPath)
        print("finished interpreting TargetScan output")
    else:
        print(f"WARNING: skipping TargetScan interpretation. No results at: {targetScanOutputPath}")
    
    # paccmitLogPath = path.join(outputFolder, "paccmit.log")
    # if path.exists(paccmitLogPath):
    #     paccmitInterpretedResultsPath = path.join(outputFolder, "PLX_paccmit_output_interpreted.txt")
    #     # check the log file to see if there were any results
    #     resultsFileName = paccmitGetResultsFileName(paccmitLogPath)
    #     if resultsFileName:
    #         paccmitOutputPath = path.join(outputFolder, resultsFileName)
    #         paccmitmRNAsToMiRNAs, paccmitMRNAsToMiRNAsToScores = parsePACCMITOutput(paccmitOutputPath)
    #     else:
    #         paccmitmRNAsToMiRNAs = {}
    #         paccmitMRNAsToMiRNAsToScores = {}
    #     interpretResults(inputFileOrderedmRNAs, paccmitmRNAsToMiRNAs, paccmitMRNAsToMiRNAsToScores, paccmitInterpretedResultsPath)
    #     print("finished interpreting PACCMIT output")
    # else:
    #     print(f"WARNING: skipping PACCMIT interpretation. No results at: {paccmitLogPath}")
    
    # get prediction overlaps between the tools
    # TODO add shared score directional changes too
    overlapsOutputPath = path.join(outputFolder, "PLX_miRNA_prediction_overlaps.txt")
    with open(overlapsOutputPath, "w") as outf:

        # truncate all miRNAs to the first bar
        if miRandamRNAsToMiRNAs:
            miRandaAllGainedTargets = set([miRNA.split("|")[0] for miRNA in miRandaAllGainedTargets])
            miRandaAllLostTargets = set([miRNA.split("|")[0] for miRNA in miRandaAllLostTargets])
        if pitamRNAsToMiRNAs:
            pitaAllGainedTargets = set([miRNA.split("|")[0] for miRNA in pitaAllGainedTargets])
            pitaAllLostTargets = set([miRNA.split("|")[0] for miRNA in pitaAllLostTargets])

        if targetScanmRNAsToMiRNAs and miRandamRNAsToMiRNAs:
            sharedGainedTargets = targetScanAllGainedTargets.intersection(miRandaAllGainedTargets)
            sharedGainedTargets = "|".join(sharedGainedTargets) if sharedGainedTargets else ""
            sharedLostTargets = targetScanAllLostTargets.intersection(miRandaAllLostTargets)
            sharedLostTargets = "|".join(sharedLostTargets) if sharedLostTargets else ""
            outf.write(f"Shared gained targets between TargetScan and miRanda: {sharedGainedTargets}\n")
            outf.write(f"Shared lost targets between TargetScan and miRanda: {sharedLostTargets}\n") 
        if targetScanmRNAsToMiRNAs and pitamRNAsToMiRNAs:
            sharedGainedTargets = targetScanAllGainedTargets.intersection(pitaAllGainedTargets)
            sharedGainedTargets = "|".join(sharedGainedTargets) if sharedGainedTargets else ""
            sharedLostTargets = targetScanAllLostTargets.intersection(pitaAllLostTargets)
            sharedLostTargets = "|".join(sharedLostTargets) if sharedLostTargets else ""
            outf.write(f"Shared gained targets between TargetScan and PITA: {sharedGainedTargets}\n")
            outf.write(f"Shared lost targets between TargetScan and PITA: {sharedLostTargets}\n")
        if miRandamRNAsToMiRNAs and pitamRNAsToMiRNAs:
            sharedGainedTargets = miRandaAllGainedTargets.intersection(pitaAllGainedTargets)
            sharedGainedTargets = "|".join(sharedGainedTargets) if sharedGainedTargets else ""
            sharedLostTargets = miRandaAllLostTargets.intersection(pitaAllLostTargets)
            sharedLostTargets = "|".join(sharedLostTargets) if sharedLostTargets else ""
            outf.write(f"Shared gained targets between PITA and miRanda: {sharedGainedTargets}\n")
            outf.write(f"Shared lost targets between PITA and miRanda: {sharedLostTargets}\n")
        if targetScanmRNAsToMiRNAs and miRandamRNAsToMiRNAs and pitamRNAsToMiRNAs:
            sharedGainedTargets = targetScanAllGainedTargets.intersection(miRandaAllGainedTargets).intersection(pitaAllGainedTargets)
            sharedGainedTargets = "|".join(sharedGainedTargets) if sharedGainedTargets else ""
            sharedLostTargets = targetScanAllLostTargets.intersection(miRandaAllLostTargets).intersection(pitaAllLostTargets)
            sharedLostTargets = "|".join(sharedLostTargets) if sharedLostTargets else ""
            outf.write(f"Shared gained targets between TargetScan, PITA, and miRanda: {sharedGainedTargets}\n")
            outf.write(f"Shared lost targets between TargetScan, PITA, and miRanda: {sharedLostTargets}\n")

    print("done.")