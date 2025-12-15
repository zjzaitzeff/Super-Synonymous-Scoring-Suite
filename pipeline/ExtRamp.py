#! /usr/bin/env python3
"""
This program is intended to extract a ramp of slowly translated codons from the beginning
of a gene sequence (DNA or RNA).
"""

from statistics import mean, median, stdev
import numpy as np
from scipy.stats import gmean, hmean, zscore
from math import exp, log, sqrt, ceil
import gzip
import argparse
import sys
import re
from multiprocessing import Pool, freeze_support
import time
import os
import traceback

VERSION = "2.0 beta 21"

def makeArgParser():
    parser = argparse.ArgumentParser(description="Extract the individual Ramp sequences from a collection of genes")
    parser.add_argument("-i", "--input", type=str, required=True, help="(input) Required FASTA file containing gene sequences (.gz or .gzip for gzipped)")
    parser.add_argument("-a", "--tAI", type=str, help="(input) csv file containing the species tAI values")
    parser.add_argument("-u", "--rscu", type=str, help="(input) FASTA file used to compute relative synonymous codon usage and relative adaptiveness of codons")
    parser.add_argument("-o", "--ramp", type=str, help="(output) FASTA file to write ramp sequences to")
    parser.add_argument("-y", "--score", default = None, type=str, help="(output) tsv file to write ramp strength and robustness scores.")
    parser.add_argument("-j", "--wij", default = None, type=str, help="(output) csv file to write calculated codon wij (efficiency) values. Note that tAI and wij cannot be used at the same time or the wij file would be a copy of the tAI file.")
    parser.add_argument("-l", "--vals", type=str, help="(output) FASTA-like file to write tAI/relative adaptiveness values for ribosome-smoothed efficiency at each window")
    parser.add_argument("-p", "--speeds", type=str, help="(output) speeds file to write tAI/relative adaptiveness values for each position in the sequence from the codon after the start codon to the codon before the stop codon. Format: Header newline list of values")
    parser.add_argument("-n", "--noRamp", type=str, help="(output) txt file to write the gene names that contained no ramp sequence")
    parser.add_argument("-z", "--removedSequences", default = None, type=str, help="(output) Write the header lines that are removed (e.g., sequence not long enough or not divisible by 3) to output file")
    parser.add_argument("-x", "--afterRamp", type=str, required=False, help="(output) FASTA file containing gene sequences after the identified ramp sequence")
    parser.add_argument("-t", "--threads", type=int, help="the number of threads you want to run, default = all")
    parser.add_argument("-w", "--window", type=int, default = 9, help="the ribosome window size in codons, default = 9 codons")
    parser.add_argument("-s", "--stdev", type=float, default = -1.0, help="the number of standard deviations below the mean the cutoff value to be included as a ramp for gmean and mean. Not used by default")
    parser.add_argument("-d", "--stdevRampLength", type=float, default = -1.0, help="the number of standard deviations used in the quality check step (lengths of the ramp sequences). Not done by default")
    parser.add_argument("-m", "--middle", type=str, default = "hmean", help="the type of statistic used to measure the middle (consensus) efficiency. Options: hmean, gmean, mean, median")
    parser.add_argument("-v", "--verbose", action="store_true", help="flag to print progress to standard error")
    parser.add_argument("-r", "--rna", action="store_true", help="flag for RNA sequences. Default: DNA")
    parser.add_argument("-f", "--determine_cutoff", action="store_true", help="flag to determine outlier percentages for mean cutoff based on species FASTA file. Default: local minimum in first 8 percent of gene")
    parser.add_argument("-c", "--cutoff", type=int, default = 8, help="Cutoff for where the local minimum must occur in the gene for a ramp to be calculated. If --determine_cutoff (-f) is used, then this value may change. Is not used if standard deviations are set.")
    parser.add_argument("-e", "--determine_cutoff_percent", type=str, default = "outlier", help="Cutoff percentile for determining percent of gene that is in an outlier region. Used in conjunction with -f. Default is true outliers. Other options include numbers from 0-99, which indicate the region of a box plot. For instance, 75 means the 75th percentile or above.")
    parser.add_argument("-q", "--seqLength", type=float, default = 300, help="Minimum nucleotide sequence length. Default is 300 nucleotides (100 amino acids)")
    parser.add_argument("-g", "--prevalidate", action="store_true", help="flag to truncate at early stop codons before testing if sequences are valid")
    args = parser.parse_args()

    if args.determine_cutoff_percent != "outlier" and not args.determine_cutoff_percent.isdigit() or (args.determine_cutoff_percent.isdigit() and (int(args.determine_cutoff_percent) < 0 or int(args.determine_cutoff_percent) > 99)):
        sys.stderr.write("Warning: determine_cutoff_percent (-e) must either be 'outlier' or an integer from 0-99. It will be set to 'outlier'.\n")
        args.determine_cutoff_percent = "outlier"
    if args.determine_cutoff != True and args.determine_cutoff_percent != "outlier":
        sys.stderr.write("Warning: determine_cutoff_percent (-e) was specified without the determine_cutoff (-f) flag. The -f flag will be set to True.\n")
        args.determine_cutoff = True
    if args.cutoff < 0 or args.cutoff > 99:
        sys.stderr.write("Warning: cutoff (-c) must be an integer from 0-99. It will be set to 8.\n")
        args.cutoff = 8
    if args.threads != None and args.threads < 0:
        sys.stderr.write("Warning: threads (-t) must be a positive number. It will be set to 'all'.\n")
        args.threads = None
    if args.window < 0:
        sys.stderr.write("Warning: window (-w) must be a positive number. It will be set to 9.\n")
        args.window = 9
    if args.stdev < 0 and args.stdev != -1:
        sys.stderr.write("Warning: stdev (-s) must be a positive number. It will be ignored.\n")
        args.stdev = -1.0
    if args.stdevRampLength < 0 and args.stdevRampLength != -1:
        sys.stderr.write("Warning: stdevRampLength (-d) must be a positive number. It will be ignored.\n")
        args.stdevRampLength = -1.0
    if args.middle == "hmean" and args.stdev >= 0:
        sys.stderr.write("Warning: hmean cannot be used with stdev (-s). -s will be ignored.\n")
        args.stdev = -1.0
    if args.stdev >= 0:
        if args.cutoff != 8:
            sys.stderr.write("Warning: cutoff (-c) cannot be used with stdev (-s). -c will be ignored.\n")
            args.cutoff = 8
        if args.determine_cutoff == True:
            sys.stderr.write("Warning: determine_cutoff (-f) cannot be used with stdev (-s). -f will be ignored.\n")
            args.determine_cutoff = None
        if args.determine_cutoff_percent != "outlier":
            sys.stderr.write("Warning: determine_cutoff_percent (-e) cannot be used with stdev (-s). -e will be ignored.\n")
            args.determine_cutoff_percent = "outlier"
    if args.determine_cutoff == True and args.cutoff != 8:
        sys.stderr.write("Warning: cutoff (-c) may change because determine_cutoff (-f) is set.\n")
    if (args.window * 3) > args.seqLength:
        args.seqLength = args.window * 3 + 1
        sys.stderr.write(f"Warning: to accommodate the window size (-w) of {args.window}, the minimum sequence length, seqLength (-q), has been set to {args.seqLength}\n")
    if args.wij is not None and args.tAI is not None:
        sys.stderr.write("Warning: wij (-j) cannot be used with tAI (-a) because the resulting wij file would be a copy of the tAI file. -j will be ignored.\n")
        args.wij = None
    return args

def fileExistsCheck(filename):
    """
    Input: a file path.
    Exits the program if the file path does not exist.
    """
    if not os.path.isfile(filename):
        sys.stderr.write(f"Error: '{filename}' does not exist\nExiting...\n")
        sys.exit()

def readSeqFile(args, inputFile):
    """
    Input: System arguments and an input file path.
    Returns an array of tuples (Name of sequence, Sequence) created from the input FASTA file.
    Sequences which are not divisible by three or have non-standard nucleotide characters are removed.
    Also returns a message to be printed if the program is in verbose mode.
    """
    fileExistsCheck(inputFile)

    seqTupleArray = []
    message = ""
    curSeqName = ""
    curSeq = ""
    inFile = ""
    countInvalidCharSeq = 0
    countInvalidDivSeq = 0
    countTooShortSeq = 0
    countTruncated = 0
    countSavedTruncated = 0
    containsU = False
    if inputFile.endswith(".gz") or inputFile.endswith(".gzip"):
        inFile = gzip.open(inputFile, "rt")
    else:
        inFile = open(inputFile, "r")
    badSeq = 0
    removedSeq = ""
    if args.removedSequences and inputFile == args.input:
        removedSeq = open(args.removedSequences, "w")

    for line in inFile:
        if line[0] == ">":
            if curSeq != "":
                seqTupleArray, curSeq, badSeq, countInvalidCharSeq, countInvalidDivSeq, countTooShortSeq, countTruncated, countSavedTruncated, containsU = readSeqFileLineHelper(seqTupleArray, curSeqName, curSeq, inputFile, removedSeq, badSeq, countInvalidCharSeq, countInvalidDivSeq, countTooShortSeq, countTruncated, countSavedTruncated, containsU)
            curSeqName = line.strip()
            curSeq = ""
        else:
            curSeq += line.strip().upper()
    seqTupleArray, curSeq, badSeq, countInvalidCharSeq, countInvalidDivSeq, countTooShortSeq, countTruncated, countSavedTruncated, containsU = readSeqFileLineHelper(seqTupleArray, curSeqName, curSeq, inputFile, removedSeq, badSeq, countInvalidCharSeq, countInvalidDivSeq, countTooShortSeq, countTruncated, countSavedTruncated, containsU)
    if args.verbose and inputFile == args.input:
        if countTruncated > 0:
            message += f"\n{countTruncated} sequences were truncated at early stop codons."
            message += f"\n\t{countSavedTruncated} sequences were saved through truncation."
        if badSeq > 0:
            message += f"\n{badSeq} sequences were invalid so they were removed!\n"
            if countInvalidCharSeq > 0:
                message += f"\t{countInvalidCharSeq} had non-standard nucleotides.\n"
            if countInvalidDivSeq > 0:
                message += f"\t{countInvalidDivSeq} were not divisible by three.\n"
            if countTooShortSeq > 0:
                message += f"\t{countTooShortSeq} were shorter than the {args.seqLength} nucleotide cutoff.\n"
    inFile.close()
    if args.removedSequences and inputFile == args.input:
        removedSeq.close()
    return seqTupleArray, message

def readSeqFileLineHelper(seqTupleArray, curSeqName, curSeq, inputFile, removedSeq, badSeq, countInvalidCharSeq, countInvalidDivSeq, countTooShortSeq, countTruncated, countSavedTruncated, containsU):
    """
    Used by the readSeqFile function to read a single sequence line from a CDS file (assumes the sequence is all on one line).
    Returns the modified seqTupleArray, the potentially modified sequence, and the updated badSeq, invalidCharSeq, invalidDivSeq, tooShortSeq, and truncated counts.
    """
    if args.rna:
        curSeq.replace("U", "T")
    chopped = False
    if args.prevalidate:
        oldSeq = curSeq
        curSeq = cutAtEarlyStop(curSeq)
        if len(oldSeq) != len(curSeq):
            chopped = True
    match = re.search("[^ATCG]", curSeq)
    invalid = False
    if match is not None:
        if "U" in curSeq and not containsU:
            containsU = True
            sys.stderr.write(f"Warning: file contains U. Did you mean to use the --rna flag?\n")
        countInvalidCharSeq += 1
        invalid = True
    elif len(curSeq) % 3 != 0:
        countInvalidDivSeq += 1
        invalid = True
    elif len(curSeq) < args.seqLength:
        countTooShortSeq += 1
        invalid = True

    if chopped:
        countTruncated += 1
    if not invalid:
        if chopped:
            countSavedTruncated += 1
        startCodon = curSeq[0:3]
        curSeq = curSeq[3:] # remove start codon
        stopCodon = curSeq[-3:]
        curSeq = curSeq[0:-3] # remove stop codon
        seqTupleArray.append((curSeqName, curSeq, startCodon, stopCodon))
    else:
        if args.removedSequences and inputFile == args.input:
            removedSeq.write(curSeqName +"\n")
        badSeq += 1

    return seqTupleArray, curSeq, badSeq, countInvalidCharSeq, countInvalidDivSeq, countTooShortSeq, countTruncated, countSavedTruncated, containsU

def cutAtEarlyStop(seq):
    """
    Returns the sequence after truncating everything after the first stop codon.
    """
    stopCodons = ["TAA", "TAG", "TGA"]
    for index in range(0, len(seq), 3):
        if seq[index:index+3] in stopCodons:
            break
    return seq[0:index+3]

def calcCodonSpeeds(seqTupleArray):
    """
    Input: an array of tuples (header, sequence) of genes.
    Returns a dictionary of codons to relative translation Speed.
    It calculates the codon proportions using all of the gene sequences and
    then uses those values to determine translation speed. Uses metrics from CAI values
    Relative Synonymous Codon Usage (RSCU) and
    Relative adaptiveness of codon (w) are calculated based on equations presented in
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC340524/pdf/nar00247-0410.pdf
    """
    codonToSpeed = {}
    totalCodons = 0
    codonGroups = AAtoCodonsDict()
    for elem in seqTupleArray:
        codonToSpeed, totalCodons = countCodons(elem[1], codonToSpeed, totalCodons)
    for aa in codonGroups:
        totalCounts = 0
        numCodons = 0
        for codon in aa:
            numCodons += 1
            if not codon in codonToSpeed:
                codonToSpeed[codon] = 0.0001 # codonToSpeed has count of codon usages
                continue
            totalCounts += codonToSpeed[codon]
        maxRSCU = 0.0
        for codon in aa:
            if numCodons == 0 or totalCounts == 0:
                rscu = 0.0001
                codonToSpeed[codon] = rscu # codonToSpeed has RSCU
                continue
            rscu = codonToSpeed[codon] / ((1.0 / numCodons) * totalCounts)
            if rscu > maxRSCU:
                maxRSCU = rscu
            codonToSpeed[codon] = rscu # codonToSpeed has RSCU
        for codon in aa:
            if maxRSCU == 0:
                w = 0.0001
                codonToSpeed[codon] = w # codonToSpeed has relative adaptiveness of codons
                continue
            w = codonToSpeed[codon] / maxRSCU # w=relative adaptiveness of codons
            codonToSpeed[codon] = w # codonToSpeed has relative adaptiveness of codons
    return codonToSpeed

def csvToDict(filename):
    """
    Input: path to a csv/tsv/space-separated file with tAI values where the first row is a list of codons
    and the second row is the list of tAI values for those codons. The file could
    also have a codon and its tAI value on the same line, each codon separated by a new line.
    Returns a dictionary of uppercase DNA codons to their tAI speed from the user provided tAI file.
    """
    fileExistsCheck(filename)

    tempDict = {}
    delimiter = re.compile(r"[,\s\t]+")
    try:
        with open(filename) as csvfile:
            for row in csvfile:
                row = row.strip().replace("\"", "").replace("'", "")
                codon, tAI = delimiter.split(row)
                if float(tAI) == 0:
                    tAI = 0.0001
                tempDict[codon.upper().replace("U", "T")] = float(tAI)
    except:
        try:
            with open(filename) as csvfile:
                count = 0
                codonArray = []
                for row in csvfile:
                    row = row.strip().replace("\"", "").replace("'", "")
                    row = delimiter.split(row)
                    if count % 2 == 0:
                        codonArray = row
                    else:
                        for i in range(len(row)):
                            if float(row[i]) == 0:
                                tempDict[codonArray[i].upper().replace("U", "T")] = 0.0001
                            else:
                                tempDict[codonArray[i].upper().replace("U", "T")] = float(row[i])
                    count += 1
        except:
            sys.stderr.write(f"Error: tAI file is not formatted correctly\nExiting...\n")
            sys.exit()
    if len(tempDict) < 61:
        sys.stderr.write("Error: tAI file is missing codon(s)\nExiting...\n")
        sys.exit()
    return tempDict

def AAtoCodonsDict():
    """
    Returns a list of a list of codons that encode each amino acid.
    Groupings are in the following order (one letter code):
    I,L,V,F,M,C,A,G,P,T,S,Y,W,Q,N,H,E,D,K,R,
    """
    return [["ATT","ATC","ATA"],["CTT","CTC","CTA","CTG","TTA","TTG"],["GTT","GTC","GTA","GTG"],["TTT","TTC"],["ATG"],["TGT","TGC"],["GCT","GCC","GCA","GCG"],["GGT","GGC","GGA","GGG"],["CCT","CCC","CCA","CCG"],["ACT","ACC","ACA","ACG"],["TCT","TCC","TCA","TCG","AGT","AGC"],["TAT","TAC"],["TGG"],["CAA","CAG"],["AAT","AAC"],["CAT","CAC"],["GAA","GAG"],["GAT","GAC"],["AAA","AAG"],["CGT","CGC","CGA","CGG","AGA","AGG"]]

def countCodons(sequence, codonToSpeed, totalCodons):
    """
    Input: A DNA sequence.
    Returns a dictionary of codons to the total number of times the codon was found
    in all provided gene sequences. These are later used to find the proportions of 
    codon usage (estimate speed of translation). Also returns the total number of codons.
    """
    for i in range(0, len(sequence), 3):
        if sequence[i:i+3] not in codonToSpeed:
            codonToSpeed[sequence[i:i+3]] = 0
        totalCodons += 1
        codonToSpeed[sequence[i:i+3]] += 1
    return codonToSpeed, totalCodons

def calcSeqSpeed(elem):
    """
    Maps the gene name to an array of the tAI/speed values for the given sequence.
    Input: tuple (header, sequence) of a FASTA record
    Output: A dictionary of seq name to a tuple (translational speed of each sequence, cutOffValue for outlier sequence)
    """
    seqToSpeed = {}
    seq = elem[1]
    speedArray = []
    stopCodons = ["TAA", "TAG", "TGA"]
    for index in range(0, len(seq), 3):
        if not seq[index:index+3] in stopCodons:
            speedArray.append(codonToSpeed[seq[index:index+3]])
        else:
            speedArray.append(0.0001)
    seqToSpeed[elem[0]] = np.array(speedArray)
    return seqToSpeed

def createSpeedsDict(result):
    """
    Returns a dictionary of all the sequence names to the array of their tAI/speed values.
    """
    seqToSpeed = {}
    for elem in result:
        seqToSpeed.update(elem)
    return seqToSpeed

def findStDev(array, mean):
    """
    Calculates geometric standard deviation from array and geometric mean.
    """
    if args.middle != "gmean": # if the user specified other mean
        return stdev(array, mean)
    # for gMean
    summation = 0
    for item in array:
        summation += (log(item/mean))**2
    return exp(sqrt(summation/len(array)))

def hmeanSliding(a, ribosomeWindowLength):
    """
    Returns the sliding harmonic mean of the array (a), given the window length (ribosomeWindowLength)
    using a convolution.
    """
    a = 1 / a
    a = np.convolve(a, np.ones(ribosomeWindowLength)/ribosomeWindowLength, mode="valid")
    return 1 / a

def gmeanSliding(a, ribosomeWindowLength):
    """
    Returns the sliding geometric mean of the array (a), given the window length (ribosomeWindowLength)
    using a convolution.
    """
    a = np.log(a)
    a = np.convolve(a, np.ones(ribosomeWindowLength)/ribosomeWindowLength, mode="valid")
    return np.exp(a)

def meanSliding(a, ribosomeWindowLength):
    """
    Returns the sliding arithmetic mean of the array (a), given the window length (ribosomeWindowLength)
    using a convolution.
    """
    return np.convolve(a, np.ones(ribosomeWindowLength)/ribosomeWindowLength, mode="valid")

def medianSliding(a, ribosomeWindowLength):
    """
    Returns the strided median of the array (a), given the window length (ribosomeWindowLength)
    using np.lib.stride.
    """
    rolling_window = np.lib.stride_tricks.as_strided(a, shape=(a.size - ribosomeWindowLength + 1, ribosomeWindowLength), strides=a.strides * 2, writeable=False)
    return np.median(rolling_window, axis=1)

def getSlidingMean(middle, a, ribosomeWindowLength):
    """
    Gets the sliding window mean of the array "a," given a middle function string name and the window length.
    """
    if middle=="mean":
        windowMeans = meanSliding(a, ribosomeWindowLength)
    elif middle=="hmean":
        windowMeans = hmeanSliding(a, ribosomeWindowLength)
    elif middle=="gmean":
        windowMeans = gmeanSliding(a, ribosomeWindowLength)
    elif middle=="median":
        windowMeans = medianSliding(a, ribosomeWindowLength)
    return windowMeans

def findWindowMeansThreaded(elem):
    """
    Calculates the window mean speeds for the specified element of seqToSpeed and returns it in tuple format.
    """
    a = np.array(seqToSpeed[elem])
    windowMeans = getSlidingMean(args.middle, a, ribosomeWindowLength)
    return (elem, windowMeans)

def calcRiboSpeed(speeds):
    """
    Returns the average of the tAI/speed values in the ribosomeWindow.
    """
    if len(speeds) == 0:
        return 0.0001
    elif args.middle == "hmean":
        return hmean(speeds)
    elif args.middle == "gmean":
        return gmean(speeds)
    elif args.middle == "mean":
        return mean(speeds)
    elif args.middle == "median":
        return median(speeds)

def writeEfficiencyFile(outPath, efficiencyVals):
    """
    Writes the efficiency values (sequence codon efficiencies or window means) to a FASTA-like file.
    """
    with open(outPath, "w") as outf:
        for header in efficiencyVals:
            outf.write(header + "\n" + ",".join(map(str, efficiencyVals[header])) + "\n")

def outputRampSeqs(rampSeqs, args):
    """
    Prints the ramp sequences to the terminal or a FASTA file
    Input: tuple of header, sequence
    """
    if args.stdevRampLength >= 0:
        rampSeqs = qualityCheck(rampSeqs, args.stdevRampLength)
    outPut = sys.stdout
    if args.ramp:
        outPut = open(args.ramp, "w")
    if args.afterRamp:
        afterRampFile = open(args.afterRamp, "w")
    if args.noRamp:
        noRampFile = open(args.noRamp, "w")
    count = 0
    for line in rampSeqs:
        if not args.afterRamp:
            if not line.startswith("None"):
                count += 1
                outPut.write(line)
            elif args.noRamp:
                noRampFile.write(line[4:])
        else:
            if not line[0].startswith("None"):
                count += 1
                outPut.write(line[0])
                afterRampFile.write(line[1])
            elif args.noRamp:
                noRampFile.write(line[0][4:])
    outPut.close()
    if args.verbose:
        sys.stderr.write(f"\n{count} Ramp Sequences found\n")
    if args.noRamp:
        noRampFile.close()

def isolateRamp(vals):
    """
    Used when stdev is used. Calculates the cutoff point for the individual ramp sequences and returns them.
    """
    record, percentThatIsRamp = vals
    mean = middleFunc[args.middle](seqToSpeed[record[0]])
    stdev = findStDev(seqToSpeed[record[0]], mean)
    cutOffVal = mean - (args.stdev*stdev)
    i = 0
    while calcRiboSpeed(seqToSpeed[record[0]][i:i+ribosomeWindowLength]) < cutOffVal and (i+ribosomeWindowLength) <= len(seqToSpeed[record[0]]):
        i += 1
    if i == 0:
        if not args.afterRamp:
            return "None" + record[0] + "\n"
        else:
            return tuple(["None" + record[0] + "\n"])

    if not args.afterRamp:
        return record[0] + "\n" + record[2] + record[1][:(i+ribosomeWindowLength)*3] + "\n"
    else:
        return tuple([record[0] + "\n" + record[2] + record[1][:(i+ribosomeWindowLength)*3] + "\n", record[0] + "\n" + record[1][(i+ribosomeWindowLength)*3:] +record[3]+ "\n"])

def isolateRampHmean(vals):
    """
    Used when stdev is not used. Calculates the cut off point for the individual ramp sequences and returns them.
    """
    record, percentThatIsRamp = vals
    speeds = seqToSpeed[record[0]]
    mean = hmean(speeds)
    windowMeans = getSlidingMean(args.middle, speeds, ribosomeWindowLength)
    bestMean = min(windowMeans)
    pos = np.where(windowMeans == bestMean)[0][0] # gets the first minimum
    perc = float(100*(pos/float(len(windowMeans)))) # percentages x 100
    if perc <= percentThatIsRamp:
        lastRampWindow = getLastRampWindowIndex(pos, windowMeans, mean)
        if lastRampWindow != -1:
            if not args.afterRamp:
                return record[0] + "\n" + record[2] + record[1][:(lastRampWindow+ribosomeWindowLength)*3] + "\n"
            else:
                return tuple([record[0] + "\n" + record[2] + record[1][:(lastRampWindow+ribosomeWindowLength)*3] + "\n", record[0] + "\n" + record[1][(lastRampWindow+ribosomeWindowLength)*3:] +record[3]+ "\n"])
    if not args.afterRamp:
        return "None" + record[0] + "\n"
    else:
        return tuple(["None" + record[0] + "\n"])

def getLastRampWindowIndex(pos, windowMeans, mean):
    """
    Gets the index of the last efficiency window of the ramp, which is the index of the first value
    after the ramp starts that is greater than the sequence's mean efficiency. If no such index exists,
    the ramp does not exist and -1 is returned.
    """
    indicesGreaterThanMean = np.where(windowMeans[pos:] >= mean) + pos
    return indicesGreaterThanMean[0][0] if len(indicesGreaterThanMean[0]) > 0 else -1

def qualityCheck(rampSeqs, numStDev):
    """
    Looks at the lengths of the ramp sequences and only use those that are within numStDev standard
    deviations of the mean length. The distribution is right-skewed so the data is log transformed.
    Returns the list of rampSeqs with the extremes removed.
    """
    lengths = []
    tempSeqs = []
    for i in range(len(rampSeqs)):
        if not rampSeqs[i].startswith("None"):
            ramp = rampSeqs[i].split("\n")
            lengths.append(log(len(ramp[1])))
            tempSeqs.append(i)
    if len(lengths) == 0:
        if args.verbose:
            sys.stderr.write("NO RAMP SEQUENCES FOUND\n")
        sys.exit()
    if len(lengths) == 1:
        return rampSeqs

    meanLen = mean(lengths)
    std = stdev(lengths, meanLen)
    removed = 0
    for index in tempSeqs:
        ramp = rampSeqs[index].split("\n")
        if log(len(ramp[1])) < meanLen - (numStDev * std) or log(len(ramp[1])) > meanLen + (numStDev * std):
            newLine = "None" + ramp[0] + "\n"
            rampSeqs[index] = newLine
            removed +=1
    if args.verbose:
        sys.stderr.write(f"{removed} sequences removed by standard deviation check of all ramp sequences (-d option).\n")
    return rampSeqs

def getBottleneck(header):
    """
    Input: header line of a sequence.
    Output: a list of tuples (percentage location of window mean minimum in gene, header, position of minimum).
    """
    allPercents = []
    speeds = seqToSpeed[header]
    windowMeans = []
    windowMeans = getSlidingMean(args.middle, speeds, args.window)
    bestMean = min(windowMeans)
    pos = [index for index, value in enumerate(windowMeans) if value == bestMean] # ensures that all local minima are included
    for p in pos:
        perc = float(ceil(100.0*((p+1)/float(len(windowMeans))))) # percentages x 100
        allPercents.append((perc, header, p))
    return allPercents

def getCutoffValue(counts):
    """
    Input: a dictionary of header lines -> tuple array of codon efficiency values for each codon.
    Output: Locations that have more local minimum codon efficiencies and are outliers, starting from the front of the gene. If no outliers exist, return 0.
    Also returns a message to be printed if the program is in verbose mode.
    """
    message = ""
    instances = list(counts.values())
    instances = np.array(instances)
    upper = 0
    if args.determine_cutoff_percent.isdigit():
        upper = np.percentile(instances, int(args.determine_cutoff_percent))
    else:
        q1= np.percentile(instances, 25)
        q3= np.percentile(instances, 75)
        upper = q3 + (1.5*(q3-q1))
    outliers = set()
    for percent in instances:
        if percent >= upper:
            outliers.add(percent)
    if args.verbose:
        message = f"The following percentages exceed the threshold ({args.determine_cutoff_percent}):"
        for x in range(1, 101):
            if counts[x] in outliers:
                message += f" {x}"

    for x in range(1, 101):
        if counts[x] in outliers:
            continue
        return x-1, message # -1 because it goes 1 past the last outlier. Divide by 100 to turn percent into decimal
    message += f"\nWarning: the cutoff percentage has been set to 0!"
    return 0, message

def init_pool(codonToSpeedParam):
    """
    Necessary for Windows multiprocessing.
    Used for calculating sequence speeds.
    """
    global codonToSpeed
    codonToSpeed = codonToSpeedParam

def init_pool2(seqToSpeedParam, ribosomeWindowLengthParam, middleFuncParam, argsParam):
    """
    Necessary for Windows multiprocessing.
    Used for calculating ramp cutoff, ramp sequences, and window means.
    """
    global seqToSpeed
    global ribosomeWindowLength
    global middleFunc
    global args
    seqToSpeed = seqToSpeedParam
    ribosomeWindowLength = ribosomeWindowLengthParam
    middleFunc = middleFuncParam
    args = argsParam

def init_pool3(codonToSpeedParam, seqToSpeedParam, percentThatIsRampParam, ribosomeWindowLengthParam, middleFuncParam, argsParam):
    """
    Necessary for Windows multiprocessing.
    Used for calculating ramp scores.
    """
    global codonToSpeed
    global seqToSpeed
    global codonToSlowerValidMutations
    global codonToFasterValidMutations
    global percentThatIsRamp
    global ribosomeWindowLength
    global middleFunc
    global args
    codonToSpeed = codonToSpeedParam
    seqToSpeed = seqToSpeedParam
    codonToValidMutations = {"AAA": ["AAC", "AAG", "AAT", "ACA", "AGA", "ATA", "CAA", "GAA"], "AAC": ["AAA", "AAG", "AAT", "ACC", "AGC", "ATC", "CAC", "GAC", "TAC"], "AAG": ["AAA", "AAC", "AAT", "ACG", "AGG", "ATG", "CAG", "GAG"], "AAT": ["AAA", "AAC", "AAG", "ACT", "AGT", "ATT", "CAT", "GAT", "TAT"], "ACA": ["AAA", "ACC", "ACG", "ACT", "AGA", "ATA", "CCA", "GCA", "TCA"], "ACC": ["AAC", "ACA", "ACG", "ACT", "AGC", "ATC", "CCC", "GCC", "TCC"], "ACG": ["AAG", "ACA", "ACC", "ACT", "AGG", "ATG", "CCG", "GCG", "TCG"], "ACT": ["AAT", "ACA", "ACC", "ACG", "AGT", "ATT", "CCT", "GCT", "TCT"], "AGA": ["AAA", "ACA", "AGC", "AGG", "AGT", "ATA", "CGA", "GGA"], "AGC": ["AAC", "ACC", "AGA", "AGG", "AGT", "ATC", "CGC", "GGC", "TGC"], "AGG": ["AAG", "ACG", "AGA", "AGC", "AGT", "ATG", "CGG", "GGG", "TGG"], "AGT": ["AAT", "ACT", "AGA", "AGC", "AGG", "ATT", "CGT", "GGT", "TGT"], "ATA": ["AAA", "ACA", "AGA", "ATC", "ATG", "ATT", "CTA", "GTA", "TTA"], "ATC": ["AAC", "ACC", "AGC", "ATA", "ATG", "ATT", "CTC", "GTC", "TTC"], "ATG": ["AAG", "ACG", "AGG", "ATA", "ATC", "ATT", "CTG", "GTG", "TTG"], "ATT": ["AAT", "ACT", "AGT", "ATA", "ATC", "ATG", "CTT", "GTT", "TTT"], "CAA": ["AAA", "CAC", "CAG", "CAT", "CCA", "CGA", "CTA", "GAA"], "CAC": ["AAC", "CAA", "CAG", "CAT", "CCC", "CGC", "CTC", "GAC", "TAC"], "CAG": ["AAG", "CAA", "CAC", "CAT", "CCG", "CGG", "CTG", "GAG"], "CAT": ["AAT", "CAA", "CAC", "CAG", "CCT", "CGT", "CTT", "GAT", "TAT"], "CCA": ["ACA", "CAA", "CCC", "CCG", "CCT", "CGA", "CTA", "GCA", "TCA"], "CCC": ["ACC", "CAC", "CCA", "CCG", "CCT", "CGC", "CTC", "GCC", "TCC"], "CCG": ["ACG", "CAG", "CCA", "CCC", "CCT", "CGG", "CTG", "GCG", "TCG"], "CCT": ["ACT", "CAT", "CCA", "CCC", "CCG", "CGT", "CTT", "GCT", "TCT"], "CGA": ["AGA", "CAA", "CCA", "CGC", "CGG", "CGT", "CTA", "GGA"], "CGC": ["AGC", "CAC", "CCC", "CGA", "CGG", "CGT", "CTC", "GGC", "TGC"], "CGG": ["AGG", "CAG", "CCG", "CGA", "CGC", "CGT", "CTG", "GGG", "TGG"], "CGT": ["AGT", "CAT", "CCT", "CGA", "CGC", "CGG", "CTT", "GGT", "TGT"], "CTA": ["ATA", "CAA", "CCA", "CGA", "CTC", "CTG", "CTT", "GTA", "TTA"], "CTC": ["ATC", "CAC", "CCC", "CGC", "CTA", "CTG", "CTT", "GTC", "TTC"], "CTG": ["ATG", "CAG", "CCG", "CGG", "CTA", "CTC", "CTT", "GTG", "TTG"], "CTT": ["ATT", "CAT", "CCT", "CGT", "CTA", "CTC", "CTG", "GTT", "TTT"], "GAA": ["AAA", "CAA", "GAC", "GAG", "GAT", "GCA", "GGA", "GTA"], "GAC": ["AAC", "CAC", "GAA", "GAG", "GAT", "GCC", "GGC", "GTC", "TAC"], "GAG": ["AAG", "CAG", "GAA", "GAC", "GAT", "GCG", "GGG", "GTG"], "GAT": ["AAT", "CAT", "GAA", "GAC", "GAG", "GCT", "GGT", "GTT", "TAT"], "GCA": ["ACA", "CCA", "GAA", "GCC", "GCG", "GCT", "GGA", "GTA", "TCA"], "GCC": ["ACC", "CCC", "GAC", "GCA", "GCG", "GCT", "GGC", "GTC", "TCC"], "GCG": ["ACG", "CCG", "GAG", "GCA", "GCC", "GCT", "GGG", "GTG", "TCG"], "GCT": ["ACT", "CCT", "GAT", "GCA", "GCC", "GCG", "GGT", "GTT", "TCT"], "GGA": ["AGA", "CGA", "GAA", "GCA", "GGC", "GGG", "GGT", "GTA"], "GGC": ["AGC", "CGC", "GAC", "GCC", "GGA", "GGG", "GGT", "GTC", "TGC"], "GGG": ["AGG", "CGG", "GAG", "GCG", "GGA", "GGC", "GGT", "GTG", "TGG"], "GGT": ["AGT", "CGT", "GAT", "GCT", "GGA", "GGC", "GGG", "GTT", "TGT"], "GTA": ["ATA", "CTA", "GAA", "GCA", "GGA", "GTC", "GTG", "GTT", "TTA"], "GTC": ["ATC", "CTC", "GAC", "GCC", "GGC", "GTA", "GTG", "GTT", "TTC"], "GTG": ["ATG", "CTG", "GAG", "GCG", "GGG", "GTA", "GTC", "GTT", "TTG"], "GTT": ["ATT", "CTT", "GAT", "GCT", "GGT", "GTA", "GTC", "GTG", "TTT"], "TAA": ["AAA", "CAA", "GAA", "TAC", "TAT", "TCA", "TTA"], "TAC": ["AAC", "CAC", "GAC", "TAT", "TCC", "TGC", "TTC"], "TAG": ["AAG", "CAG", "GAG", "TAC", "TAT", "TCG", "TGG", "TTG"], "TAT": ["AAT", "CAT", "GAT", "TAC", "TCT", "TGT", "TTT"], "TCA": ["ACA", "CCA", "GCA", "TCC", "TCG", "TCT", "TTA"], "TCC": ["ACC", "CCC", "GCC", "TAC", "TCA", "TCG", "TCT", "TGC", "TTC"], "TCG": ["ACG", "CCG", "GCG", "TCA", "TCC", "TCT", "TGG", "TTG"], "TCT": ["ACT", "CCT", "GCT", "TAT", "TCA", "TCC", "TCG", "TGT", "TTT"], "TGA": ["AGA", "CGA", "GGA", "TCA", "TGC", "TGG", "TGT", "TTA"], "TGC": ["AGC", "CGC", "GGC", "TAC", "TCC", "TGG", "TGT", "TTC"], "TGG": ["AGG", "CGG", "GGG", "TCG", "TGC", "TGT", "TTG"], "TGT": ["AGT", "CGT", "GGT", "TAT", "TCT", "TGC", "TGG", "TTT"], "TTA": ["ATA", "CTA", "GTA", "TCA", "TTC", "TTG", "TTT"], "TTC": ["ATC", "CTC", "GTC", "TAC", "TCC", "TGC", "TTA", "TTG", "TTT"], "TTG": ["ATG", "CTG", "GTG", "TCG", "TGG", "TTA", "TTC", "TTT"], "TTT": ["ATT", "CTT", "GTT", "TAT", "TCT", "TGT", "TTA", "TTC", "TTG"]}
    codonToSlowerValidMutations, codonToFasterValidMutations = getCodonToValidMutationDicts(codonToValidMutations, codonToSpeed)
    percentThatIsRamp = percentThatIsRampParam
    ribosomeWindowLength = ribosomeWindowLengthParam
    middleFunc = middleFuncParam
    args = argsParam

def constructDict(speedSeqs):
    """
    Creates a dictionary from a list of tuples. Used to convert the results of the findWindowMeansThreaded function into a dictionary.
    """
    windowMeans = {}
    for tuple in speedSeqs:
        windowMeans[tuple[0]] = tuple[1]
    return windowMeans

def writeCodonSpeedsFile(codonToSpeed):
    """
    Writes the calculated codon speeds (wij) to a csv file.
    """
    stopCodons = ["TAA", "TAG", "TGA"]
    with open(args.wij, "w", newline="") as csvfile:
        for codon, speed in codonToSpeed.items():
            if codon not in stopCodons:
                csvfile.write(f"{codon},{speed}\n")

def outputRampScores(rampScores, outPath):
    """
    Writes the ramp scores to the outPath file.
    """
    with open(outPath, "w") as outf:
        outf.write("header\tCDS length\tramp length\tnon-ramp min mean\tramp min mean\tramp score\tnum ramp region status changing mutations\tnum non-ramp region status changing mutations\tnum reasonably valid ramp mutations\tnum reasonably valid non-ramp mutations\tramp robustness\n")
        for rampScore in rampScores:
            if rampScore != "None":
                outf.write("\t".join(rampScore) + "\n")

def determineNewCutoff(pool, seqToSpeed):
    """
    Determines the cutoff percentile for ramp sequences using outliers in bottleneck frequencies for each gene percentile.
    See getCutoffValue function for more information.
    """
    seqHeaders = list(seqToSpeed.keys())
    poolArgs = [(getBottleneck, seqHeader) for seqHeader in seqHeaders]
    bottlenecks = pool.map(wrapperFunction, poolArgs)
    validatePoolResults(bottlenecks, pool)
    counts = {}
    for percents in bottlenecks:
        for p in percents:
            place = int(p[0])
            if not place in counts:
                counts[place] = 0
            counts[place] +=1
    for perc in range(1, 101):
        if not perc in counts:
            counts[perc] = 0
    percentThatIsRamp, message = getCutoffValue(counts)
    return percentThatIsRamp, message

def wrapperFunction(funcAndArgs):
    """
    Wrapper function for catching and returning exceptions in pool multiprocessing functions so they can be handled in the main process.
    """
    func, args = funcAndArgs
    try:
        return func(args)
    except Exception as e:

        excType, excValue, excTraceback = sys.exc_info()
        errorMessage = traceback.format_exception(excType, excValue, excTraceback)
        errorMessage = "".join(errorMessage)
        errorMessage = f"\nUNEXPECTED ERROR:\nExtRamp has crashed. Please send us the command you ran, the error below, and any other pertinent information.\n\n{errorMessage}"

        return Exception(errorMessage)

def validatePoolResults(results, pool=None):
    """
    If any of the workers return an exception, closes the pool and exits the program.
    """    
    if any(isinstance(result, Exception) for result in results):
        # get the first exception and print its message
        for result in results:
            if isinstance(result, Exception):
                sys.stderr.write(f"{result}\n")
                break
        if pool:
            pool.close()
            pool.join()
        sys.stderr.write(f"Exiting...\n")
        sys.exit()

def getRampStrength(speeds, percentThatIsRamp):
    """
    Gets the strength for the ramp for the speeds. A positive score indicates the presence of a ramp while 
    a negative score indicates an absence of a ramp. The score is calculated using the following formula
    zscore non-ramp region minimum window mean - zscore ramp region minimum window mean
    Note that the ramp region is up to the ramp cutoff, even if the ramp is shorter or longer than the ramp cutoff.
    Using the cutoff instead of the actual length makes comparisons between wild/mutant sequences easier.
    Returns the non-ramp min mean, ramp min mean, and ramp strength
    """
    windowMeans = getSlidingMean(args.middle, speeds, ribosomeWindowLength)

    # if the window means are all about same, zscores can't be calculated, and neither can a strength score
    # this typically triggered by long repeat sequences
    if np.std(windowMeans) < 1e-8:
        return "NA", "NA", "NA"
    
    zWindowMeans = zscore(windowMeans)

    rampCutoff = int(len(zWindowMeans) * percentThatIsRamp / 100) + 1 # plus 1 to include the last window of the ramp region

    rampMinMean = min(zWindowMeans[:rampCutoff])
    # get the non-ramp min mean
    if rampCutoff < len(zWindowMeans):
        nonRampMinMean = min(zWindowMeans[rampCutoff:])
        rampStrength = nonRampMinMean - rampMinMean
    elif rampCutoff > len(zWindowMeans):
        # if the ramp is most of the sequence except for a piece smaller than a window, get the mean of the codon speeds not in the ramp
        nonRampMinMean = middleFunc[args.middle](speeds[rampCutoff:])
        rampStrength = nonRampMinMean - rampMinMean
    else:
        # if the ramp is the entire sequence, then a strength score can't be calculated
        return "NA", "NA", "NA"
    
    return nonRampMinMean, rampMinMean, rampStrength

def getCodonToValidMutationDicts(codonToValidMutations, codonToSpeed):
    """
    Returns dictionaries of codons to all codons that are slower (codonToSlowerValidMutations) or
    faster (codonToFasterValidMutations) than the codon.
    """
    # get the possible faster and slower mutations for each codon
    codonToFasterValidMutationSpeeds = {}
    codonToSlowerValidMutationSpeeds = {}
    codonToFasterValidMutations = {}
    codonToSlowerValidMutations = {}
    for codon in codonToValidMutations:
        codonToFasterValidMutationSpeeds[codon] = {}
        codonToSlowerValidMutationSpeeds[codon] = {}
        for mutantCodon in codonToValidMutations[codon]:
            if mutantCodon in codonToSpeed:
                codonSpeed = codonToSpeed[codon] if codon in codonToSpeed else 0.0001
                if codonToSpeed[mutantCodon] > codonSpeed:
                    codonToFasterValidMutationSpeeds[codon][mutantCodon] = codonToSpeed[mutantCodon]
                elif codonToSpeed[mutantCodon] < codonSpeed:
                    codonToSlowerValidMutationSpeeds[codon][mutantCodon] = codonToSpeed[mutantCodon]
        # sort the faster and slower mutations by speed
        codonToFasterValidMutations[codon] = sorted(codonToFasterValidMutationSpeeds[codon], key=codonToFasterValidMutationSpeeds[codon].get, reverse=True)
        codonToSlowerValidMutations[codon] = sorted(codonToSlowerValidMutationSpeeds[codon], key=codonToSlowerValidMutationSpeeds[codon].get, reverse=False)
    return codonToSlowerValidMutations, codonToFasterValidMutations

def setUpGetNumMutations(speeds, windowMeans, percentThatIsRamp):
    """
    Returns the reciprocal speeds, the ramp region end index, the ramp minimum efficiency, and the non-ramp minimum efficiency.
    """
    reciprocalSpeeds = 1 / speeds
    rampRegionEnd = int(len(windowMeans) * percentThatIsRamp / 100)
    rampMin = np.min(windowMeans[:rampRegionEnd])
    nonRampMin = np.min(windowMeans[rampRegionEnd:])

    return reciprocalSpeeds, rampRegionEnd, rampMin, nonRampMin

def evaluateMutations(seqCodons, reciprocalSpeeds, rampWindowsToCheck, ribosomeWindowLength, comparisonValue, codonToReciprocalSpeed, codonToReasonablyValidMutations, mutationCriteria, comparisonType):
    """
    Checks each codon corresponding to each window index in rampWindowsToCheck to for single nucleotide mutations that alter 
    the ramp status. Mutations are only considered if the resulting mutant codon's efficiency is above (if comparisonType == "above") 
    or below (if comparisonType == "below") the wild-type codon efficiency. Valid mutant codons are determined by the 
    codonToReasonablyValidMutations dictionary. If mutationCriteria is "set", the "numMutations" that can alter the ramp status are
    obtained using a dictionary of positions to sets of ramp status altering mutations (reasonablyValidMutations). If mutationCriteria
    is "dict", "numMutations" that can alter the ramp status are obtained using a dictionary of positions to mutant codons to the number
    of windows they affect and only the mutations that affect all ramp windows (affectWindowsCriteria) are counted in "numMutations".
    """
    reasonablyValidMutations = {} # {seqPosition: set(mutantCodon)}
    
    affectWindowsCriteria = len(rampWindowsToCheck) # used to determine if a mutation must affect all ramp windows
    
    numMutations = 0  # Counts mutations meeting the criteria

    mutationPosToCodon = {} # seqPosition to set of mutantCodons or dict of mutantCodon to numMutations

    for windowIndex in rampWindowsToCheck:
        codons = seqCodons[windowIndex:windowIndex + ribosomeWindowLength]
        reciprocalSpeedsWindow = reciprocalSpeeds[windowIndex:windowIndex + ribosomeWindowLength]
        # pre-calculate the reciprocal sum for the initial array
        initialReciprocalSum  = np.sum(reciprocalSpeedsWindow)
        for codonIndex, codon in enumerate(codons):
            # the position of the codon in the dna sequence; prevents double counting codons in multiple windows
            seqPosition = windowIndex * 3 + codonIndex * 3
            if seqPosition not in reasonablyValidMutations:
                reasonablyValidMutations[seqPosition] = set()
            # get a dictionary of all the mutant codons that are higher/lower (based on comparisonType) than the current speed
            mutantCodonsToCheck = codonToReasonablyValidMutations[codon]
            reasonablyValidMutations[seqPosition].update(mutantCodonsToCheck)
            for mutantCodon in mutantCodonsToCheck:
                mutantCodonReciprocalSpeed = codonToReciprocalSpeed[mutantCodon]
                # get new windowMean with new speed inserted
                newReciprocalSum = initialReciprocalSum + mutantCodonReciprocalSpeed - reciprocalSpeedsWindow[codonIndex]
                newWindowMean = ribosomeWindowLength / newReciprocalSum
                if (comparisonType == "above" and newWindowMean > comparisonValue) or (comparisonType == "below" and newWindowMean < comparisonValue):
                    if mutationCriteria == "set":
                        if seqPosition not in mutationPosToCodon:
                            mutationPosToCodon[seqPosition] = set()
                        mutationPosToCodon[seqPosition].add(mutantCodon)
                    else:  # For dict, used to track mutations that need to affect a certain number of windows
                        if seqPosition not in mutationPosToCodon:
                            mutationPosToCodon[seqPosition] = {}
                        if mutantCodon not in mutationPosToCodon[seqPosition]:
                            mutationPosToCodon[seqPosition][mutantCodon] = 0
                        mutationPosToCodon[seqPosition][mutantCodon] += 1
                        # if the number of mutations at a particular position reaches the rampWindowsToAffect threshold, increment numMutations                        
                        if affectWindowsCriteria and mutationPosToCodon[seqPosition][mutantCodon] == affectWindowsCriteria:
                            numMutations += 1
                else:
                    break  # No need to check further if one condition fails, assuming sorted by relevance

    if mutationCriteria == "set":
        # Calculate the total number of unique mutations across all positions for set criteria
        for mutations in mutationPosToCodon.values():
            numMutations += len(mutations)

    return numMutations, reasonablyValidMutations

def calculateRampRobustness(rampRobustnessMultiplier, numRampRegionChangingMutations, numNonRampRegionChangingMutations, reasonablyValidRampMutations, reasonablyValidNonRampMutations):
    """
    Ramp robustness is the fraction of mutations that are reasonably valid that do not change the ramp status of the sequence, multiplied 
    by -1 if the sequence does not have a ramp. A single nucleotide mutation is reasonably valid if the mutant codon's corresponding efficiency
    is above () or below () the wild-type codon efficiency. For example, if a sequence without a ramp has a ramp region with 10 mutations 
    out of 100 that are reasonably valid and the non-ramp region has 20 mutations out of 100 that are reasonably valid, the ramp robustness 
    would be -1 * (1 - 30/100) = 0.7. 0.7 means that 70% of the mutations that are reasonably valid do not create a ramp in the sequence.
    Returns the number of reasonably valid ramp mutations, the number of reasonably valid non-ramp mutations, and the ramp robustness score.
    """
    numChangingMutations = numRampRegionChangingMutations + numNonRampRegionChangingMutations
    numReasonablyValidRampMutations = sum(len(mutations) for mutations in reasonablyValidRampMutations.values())
    numReasonablyValidNonRampMutations = sum(len(mutations) for mutations in reasonablyValidNonRampMutations.values())
    totalReasonablyValidMutations = numReasonablyValidRampMutations + numReasonablyValidNonRampMutations
    rampRobustness = rampRobustnessMultiplier * (1 - (numChangingMutations / totalReasonablyValidMutations))
    
    return numReasonablyValidRampMutations, numReasonablyValidNonRampMutations, rampRobustness

def getNumDestroyingMutations(seqCodons, speeds, windowMeans, codonToReciprocalSpeed, codonToSlowerValidMutations, codonToFasterValidMutations, percentThatIsRamp, ribosomeWindowLength):
    """
    Gets the number of reasonably valid mutations that destroy the ramp of the sequence. Ramp region mutations that destroy 
    the ramp must increase ALL ramp windows above the non-ramp minimum efficiency. Non-ramp region mutations that destroy 
    the ramp must decrease ANY non-ramp window below the ramp minimum efficiency.
    Returns the number of ramp region ramp destroying mutations, the number of non-ramp region ramp destroying mutations, the number of
    reasonably valid ramp mutations, the number of reasonably valid non-ramp mutations, and the ramp robustness score.
    """
    reciprocalSpeeds, rampRegionEnd, rampMin, nonRampMin = setUpGetNumMutations(speeds, windowMeans, percentThatIsRamp)

    # check if all of the ramp windows can be mutated above the nonRampMin by one mutation
    rampWindowsToCheck = np.where(windowMeans[:rampRegionEnd] < nonRampMin)[0].tolist()
    # check if the windows are close enough for a single mutation to affect all of them
    windowsAllOverlap = True
    for index in range(len(rampWindowsToCheck)-1):
        if rampWindowsToCheck[index] + ribosomeWindowLength < rampWindowsToCheck[index+1]:
            windowsAllOverlap = False
            break
    if windowsAllOverlap:
        numRampRegionDestroyingMutations, reasonablyValidRampMutations = evaluateMutations(seqCodons, reciprocalSpeeds, rampWindowsToCheck, ribosomeWindowLength, nonRampMin, codonToReciprocalSpeed, codonToFasterValidMutations, "dict", "above")
    else:
        numRampRegionDestroyingMutations = 0
        reasonablyValidRampMutations = {}
    # check if any of the non-ramp windows can be mutated below the rampMin
    numNonRampRegionDestroyingMutations, reasonablyValidNonRampMutations = evaluateMutations(seqCodons, reciprocalSpeeds, range(rampRegionEnd,len(windowMeans)), ribosomeWindowLength, rampMin, codonToReciprocalSpeed, codonToSlowerValidMutations, "set", "below")

    numReasonablyValidRampMutations, numReasonablyValidNonRampMutations, rampRobustness = calculateRampRobustness(1, numRampRegionDestroyingMutations, numNonRampRegionDestroyingMutations, reasonablyValidRampMutations, reasonablyValidNonRampMutations)
    
    return numRampRegionDestroyingMutations, numNonRampRegionDestroyingMutations, numReasonablyValidRampMutations, numReasonablyValidNonRampMutations, rampRobustness

def getNumCreatingMutations(seqCodons, speeds, windowMeans, codonToReciprocalSpeed, codonToSlowerValidMutations, codonToFasterValidMutations, percentThatIsRamp, ribosomeWindowLength):
    """
    Gets the number of reasonably valid mutations that create a ramp in the sequence. Ramp region mutations that create 
    a ramp must decrease ANY ramp windows below the non-ramp minimum efficiency. Non-ramp region mutations that
    create a ramp must increase ALL non-ramp windows below the ramp minimum efficiency.
    Returns the number of ramp region ramp creating mutations, the number of non-ramp region ramp creating mutations, the number of
    reasonably valid ramp mutations, the number of reasonably valid non-ramp mutations, and the ramp robustness score.
    """
    reciprocalSpeeds, rampRegionEnd, rampMin, nonRampMin = setUpGetNumMutations(speeds, windowMeans, percentThatIsRamp)

    # check any of the ramp window means can be mutated below the nonRampMin
    numRampRegionCreatingMutations, reasonablyValidRampMutations = evaluateMutations(seqCodons, reciprocalSpeeds, range(0,rampRegionEnd), ribosomeWindowLength, nonRampMin, codonToReciprocalSpeed, codonToSlowerValidMutations, "set", "below")
    # check if all of the non-ramp windows can be mutated above the rampMin by one mutation
    nonRampWindowsToCheck = np.where(windowMeans[rampRegionEnd:] < rampMin)[0].tolist()
    # check if the windows are close enough for a single mutation to affect all of them
    windowsAllOverlap = True
    for index in range(len(nonRampWindowsToCheck)-1):
        if nonRampWindowsToCheck[index] + ribosomeWindowLength < nonRampWindowsToCheck[index+1]:
            windowsAllOverlap = False
            break        
    if windowsAllOverlap:
        numNonRampRegionCreatingMutations, reasonablyValidNonRampMutations = evaluateMutations(seqCodons, reciprocalSpeeds, nonRampWindowsToCheck, ribosomeWindowLength, rampMin, codonToReciprocalSpeed, codonToFasterValidMutations, "dict", "above")
    else:
        numNonRampRegionCreatingMutations = 0
        reasonablyValidNonRampMutations = {}

    numReasonablyValidRampMutations, numReasonablyValidNonRampMutations, rampRobustness = calculateRampRobustness(-1, numRampRegionCreatingMutations, numNonRampRegionCreatingMutations, reasonablyValidRampMutations, reasonablyValidNonRampMutations)
    
    return numRampRegionCreatingMutations, numNonRampRegionCreatingMutations, numReasonablyValidRampMutations, numReasonablyValidNonRampMutations, rampRobustness

def getRampRobustness(seq, speeds, codonToSpeed, codonToSlowerValidMutations, codonToFasterValidMutations, percentThatIsRamp, ribosomeWindowLength):
    """
    A wrapper for the getNumDestroyingMutations and getNumCreatingMutations functions. Determines if the sequence has a ramp, 
    then calls the appropriate function to get the number of mutations that change the ramp status of the sequence.
    Returns the number of ramp region ramp status changing mutations, the number of non-ramp region ramp status changing mutations, the number of
    reasonably valid ramp mutations, the number of reasonably valid non-ramp mutations, and the ramp robustness score.
    """
    seqCodons = np.array([seq[i:i+3] for i in range(0, len(seq), 3)])
    windowMeans = hmeanSliding(speeds, ribosomeWindowLength)
    minIndex = np.argmin(windowMeans)
    minWindowMeanPercent = 100 * (minIndex / len(windowMeans))
    # get codon to reciprocal speed
    codonToReciprocalSpeed = {codon: 1 / speed for codon, speed in codonToSpeed.items()}
    if minWindowMeanPercent <= percentThatIsRamp: # if the sequence has a ramp
        numRampRegionStatusChangingMutations, numNonRampRegionStatusChangingMutations, numReasonablyValidRampMutations, numReasonablyValidNonRampMutations, rampRobustness = getNumDestroyingMutations(seqCodons, speeds, windowMeans, codonToReciprocalSpeed, codonToSlowerValidMutations, codonToFasterValidMutations, percentThatIsRamp, ribosomeWindowLength)
    else:
        numRampRegionStatusChangingMutations, numNonRampRegionStatusChangingMutations, numReasonablyValidRampMutations, numReasonablyValidNonRampMutations, rampRobustness = getNumCreatingMutations(seqCodons, speeds, windowMeans, codonToReciprocalSpeed, codonToSlowerValidMutations, codonToFasterValidMutations, percentThatIsRamp, ribosomeWindowLength)

    return numRampRegionStatusChangingMutations, numNonRampRegionStatusChangingMutations, numReasonablyValidRampMutations, numReasonablyValidNonRampMutations, rampRobustness

def getRampLength(speeds):
    windowMeans = getSlidingMean(args.middle, speeds, ribosomeWindowLength)
    mean = hmean(speeds)
    bestMean = min(windowMeans)
    pos = np.where(windowMeans == bestMean)[0][0] # gets the first minimum
    perc = float(100*(pos/float(len(windowMeans)))) # percentages x 100
    if perc <= percentThatIsRamp: # if the sequence has a ramp
        lastRampWindow = getLastRampWindowIndex(pos, windowMeans, mean)
        rampLen = (lastRampWindow+ribosomeWindowLength)*3 + 3 # add the start codon
        return rampLen 
    else: # if the sequence doesn't have a ramp, return -1
        return -1

def scoreRamps(seqTuple):
    """
    Gets the strength and robustness for an input sequence tuple, calling getRampStrength and getRampRobustness.
    Returns the header, length of the input CDS (including start and stop codons), length of calculated ramp (including start codon), 
    minimum of the non-ramp region, minimum of the ramp region, ramp strength, number of ramp region status changing mutations, 
    number of non-ramp region status changing mutations, number of reasonably valid ramp mutations, number of reasonably 
    valid non-ramp mutations, and the ramp robustness score.
    """
    header = seqTuple[0]
    seq = seqTuple[1]
    speeds = seqToSpeed[header]

    nonRampMinMean, rampMinMean, rampStrength = getRampStrength(speeds, percentThatIsRamp)
    numRampRegionStatusChangingMutations, numNonRampRegionStatusChangingMutations, numReasonablyValidRampMutations, numReasonablyValidNonRampMutations, rampRobustness = getRampRobustness(seq, speeds, codonToSpeed, codonToSlowerValidMutations, codonToFasterValidMutations, percentThatIsRamp, ribosomeWindowLength)

    cdsLength = len(speeds) * 3 + 6 # x3 to get nucleotides and add start and stop back 
    rampLength = getRampLength(speeds)
    
    return tuple([header, str(cdsLength), str(rampLength), str(nonRampMinMean), str(rampMinMean), str(rampStrength), 
                str(numRampRegionStatusChangingMutations), str(numNonRampRegionStatusChangingMutations), 
                str(numReasonablyValidRampMutations), str(numReasonablyValidNonRampMutations), str(rampRobustness)])

if __name__ == "__main__":
    freeze_support()
    
    # parse arguments
    args = makeArgParser()
    if args.verbose:
        sys.stderr.write(f"Running ExtRamp version {VERSION}\n\n")
        initstart = time.time()
    ribosomeWindowLength = args.window
    middleFunc = {"hmean": hmean, "gmean": gmean, "mean": mean, "median": median}
    if args.middle not in middleFunc:
        sys.stderr.write("Error: args.middle must be one of the following options:\n\thmean\n\tgmean\n\tmean\n\tmedian\nExiting...\n")
        sys.exit()

    # read sequences
    if args.verbose:
        sys.stderr.write("Reading Sequences... ")
        start = time.time()
    seqTupleArray, message = readSeqFile(args, args.input)
    if args.verbose:
        sys.stderr.write(f"took: {time.time() - start}\n")
        sys.stderr.write(f"{message}\n")
        sys.stderr.write(f"Total valid sequences: {len(seqTupleArray)}\n\n")
    if len(seqTupleArray) == 0:
        sys.stderr.write("No sequences passed the initial filter. Ramp sequences were unable to be calculated.\n")
        sys.exit()

    # calculate codon speeds
    codonToSpeed = {}
    if args.verbose:
        sys.stderr.write("Calculating codon speeds... ")
        start = time.time()
    if args.tAI:
        codonToSpeed = csvToDict(args.tAI)
    else:
        if args.rscu:
            seqTupleArray_rscu, message = readSeqFile(args, args.rscu)
            if args.verbose:
                rscuTime = time.time() - start
            codonToSpeed = calcCodonSpeeds(seqTupleArray_rscu)
        else:
            codonToSpeed = calcCodonSpeeds(seqTupleArray)
    if args.verbose:
        sys.stderr.write(f"took: {time.time() - start}\n")
        if args.rscu and message != "":
            message = message.replace("\n","\n\t")
            sys.stderr.write(f"\tReading RSCU file took: {rscuTime}:{message}\n")
    if args.wij:
        if args.verbose:
            sys.stderr.write("Writing codon speeds to file... ")
            start = time.time()
        writeCodonSpeedsFile(codonToSpeed)
        if args.verbose:
            sys.stderr.write(f"took: {time.time() - start}\n")

    # calculate sequence speeds
    if args.verbose:
        sys.stderr.write("Calculating sequence speeds... ")
        start = time.time()
    pool = Pool(processes = args.threads, initializer = init_pool, initargs = (codonToSpeed,))
    poolArgs = [(calcSeqSpeed, record) for record in seqTupleArray]
    seqToSpeed = pool.map(wrapperFunction, poolArgs)
    validatePoolResults(seqToSpeed, pool)
    seqToSpeed = createSpeedsDict(seqToSpeed) # header -> tuple of array of codon efficiency values
    pool.close()
    pool.join()
    if args.verbose:
        sys.stderr.write(f"took: {time.time() - start}\n")
    if args.speeds:
        if args.verbose:
            sys.stderr.write("Writing speeds data to file... ")
            start = time.time()
        writeEfficiencyFile(args.speeds, seqToSpeed)
        if args.verbose:
            sys.stderr.write(f"took: {time.time() - start}\n")

    # create ramp pool
    if args.verbose:
        sys.stderr.write("Creating ramp pool... ")
        start = time.time()
    pool = Pool(processes = args.threads, initializer = init_pool2, initargs = (seqToSpeed, ribosomeWindowLength, middleFunc, args))
    if args.verbose:
        sys.stderr.write(f"took: {time.time() - start}\n")

    # determine ramp cutoff
    percentThatIsRamp = args.cutoff
    if args.determine_cutoff:
        if args.verbose:
            sys.stderr.write("Determining ramp cutoff percentage... ")
            start = time.time()
        percentThatIsRamp, message = determineNewCutoff(pool, seqToSpeed)
        if args.verbose:
            sys.stderr.write(f"took: {time.time() - start}\n")
            sys.stderr.write(f"{message}\n")
            sys.stderr.write(f"\nThe cutoff percentage is {percentThatIsRamp}%\n")

    # isolate ramp sequences
    if args.verbose:
        sys.stderr.write("\nIsolating ramp sequences... ")
        start = time.time()
    if args.stdev >= 0:
        rampFunction = isolateRamp
    else:
        rampFunction = isolateRampHmean
    poolArgs = [(rampFunction, (seqTuple, percentThatIsRamp)) for seqTuple in seqTupleArray]
    rampSeqs = pool.map(wrapperFunction, poolArgs)
    validatePoolResults(rampSeqs, pool)
    if args.verbose:
        sys.stderr.write(f"took: {time.time() - start}\n")
        start = time.time()

    # write ramp sequences to a FASTA file
    outputRampSeqs(rampSeqs, args)
    if args.verbose:
        sys.stderr.write(f"\nOutputting ramps took: {time.time() - start}\n")

    # output window mean speed values in a FASTA-like file if indicated
    if args.vals:
        if args.verbose:
            sys.stderr.write("Calculating window mean speeds for file... ")
            start = time.time()
        poolArgs = [(findWindowMeansThreaded, seq) for seq in seqToSpeed]
        windowMeans = pool.map(wrapperFunction, poolArgs)
        validatePoolResults(windowMeans, pool)
        windowMeans = constructDict(windowMeans)
        if args.verbose:
            sys.stderr.write(f"took: {time.time() - start}\n")
            sys.stderr.write("Writing window mean speeds file... ")
            start = time.time()
        writeEfficiencyFile(args.vals, windowMeans)
        if args.verbose:
            sys.stderr.write(f"took: {time.time() - start}\n")

    # output ramp scores to a file
    if args.score:
        pool.close()
        pool.join()

        # create score pool
        if args.verbose:
            sys.stderr.write("Creating score pool... ")
            start = time.time()
        pool = Pool(processes = args.threads, initializer = init_pool3, initargs = (codonToSpeed, seqToSpeed, percentThatIsRamp, ribosomeWindowLength, middleFunc, args))
        if args.verbose:
            sys.stderr.write(f"took: {time.time() - start}\n")

        if args.verbose:
            sys.stderr.write("Calculating ramp scores... ")
            start = time.time()
        poolArgs = [(scoreRamps, seqTuple) for seqTuple in seqTupleArray]
        rampScores = pool.map(wrapperFunction, poolArgs)
        validatePoolResults(rampScores, pool)
        if args.verbose:
            sys.stderr.write(f"took: {time.time() - start}\n")
            sys.stderr.write("Writing ramp scores file... ")
            start = time.time()
        outputRampScores(rampScores, args.score)
        if args.verbose:
            sys.stderr.write(f"took: {time.time() - start}\n")

    # close the last pool used and finish the program
    pool.close()
    pool.join()
    if args.verbose:
        sys.stderr.write(f"\nTotal time: {time.time() - initstart}\n")
