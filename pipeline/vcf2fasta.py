# This backup was taken before optimizing reading the gff to only read chromosomes found in the VCF

import argparse
import sys
import os
import re
import time
import pysam
import collections

# helper.py
import helper

MANE_SELECT_FEATURES = ["CDS", "RNase_MRP_RNA", "antisense_RNA", "exon", "lnc_RNA", "mRNA", "snRNA", "snoRNA", "telomerase_RNA"]


def get_parser():
    # parse arguments
    parser = argparse.ArgumentParser(
        prog="vcf2fasta.py",
        formatter_class=argparse.RawTextHelpFormatter,
        description="""
        Converts regions/intervals in the genome into FASTA alignments
        provided a VCF file, a GFF/GTF file, and FASTA reference.\n""",
        epilog="""
        All files must be indexed. So before running the code make sure
        that your reference FASTA file is indexed:

        samtools faidx genome.fas

        BGZIP compress and TABIX index your VCF file:

        bgzip variants.vcf
        tabix variants.vcf.gz

        The GFF/GTF file does not need to be indexed.

        examples:
        python vcf2fasta.py -f genome.fas -v variants.vcf.gz -g intervals.gff -e CDS
        \n""",
    )
    parser.add_argument("--fasta", "-f", metavar="GENOME", type=str, required=True, help="FASTA file with the reference genome.")
    parser.add_argument("--vcf", "-v", metavar="VCF", type=str, required=True, help="a tabix-indexed VCF file.")
    parser.add_argument("--gff", "-g", metavar="GFF/GTF", type=str, required=False, help="GFF/GTF file.")
    parser.add_argument("--bed", metavar="BED", type=str, required=False, help="BED file.")
    parser.add_argument("--feat", "-e", default="gene", metavar="FEAT", type=str, required=False, help="feature/annotation in the GFF file (i.e. gene, CDS, intron). Use multiple features by separating them with commas; default: 'gene'.\n")
    parser.add_argument("--blend", "-b", action="store_true", default=False, help="concatenate GFF entries of FEAT into a single alignment. Useful for CDS. (default: False)")
    parser.add_argument("--inframe", "-i", action="store_true", default=False, help="force the first codon of the sequence to be inframe. Useful for incomplete CDS. (default: False)")
    parser.add_argument("--out", "-o", metavar="OUT", type=str, default="vcf2fasta", help="provide a name for the output directory (optional)")
    parser.add_argument("--addref", "-r", action="store_true", default=False, help="include the reference sequence in the FASTA alignment (default: False)")
    parser.add_argument("--skip", "-s", action="store_true", default=False, help="skips features without variants (default: False)")
    parser.add_argument("--force", action="store_true", default=False, help="no prompt for file overwritting (default: False)")
    parser.add_argument("--skipscaffoldcontigs", action="store_true", default=False, help="skips reading in scaffolds and contigs (NT and NW sequences) (default: False)")
    parser.add_argument("--mane", action="store_true", default=False, help="only writes out results for MANE select transcripts (default: False)")
    parser.add_argument("--byvariant", action="store_true", default=False, help="writes out separate FASTA files for each variant/position. If selected, 'addref' is automatically set to True (default: False)")
    
    return parser


# this function gets everything processed to produce a set
# of sequences in a dictionary, along with other info
def getSequences(intervals, gene, feat, ref, vcf, ploidy, phased, samples, ref_chrom_format, vcf_chrom_format, args):
    seqs = collections.defaultdict()
    variants_inserted = set()
    feat_str = "_" if feat == "" else "_" + feat + "_"
    # prep sequence dictionary # TODO can we remove prep?
    if args.blend:
        seqs[gene] = collections.defaultdict()
    else:
        feat_ind = 0
        for rec in intervals[gene]:
            featname = gene + feat_str + str(feat_ind)
            seqs[featname] = collections.defaultdict()
            feat_ind += 1

    if phased:
        for key in seqs.keys():
            for sample in samples:
                for i in range(ploidy):
                    seqs[key][sample + "_" + str(i)] = ""
            # add a dummy reference sequence if requested
            if args.addref:
                seqs[key]["REF_0"] = "" # TODO consider removing _0 since there is only one REF sequence
    else:
        for key in seqs.keys():
            for sample in samples:
                seqs[key][sample] = ""
            # add a dummy reference sequence if requested
            if args.addref:
                seqs[key]["REF"] = ""
    # initial values before looping through VCF slice
    varsites = 0
    feat_ind = 0
    featname = gene
    codon_start = collections.defaultdict(list)
    full_ref_seq = ""
    for rec in intervals[gene][feat]:
        if not args.blend:
            featname = gene + feat_str + str(feat_ind)
            feat_ind += 1
        # get a new copy of seqs for every feature
        tmpseqs = collections.defaultdict()
        # extract relevant info from GFF/BED
        if args.bed:
            chrom, start, end, strand, cs = rec[0], int(rec[1]), int(rec[2]), "+", 0
        else:
            chrom, start, end, strand, cs = (
                rec[0],
                int(rec[3]) - 1,
                int(rec[4]),
                rec[6],
                rec[7],
            )
        # add cs to codon_start (cs is the frame: 0,1,2 or ".")
        codon_start[featname].append(cs)
        # extract sequence from reference
        refseq = ref.fetch(chrom, start, end).upper()
        # propagate reference sequence to all samples
        for sample in seqs[featname].keys():
            tmpseqs[sample] = refseq
        # initialize positive or negative positions to extend
        # in case of indels
        posadd = 0
        # convert chrom names if needed
        vcf_chrom = chrom
        if ref_chrom_format != vcf_chrom_format:
            vcf_chrom = helper.convert_chrom_format(chrom, ref_chrom_format, vcf_chrom_format)
        vrecs = vcf.fetch(vcf_chrom, start, end) if vcf_chrom in vcf.header.contigs else []
        vrecs = [vrec for vrec in vrecs]
        for vrec in vrecs:
            # count variants within feature
            varsites += 1
            variants_inserted.add((vcf_chrom, vrec.pos))
            # get seq position of variant
            pos = vrec.pos - start - 1 + posadd
            if args.byvariant:
                # add sample corresponding to the variant
                seqs[featname][(vcf_chrom, vrec.pos)] = seqs[featname]["REF"]
                tmpseqs[(vcf_chrom, vrec.pos)] = tmpseqs["REF"]
            # get a dict of extended alleles, including indels
            # TODO make sure only the byvariant sample just added get the alt allele while the rest get the ref
            alleles, max_len = getAlleles(vrec, ploidy, phased, seqs[featname].keys(), (vcf_chrom, vrec.pos), args.addref, args.byvariant)
            ref_len = len(vrec.ref)
            for sample in alleles.keys():
                 # TODO for byvariant, how do I update the variants
                tmpseqs[sample] = UpdateSeq(alleles, sample, pos, ref_len, tmpseqs[sample])
            # this is the cumulative number of positions to add if
            # positions are taken or added to the sequence by indels
            posadd += max_len - ref_len
        for sample in seqs[featname].keys():
            if strand == "-":
                # reverse complement sequence if needed
                seqs[featname][sample] = seqs[featname][sample] + revcomp(tmpseqs[sample])
            else:
                seqs[featname][sample] = seqs[featname][sample] + tmpseqs[sample]    
    # adjust if inframe is on
    if args.inframe:
        if args.blend and codon_start[featname][0] != ".":
            if strand == "+" and codon_start[featname][0] != "0":
                for key in seqs[featname].keys():
                    seqs[featname][key] = seqs[featname][key][
                        int(codon_start[featname][0]) :
                    ]
            elif strand == "-" and codon_start[featname][-1] != "0":
                for key in seqs[featname].keys():
                    seqs[featname][key] = seqs[featname][key][
                        int(codon_start[featname][-1]) :
                    ]
        elif not args.blend:
            for featname in seqs.keys():
                for key in seqs[featname].keys():
                    seqs[featname][key] = seqs[featname][key][
                        int(codon_start[featname][0]) :
                    ]
    
    # remove "-" from sequences
    for featname in seqs.keys():
        for key in seqs[featname].keys():
            seqs[featname][key] = seqs[featname][key].replace("-", "")

    return seqs, varsites, variants_inserted


# main algorithm of vcf2fasta
# extracts alleles for each sample and expands them if indels
# deals with both phased and unphased data
def getAlleles(rec, ploidy, phased, samples, chrom_pos, addref, byvariant):
    # extract all alleles for a given SNP/var. pos.
    alleles = {i[0]: i[1].alleles for i in rec.samples.items()}
    # if byvariant, set all alt samples to the first alt allele # TODO what about multiple alts?
    if byvariant:
        for sample in samples:
            # if the sample is a tuple
            if isinstance(sample, tuple):
                if sample == chrom_pos:
                    alleles[sample] = (rec.alts[0],)
                else: # all other samples get the ref allele
                    alleles[sample] = (rec.ref,)
    else:
        # if no samples, create an "ALT" sample with the first alt allele # TODO what about multiple alts?
        if len(alleles) == 0 or byvariant:
            alleles = {"ALT": (rec.alts[0],)}
    # collapse list into alleles that are segregating and are not missing
    segregating = list(set(sum([[x for x in alleles[i]] for i in alleles.keys()], [])))
    # add ref allele if addref is true
    if addref:
        segregating = list(set(segregating + [rec.ref]))
    # get the length of the longest allele
    max_len = max([len(i) for i in segregating if i is not None])
    # make a dictionary of expanded alleles
    dict_expanded = {
        i: (i + "-" * (max_len - len(i))) for i in segregating if i is not None
    }
    # replace short alleles with expanded alleles for samples without missing data
    alleles_expanded = {
        i: [dict_expanded[j] for j in alleles[i]]
        for i in alleles.keys()
        if alleles[i][0] is not None
    }
    # add ref allele if addref is true
    if addref:
        alleles_expanded["REF"] = [dict_expanded[rec.ref]]
    # add one for missing data if any, and incorporate to alleles_expanded dict
    if None in segregating:
        dict_expanded[""] = "?" * max_len
        alleles_missing = {
            i: [dict_expanded[""] for j in range(ploidy)]
            for i in alleles.keys()
            if alleles[i][0] is None
        }
        for i in alleles_missing.keys():
            alleles_expanded[i] = alleles_missing[i]
    if phased:
        alleles_expanded = makePhased(alleles_expanded)  # make them phased
    else:
        alleles_expanded = {
            key: getIUPAC(alleles_expanded[key]) for key in alleles_expanded.keys()
        }
    return alleles_expanded, max_len


# takes a genotype dict from getAlleles and turns it into
# a phased dict with each haplotype as the sample + a zero-based index
def makePhased(alleles):
    alleles_phased = {}
    for samp in alleles.keys():
        for i in range(len(alleles[samp])):
            # append zero-based index
            alleles_phased[samp + "_" + str(i)] = alleles[samp][i]
    return alleles_phased


# function to update sequence dict based on alleles or collapsed alleles
def UpdateSeq(alleles, samp, pos, ref_len, seq):
    return seq[:pos] + alleles[samp] + seq[pos + ref_len :]


# same as UpdateSeqPhased but returns a collapsed genotype using IUPAC codes
# def UpdateSeqIUPAC(alleles, samp, pos, ref_len, seq):
#     return seq[:pos] + getIUPAC(alleles[samp]) + seq[pos+ref_len:]


# collapses list into string of IUPAC codes
### ISSUE
# produces shorter alignments sometimes
# needs fixing
###
def getIUPAC(x):
    """
    Collapses two or more alleles into a single
    IUPAC string
    """
    # first check if data is missing
    if len(x) == 1 or x[0][0] == "?":
        return x[0]
    elif len(list(set(x))) == 1:
        return x[0]
    else:
        iupacd = ""
        # loop through positions
        for i in range(len(x[0])):
            # keep nucleotides only
            nuc = list(set([y[i] for y in x]))
            # if single nuc
            if len(nuc) == 1:
                iupacd += nuc[0]
            else:
                nuc = [j for j in nuc if j != "-"]
                if len(nuc) == 1:
                    iupacd += nuc[0]
                # if multiple nucleotides per position
                if len(nuc) == 2:
                    if "A" in nuc and "G" in nuc:
                        iupacd += "R"
                    elif "A" in nuc and "T" in nuc:
                        iupacd += "W"
                    elif "A" in nuc and "C" in nuc:
                        iupacd += "M"
                    elif "C" in nuc and "T" in nuc:
                        iupacd += "Y"
                    elif "C" in nuc and "G" in nuc:
                        iupacd += "S"
                    elif "G" in nuc and "T" in nuc:
                        iupacd += "K"
                    else:
                        iupacd += "?"
                elif len(nuc) == 3:
                    if "A" in nuc and "T" in nuc and "C":
                        iupacd += "H"
                    elif "A" in nuc and "T" in nuc and "G":
                        iupacd += "D"
                    elif "G" in nuc and "T" in nuc and "C":
                        iupacd += "B"
                    elif "A" in nuc and "G" in nuc and "C":
                        iupacd += "V"
                    else:
                        iupacd += "?"
                elif len(nuc) == 4:
                    if "A" in nuc and "G" in nuc and "C" and "T" in nuc:
                        iupacd += "N"
                    else:
                        iupacd += "?"
        return iupacd


def printFasta(seqs, header, out, variant=None):
    """
    writes FASTA to file
    """
    # If the REF sequence is present, print it first
    if "REF" in seqs.keys():
        out.write(f">REF {header}\n{seqs['REF']}\n")
    elif "REF_0" in seqs.keys():
        out.write(f">REF_0 {header}\n{seqs['REF_0']}\n")
    # print the rest of the sequences
    for sample in seqs.keys():
        if variant is not None:
            if sample == variant:
                # print only the sample corresponding to the variant (used by args.byvariant)
                out.write(f">ALT {header}\n{seqs[sample]}\n")
        else:
            if sample != "REF" and sample != "REF_0":
                out.write(f">{sample} {header}\n{seqs[sample]}\n")


def processGeneNameGFF(lastfield):
    """
    Makes a list of all the annotation fields in the last column [8]
    delimited by ";"
    Input is the string of the last field in GFF
    """
    last = collections.defaultdict()
    for i in lastfield.split(";"):
        if "=" in i:
            x = i.split("=")
            last[re.sub('"| ', "", x[0])] = re.sub('"| ', "", x[1])
    return last


def processGeneNameGTF(lastfield):
    """
    Makes a list of all the annotation fields in the last column [8]
    delimited by ";"
    Input is the string of the last field in GFF
    """
    last = collections.defaultdict()
    for i in lastfield.split(";"):
        if " " in i:
            x = i.split(" ")
            last[x[0]] = re.sub('"| ', "", x[1])
    return last


def ReadBED(file):
    """
    returns a dictionary of with 0-based genomic intervals
    """
    with open(file) as f:
        lines = f.read().splitlines()
        bed = {"g" + str(i + 1): [lines[i].split("\t")] for i in range(len(lines))}
    return bed


def count_num_lines(filepath, chunk_size=65536):
    with open(filepath, 'rb') as inf:
        num_lines = 0
        while True:
            chunk = inf.read(chunk_size)
            if not chunk:
                break
            num_lines += chunk.count(b'\n')
        # Handle case where file doesn't end with a newline
        if chunk and not chunk.endswith(b'\n'):
            num_lines += 1
    return num_lines


def update_progress(line_counter, num_lines, t1):
    line_counter += 1
    if line_counter % 10000 == 0 or line_counter == num_lines:
        progress = make_progress_bar(line_counter, num_lines, t1, 70)
        print("\r", progress[0] % progress[1:], end="", flush=True)
    return line_counter


def ReadGFF(parser, vcf_chroms, vcf_chrom_format, args):
    """
    returns a nested dictionary named after every feature name
    as well as every feature name [3]
    """
    intervals_chrom_format = ""
    file = args.gff
    # find out if GTF or GFF
    if file.split(".")[-1].lower() == "gff" or file.split(".")[-1].lower() == "gff3":
        format = "gff"
    elif file.split(".")[-1].lower() == "gtf":
        format = "gtf"
    else:
        print("Cannot figure out GFF/GTF format. File should end with .gff or .gtf")
        sys.exit(parser.print_help())
    # start counting time
    t1 = time.time()
    # get the number of lines in the file using wc -l
    num_lines = count_num_lines(file)
    gff = {}
    line_counter = 0
    chrom = ""
    with open(file, "r") as f:
        if format == "gff":
            for line in f:
                line_counter = update_progress(line_counter, num_lines, t1)
                if args.skip and chrom and chrom not in vcf_chroms: # TODO fix this
                    # skip chromosomes not in VCF
                    while line.startswith(chrom):
                        line = f.readline()
                        line_counter = update_progress(line_counter, num_lines, t1)
                if line[0] != "#":
                    fields = line.rstrip().split("\t")
                    last = processGeneNameGFF(fields[8])
                    feature = fields[2]
                    if feature not in args.feat:
                        continue
                    chrom = fields[0]
                    if args.skipscaffoldcontigs:
                        # skip scaffold and contig sequences
                        if chrom.startswith("NT_") or chrom.startswith("NW_"):
                            continue
                    if intervals_chrom_format == "":
                        intervals_chrom_format = helper.get_chrom_format(chrom)
                        # convert all vcf_chroms to intervals format if necessary
                        if args.skip and intervals_chrom_format != vcf_chrom_format:
                            vcf_chroms = {helper.convert_chrom_format(c, vcf_chrom_format, intervals_chrom_format) for c in vcf_chroms}
                    if args.skip and chrom not in vcf_chroms:
                        continue
                    if args.mane and feature in MANE_SELECT_FEATURES: # TODO
                        # TODO if no lines have MANESelect tag, ignore this filter
                        if "tag" not in last or "MANESelect" not in last["tag"]:
                            continue
                    if last.get("Name"):
                        if last["Name"] not in gff:
                            gff[last["Name"]] = {}
                        if feature not in gff[last["Name"]]:
                            gff[last["Name"]][feature] = []
                        gff[last["Name"]][feature].append(fields)
                    elif last.get("Parent"):
                        if last["Parent"] not in gff:
                            gff[last["Parent"]] = {}
                        if feature not in gff[last["Parent"]]:
                            gff[last["Parent"]][feature] = []
                        gff[last["Parent"]][feature].append(fields)
                    elif last.get("ID"):
                        if last["ID"] not in gff:
                            gff[last["ID"]] = {}
                        if feature not in gff[last["ID"]]:
                            gff[last["ID"]][feature] = []
                        gff[last["ID"]][feature].append(fields)
                if line_counter % 1000 == 0 or line_counter == num_lines:
                    progress = make_progress_bar(line_counter, num_lines, t1, 70)
                    print("\r", progress[0] % progress[1:], end="", flush=True)
        elif format == "gtf":
            for line in f:
                line_counter += 1
                if line[0] != "#":
                    fields = line.rstrip().split("\t")
                    last = processGeneNameGTF(fields[8])
                    feature = fields[2]
                    if feature not in args.feat:
                        continue
                    # TODO add mane select filtering for gtf
                    if args.skipscaffoldcontigs:
                        chrom = fields[0] # TODO make sure this is the right field
                        if chrom.startswith("NT_") or chrom.startswith("NW_"):
                            continue
                    if last.get("transcript_id"):
                        if last["transcript_id"] not in gff:
                            gff[last["transcript_id"]] = {}
                        if feature not in gff[last["transcript_id"]]:
                            gff[last["transcript_id"]][feature] = []
                        gff[last["transcript_id"]][feature].append(fields)
                    elif last.get("gene_id"):
                        if last["gene_id"] not in gff:
                            gff[last["gene_id"]] = {}
                        if feature not in gff[last["gene_id"]]:
                            gff[last["gene_id"]][feature] = []
                        gff[last["gene_id"]][feature].append(fields)
    print("")  # for new line after progress bar
    return gff, intervals_chrom_format


def filterFeatureInGFF(gff, feat):
    """
    Keep/filters GFF records that include specified feature
    and returns a new GFF
    """
    filtered_gff = collections.defaultdict()
    for gene in gff.keys():
        if len(gff[gene][feat]) != 0:
            filtered_gff[gene] = gff[gene][feat]
    return filtered_gff


def getPloidy(vcf):
    var = [y for x, y in next(vcf.fetch()).samples.items()]
    p = [len(v.get("GT")) for v in var if v.get("GT")[0] is not None]
    return p[0] if p else 1


def getPhased(vcf):
    var = [y for x, y in next(vcf.fetch()).samples.items()]
    p = any([not v.phased for v in var])
    return not p if var else False


def revcomp(seq):
    tt = seq.maketrans("ACGT?N-", "TGCA?N-")
    return seq[::-1].translate(tt)


def make_progress_bar(rec, total, t1, width):
    i = (rec / total * 100) % 100
    if i != 0:
        plus = "+" * int(i * (width / 100))
        dots = "." * (width - int(i * width / 100))
    else:
        plus = "+" * width
        dots = ""
        i = 100
    t2 = time.time()
    elapsed = t2 - t1
    return "[" + plus + dots + "] " + "%5.2f%% %7.2f s", i, elapsed


def get_samples(vcf):
    # Function to get a list of sample names from the VCF file using pysam
    if hasattr(vcf, "header") and hasattr(vcf.header, "samples"):
        return list(vcf.header.samples)
    return []


def get_variants(vcf):
    # Function to get a list of variant positions from the VCF file using pysam
    variants = set()
    for rec in vcf.fetch():
        for alt in rec.alts:
            variants.add((rec.chrom, rec.pos, rec.ref, alt))
    return variants


def generate_header(intervals, region, feat, variants_inserted, varsites, vcf):
    header = intervals[region][feat][0][-1]
    if varsites > 0:
        header += "    "
        for var in variants_inserted:
            var_info = vcf.fetch(var[0], var[1] - 1, var[1])
            rs_id = ""
            ref_allele = ""
            alt_alleles = set()
            for v in var_info:
                rs_id = v.id if v.id is not None else "."
                ref_allele = v.ref
                alt_alleles.update(set(v.alts))
            header += f"{var[0]}:{var[1]}:{rs_id}:{ref_allele}:{','.join(alt_alleles)},"
        header = header.rstrip(",")
    return header


def main():
    # TODO skip reading GFF chromosomes not in VCF
    # TODO check that vcf is indexed before proceeding
    # TODO add option for treating VCF as if no samples are present (might need to write a brand new file for that to work, then bgzip and index it)

    parser = get_parser()
    args = parser.parse_args()
    if args.byvariant:
        args.addref = True

    # check that the files exist
    should_exit = False
    if not os.path.exists(args.vcf):
        print(f"VCF file '{args.vcf}' does not exist.")
    if not os.path.exists(args.fasta):
        print(f"FASTA file '{args.fasta}' does not exist.")
    if args.gff and not os.path.exists(args.gff):
        print(f"GFF file '{args.gff}' does not exist.")
    if args.bed and not os.path.exists(args.bed):
        print(f"BED file '{args.bed}' does not exist.")
    if should_exit:
        print("Exiting ...")
        sys.exit(1)

    if args.feat:
        args.feat = [f.strip() for f in args.feat.split(",")]

    print("Finished parsing arguments")

    # check that either gff or bed is provided, but not both
    if args.gff and args.bed:
        print("Only --gff or --bed are allowed. Exiting ...")
        sys.exit(parser.print_help())
    elif not args.gff and not args.bed:
        print("Either --gff or --bed are required. Exiting ...")
        sys.exit(parser.print_help())
    if args.feat == "": # TODO if user intentionally uses "", select all features instead of defaulting to "gene"
        print("No feature specified with --feat. Using 'gene' as default.")
        args.feat = ["gene"]

    print("Reading VCF file [", args.vcf, "] ... ", end="", sep="")
    vcf = pysam.VariantFile(args.vcf)
    vcf_chrom_format = helper.get_chrom_format(list(vcf.header.contigs)[0])
    print("done")
    if args.byvariant:
        samples = []
    else:
        samples = get_samples(vcf)
        # if no samples are present, create an "ALT" sample that is haploid and has the first alt allele for all variants
        if len(samples) == 0:
            print("No samples were found in VCF. Creating an 'ALT' sample.")
            samples.append("ALT")
    vcf_chroms = set(vcf.header.contigs)

     # TODO consider reading VCF first, then if args.skip is set, only read in chromosomes found in the vcf
    if args.bed:
        print("Reading BED file [", args.bed, "] ... ", end="", sep="")
        # TODO make sure bed files can skip scaffold/contig names if desired
        intervals = ReadBED(args.bed)
    else:
        # read GFF file
        print("Reading GFF file [", args.gff, "] ... ")
        # TODO make sure gff files can skip scaffold/contig names if desired
        intervals, intervals_chrom_format = ReadGFF(parser, vcf_chroms, vcf_chrom_format, args)
    print("done")
    # get the intervals chromosome format
    # intervals_chrom_format = helper.get_chrom_format(list(intervals.values())[0][0][0])
    intervals_chrom_format = helper.get_chrom_format(list(list(intervals.values())[0].values())[0][0][0])

    print("Reading FASTA reference file [", args.fasta, "] ... ", end="", sep="")
    ref = pysam.FastaFile(args.fasta)
    ref_chrom_format = helper.get_chrom_format(ref.references[0])
    print("done")

    # handling for if the chromosome formats between any of the files are incompatible
    if intervals_chrom_format != ref_chrom_format:
        # if the GFF/BED and FASTA chromosome formats don't match, the user must fix it
        print(f"Error: Chromosome format between the GFF/BED ('{intervals_chrom_format}') and reference FASTA ('{ref_chrom_format}') do not match.")
        print("Exiting ...")
        sys.exit(1)

    if args.byvariant:
        ploidy = 1
        phased = False
        print("By-variant mode selected. Setting ploidy to 1 and treating as unphased.")
    else:
        ploidy = getPloidy(vcf)
        print("Ploidy is:", ploidy)
        phased = getPhased(vcf)
        if not phased:
            print('No phased genotypes found on first variant. Treating as "unphased"')
        else:
            print('Phased genotypes found on first variant. Treating as "phased"')

    # if args.feat:
    #     outdir = args.out + "_" + args.feat
    # else:
    outdir = args.out.rstrip("/")

    if args.blend:
        print("Concatenating all [", args.feat, "]") # TODO write this better now that it is a list
    else:
        print("Writing all [", "intervals" if args.feat == "" else args.feat, "] separately") # TODO write this better now that it is a list
    print("Setting output directory to:", outdir)

    if not os.path.exists(outdir) or args.force:
        os.makedirs(outdir, exist_ok=True)
    else:
        proceed = input(
            outdir + " exists. Files will be replaced. Do you want to proceed? [y|n]: "
        )
        if not re.match("[Yy][EEs]*", proceed):
            print("Exiting ...")
            sys.exit(parser.print_help())

    # get region keys from GFF
    regions = list(intervals.keys())
    if len(regions) == 0:
        print("No regions found. Exiting ...")
        sys.exit(0)
    else:
        print("Total number of regions found:", len(regions))

    # check multiple records per feat
    if args.gff:
        single = [len(intervals[i]) == 1 for i in regions]
        if all(single):
            args.blend = True
            print("Found all regions with single records. Treating as --blend")

    # start counting time
    t1 = time.time()
    # start counter
    feature_counter = 0
    # count skipped regions
    withdata = 0

    variants_not_inserted = get_variants(vcf)
    num_total_variants = len(variants_not_inserted)

    for region in regions:
        for feat in intervals[region]:
            sequences, varsites, variants_inserted = getSequences(intervals, region, feat, ref, vcf, ploidy, phased, samples, ref_chrom_format, vcf_chrom_format, args)
            # remove variants_inserted from variants_not_inserted
            variants_not_inserted = variants_not_inserted - variants_inserted
            if args.skip and varsites != 0:
                withdata += 1
            if (args.skip and varsites != 0) or not args.skip:
                if args.byvariant: # TODO remove duplicate code
                    for var in variants_inserted:
                        header = generate_header(intervals, region, feat, {var}, 1, vcf)
                        seqData = processGeneNameGFF(header) if args.gff else processGeneNameGTF(header) # TODO make sure this works for bed files too
                        is_MANE = "tag" in seqData and "MANESelect" in seqData["tag"]
                        geneName = seqData["gene"] if "gene" in seqData else region
                        for featname in sequences.keys(): # TODO add the ability to put all sequence pairs in the same file
                            base_dir = f"{outdir}/{geneName}-{var[0]}-{var[1]}"
                            if not os.path.exists(base_dir):
                                os.makedirs(base_dir, exist_ok=True)
                            output_dir = f"{base_dir}/{feat}-{featname}-{var[0]}-{var[1]}.fa"
                            if is_MANE:
                                output_dir = f"{base_dir}/{feat}-{featname}-{var[0]}-{var[1]}-MANE.fa"
                            with open(output_dir, "w") as out:
                                printFasta(sequences[featname], header, out, var)
                else:
                    header = generate_header(intervals, region, feat, variants_inserted, varsites, vcf)
                    seqData = processGeneNameGFF(header) if args.gff else processGeneNameGTF(header) # TODO make sure this works for bed files too
                    is_MANE = "tag" in seqData and "MANESelect" in seqData["tag"]
                    geneName = seqData["gene"] if "gene" in seqData else region
                    for featname in sequences.keys(): # TODO add the ability to put all sequence pairs in the same file
                        base_dir = f"{outdir}/{geneName}"
                        if not os.path.exists(base_dir):
                            os.makedirs(base_dir, exist_ok=True)
                        output_dir = f"{base_dir}/{feat}-{featname}.fa"
                        if is_MANE:
                            output_dir = f"{base_dir}/{feat}-{featname}-MANE.fa"
                        with open(output_dir, "w") as out:
                            printFasta(sequences[featname], header, out)
        feature_counter += 1
        progress = make_progress_bar(feature_counter, len(regions), t1, 70)
        print("\r", progress[0] % progress[1:], end="", flush=True)
    print("")

    if len(variants_not_inserted) > 0:
        with open(f"{outdir}/variants_not_inserted.txt", "w") as outvar:
            outvar.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for chrom, pos in sorted(variants_not_inserted):
                infos = vcf.fetch(chrom, pos - 1, pos)
                for info in infos:
                    chrom = info.chrom
                    pos = info.pos
                    rs_id = info.id if info.id is not None else "."
                    ref = info.ref
                    alts = ",".join(info.alts) if info.alts is not None else "."
                    qual = str(info.qual) if info.qual is not None else "."
                    filt = ";".join(info.filter.keys()) if info.filter.keys() else "PASS"
                    info_str = ";".join([f"{k}={v}" for k, v in info.info.items()]) if info.info else "."
                    outvar.write(f"{chrom}\t{pos}\t{rs_id}\t{ref}\t{alts}\t{qual}\t{filt}\t{info_str}\n")
        print(f"Warning: {len(variants_not_inserted)} of {num_total_variants} variants were not inserted into any sequence because they did not overlap with any relevant annotated region.")
        print("They were written to 'variants_not_inserted.txt' in the output directory.")
    else:
        print(f"All {num_total_variants} variants were inserted into sequences.")
    
    if args.skip:
        # TODO fix this when byvariant option is used- currently counts features, not files written
        print(f"Skipped {feature_counter - withdata} regions with no variants. Wrote {withdata} regions with variants.")


if __name__ == "__main__":
    main()
