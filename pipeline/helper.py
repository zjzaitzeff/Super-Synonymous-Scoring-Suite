# This file contains helper functions related to converting chromosome formats

import re

NCBI_TO_SINGLE_LETTER = {"NC_000001.11": "1", "NC_000002.12": "2", "NC_000003.12": "3", "NC_000004.12": "4", "NC_000005.10": "5", "NC_000006.12": "6", "NC_000007.14": "7", "NC_000008.11": "8", "NC_000009.12": "9", "NC_000010.11": "10", "NC_000011.10": "11", "NC_000012.12": "12", "NC_000013.11": "13", "NC_000014.9": "14", "NC_000015.10": "15", "NC_000016.10": "16", "NC_000017.11": "17", "NC_000018.10": "18", "NC_000019.10": "19", "NC_000020.11": "20", "NC_000021.9": "21", "NC_000022.11": "22", "NC_000023.11": "X", "NC_000024.10": "Y", "NC_012920.1": "12920"}
SINGLE_LETTER_TO_NCBI = {v: k for k, v in NCBI_TO_SINGLE_LETTER.items()}


def get_chrom_format(chrom):
    """
    Returns "chr", "single_letter", "ncbi", or "unknown"
    """ 
    if chrom.startswith("chr"):
        return "chr"
    elif re.match(r"^[0-9XYM]+$", chrom):
        return "single_letter"
    elif chrom in NCBI_TO_SINGLE_LETTER:
        return "ncbi"
    else:
        return "unknown"


def is_chrom_valid(chrom):
    chrom = chrom.lower()
    if chrom.startswith('chr'):
        chrom = chrom[3:]
    if chrom in [str(i) for i in range(1, 23)] + ['x', 'y', 'm', 'mt']:
        return True
    # check ncbi format
    if chrom.startswith('nc_') or chrom.startswith('nt_') or chrom.startswith('nw_'):
        return True
    return False


def convert_chrom_format(chrom, from_format, to_format):
    """
    Converts the chromosome format of the given chromosome from one format to another.
    Supported formats: 'chr', 'single_letter', 'ncbi',
    """
    if from_format == to_format:
        return chrom  # no conversion needed
    if from_format not in ["chr", "single_letter", "ncbi"] or to_format not in ["chr", "single_letter", "ncbi"]:
        return chrom  # unsupported format, return original
    
    if from_format == "chr":
        chrom = chrom[3:]  # remove 'chr' prefix for easier mapping
        from_format = "single_letter"  # treat as single_letter for mapping
    if from_format == "single_letter":
        if to_format == "chr":
            return "chr" + chrom
        elif to_format == "ncbi" and chrom in SINGLE_LETTER_TO_NCBI:
            return SINGLE_LETTER_TO_NCBI[chrom]
        else:
            return chrom  # no mapping found, return original
    if from_format == "ncbi":
        if to_format == "single_letter":
            if chrom in NCBI_TO_SINGLE_LETTER:
                return NCBI_TO_SINGLE_LETTER[chrom]
            else:
                return chrom  # no mapping found, return original
        elif to_format == "chr":
            if chrom in NCBI_TO_SINGLE_LETTER:
                return "chr" + NCBI_TO_SINGLE_LETTER[chrom]
            else:
                return chrom  # no mapping found, return original
    return chrom  # fallback, return original
