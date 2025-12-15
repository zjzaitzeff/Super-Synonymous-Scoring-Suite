# This script analyzes the master summary data and generates a report.

import sys
import numpy as np
import matplotlib.pyplot as plt

def get_number(value):
    try:
        return float(value)
    except ValueError:
        return value

def is_number(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

def read_file(file_path):
    variant_to_data = {}
    with open(file_path) as inf:
        header = inf.readline()
        for line in inf:
            items = line.strip().split("\t")
            variant_id = items[0]
            chrom = items[1]
            pos = items[2]
            rsID = items[3]
            ref = items[4]
            alt = items[5]
            info = items[6]
            targetscan = items[7]
            pita = items[8]
            miranda = items[9]
            miRNA_2_tool_agreement = int(items[10])
            miRNA_3_tool_agreement = int(items[11])
            sei = get_number(items[12])
            spliceai = get_number(items[13])
            vienna_scores = items[14].split("|")
            vienna_scores = [float(x[x.index("=")+1:]) for x in vienna_scores if x != "NA"]
            codon_efficiency = get_number(items[15])
            extramp = items[16]
            clinvar_sig = items[17]
            variant_to_data[variant_id] = {
                "chrom": chrom,
                "pos": pos,
                "rsID": rsID,
                "ref": ref,
                "alt": alt,
                "info": info,
                "targetscan": targetscan,
                "pita": pita,
                "miranda": miranda,
                "miRNA_2_tool_agreement": miRNA_2_tool_agreement,
                "miRNA_3_tool_agreement": miRNA_3_tool_agreement,
                "sei": sei,
                "spliceai": spliceai,
                "vienna_scores": vienna_scores,
                "codon_efficiency": codon_efficiency,
                "extramp": extramp
            }

    return variant_to_data

def get_cutoff(values):
    mean_val = np.mean(values)
    std_val = np.std(values)
    cutoff_high = mean_val + 2.7*std_val
    cutoff_low = mean_val - 2.7*std_val
    return cutoff_low, cutoff_high

def get_vienna_data(variant_to_data):
    delta_mfes = []
    delta_efes = []
    delta_cfes = []
    delta_cdes = []
    delta_meafes = []
    delta_ends = []
    mfeeds = []
    cfeeds = []
    meads = []
    efeeds = []
    for variant, data in variant_to_data.items():
        delta_mfes.append(data["vienna_scores"][0])
        delta_efes.append(data["vienna_scores"][1])
        delta_cfes.append(data["vienna_scores"][2])
        delta_cdes.append(data["vienna_scores"][3])
        delta_meafes.append(data["vienna_scores"][4])
        delta_ends.append(data["vienna_scores"][5])
        mfeeds.append(data["vienna_scores"][6])
        cfeeds.append(data["vienna_scores"][7])
        meads.append(data["vienna_scores"][8])
        efeeds.append(data["vienna_scores"][9])
    return delta_mfes, delta_efes, delta_cfes, delta_cdes, delta_meafes, delta_ends, mfeeds, cfeeds, meads, efeeds

def get_vienna_cutoffs(delta_mfes, delta_efes, delta_cfes, delta_cdes, delta_meafes, delta_ends, mfeeds, cfeeds, meads, efeeds):
    delta_mfe_cutoffs = get_cutoff(delta_mfes)
    delta_efe_cutoffs = get_cutoff(delta_efes)
    delta_cfe_cutoffs = get_cutoff(delta_cfes)
    delta_cde_cutoffs = get_cutoff(delta_cdes)
    delta_meafe_cutoffs = get_cutoff(delta_meafes)
    delta_end_cutoffs = get_cutoff(delta_ends)
    mfeed_cutoffs = get_cutoff(mfeeds)
    cfeed_cutoffs = get_cutoff(cfeeds)
    mead_cutoffs = get_cutoff(meads)
    efeed_cutoffs = get_cutoff(efeeds)

    return delta_mfe_cutoffs, delta_efe_cutoffs, delta_cfe_cutoffs, delta_cde_cutoffs, delta_meafe_cutoffs, delta_end_cutoffs, mfeed_cutoffs, cfeed_cutoffs, mead_cutoffs, efeed_cutoffs

def is_outlier(vienna_scores, delta_mfe_cutoffs, delta_efe_cutoffs, delta_cfe_cutoffs, delta_cde_cutoffs, delta_meafe_cutoffs, delta_end_cutoffs, mfeed_cutoffs, cfeed_cutoffs, mead_cutoffs, efeed_cutoffs):
    delta_mfe, delta_efe, delta_cfe, delta_cde, delta_meafe, delta_end, mfeed, cfeed, mead, efeed = vienna_scores
    if (delta_mfe < delta_mfe_cutoffs[0] or delta_mfe > delta_mfe_cutoffs[1] or
        delta_efe < delta_efe_cutoffs[0] or delta_efe > delta_efe_cutoffs[1] or
        delta_cfe < delta_cfe_cutoffs[0] or delta_cfe > delta_cfe_cutoffs[1] or
        delta_cde < delta_cde_cutoffs[0] or delta_cde > delta_cde_cutoffs[1] or
        delta_meafe < delta_meafe_cutoffs[0] or delta_meafe > delta_meafe_cutoffs[1] or
        delta_end < delta_end_cutoffs[0] or delta_end > delta_end_cutoffs[1] or
        mfeed > mfeed_cutoffs[1] or
        cfeed > cfeed_cutoffs[1] or
        mead > mead_cutoffs[1] or
        efeed > efeed_cutoffs[1]):
        return True
    return False

def plot_histogram(data, title, xlabel, ylabel, output_file):
    plt.hist(data, bins=50, color='blue', alpha=0.7)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(axis='y', alpha=0.75)
    plt.savefig(output_file)
    plt.close()

def plot_vienna_histograms(delta_mfes, delta_efes, delta_cfes, delta_cdes, delta_meafes, delta_ends, mfeeds, cfeeds, meads, efeeds):
    plot_histogram(delta_mfes, 'Distribution of Delta MFE Values', 'Delta MFE (kcal/mol)', 'Number of Variants', 'delta_mfe_distribution.png')
    plot_histogram(delta_efes, 'Distribution of Delta EFE Values', 'Delta EFE (kcal/mol)', 'Number of Variants', 'delta_efe_distribution.png')
    plot_histogram(delta_cfes, 'Distribution of Delta CFE Values', 'Delta CFE (kcal/mol)', 'Number of Variants', 'delta_cfe_distribution.png')
    plot_histogram(delta_cdes, 'Distribution of Delta CDE Values', 'Delta CDE (kcal/mol)', 'Number of Variants', 'delta_cde_distribution.png')
    plot_histogram(delta_meafes, 'Distribution of Delta MEAFE Values', 'Delta MEAFE (kcal/mol)', 'Number of Variants', 'delta_meafe_distribution.png')
    plot_histogram(delta_ends, 'Distribution of Delta END Values', 'Delta END (kcal/mol)', 'Number of Variants', 'delta_end_distribution.png')
    plot_histogram(mfeeds, 'Distribution of MFEED Values', 'MFEED (kcal/mol)', 'Number of Variants', 'mfeed_distribution.png')
    plot_histogram(cfeeds, 'Distribution of CFEED Values', 'CFEED (kcal/mol)', 'Number of Variants', 'cfeed_distribution.png')
    plot_histogram(meads, 'Distribution of MEAD Values', 'MEAD (kcal/mol)', 'Number of Variants', 'mead_distribution.png')
    plot_histogram(efeeds, 'Distribution of EFEED Values', 'EFEED (kcal/mol)', 'Number of Variants', 'efeed_distribution.png')

if __name__ == "__main__":
    # input_file = sys.argv[1]
    input_file = "outputs_jag1/master_summary.tsv"

    variant_to_data = read_file(input_file)

    total_variants = 0
    explained_by_miRNA = 0
    explained_by_splicing = 0
    explained_by_tfbs = 0
    explained_by_rna_structure = 0
    explained_by_codon_efficiency = 0
    explained_by_extramp = 0
    not_explained = 0
    multiple_explanations = 0

    delta_mfes, delta_efes, delta_cfes, delta_cdes, delta_meafes, delta_ends, mfeeds, cfeeds, meads, efeeds = get_vienna_data(variant_to_data)
    delta_mfe_cutoffs, delta_efe_cutoffs, delta_cfe_cutoffs, delta_cde_cutoffs, delta_meafe_cutoffs, delta_end_cutoffs, mfeed_cutoffs, cfeed_cutoffs, mead_cutoffs, efeed_cutoffs = get_vienna_cutoffs(delta_mfes, delta_efes, delta_cfes, delta_cdes, delta_meafes, delta_ends, mfeeds, cfeeds, meads, efeeds)
    # plot_vienna_histograms(delta_mfes, delta_efes, delta_cfes, delta_cdes, delta_meafes, delta_ends, mfeeds, cfeeds, meads, efeeds)
    
    for variant, data in variant_to_data.items():
        total_variants += 1
        targetscan = data["targetscan"]
        pita = data["pita"]
        miranda = data["miranda"]
        miRNA_2_tool_agreement = data["miRNA_2_tool_agreement"]
        miRNA_3_tool_agreement = data["miRNA_3_tool_agreement"]
        sei = data["sei"]
        spliceai = data["spliceai"]
        vienna_scores = data["vienna_scores"]
        codon_efficiency = data["codon_efficiency"]
        extramp = data["extramp"]
        explained = False
        if is_number(sei):
            explained_by_tfbs += 1
            explained = True
        if miRNA_3_tool_agreement >= 1:
            explained_by_miRNA += 1
            explained = True
        if is_number(spliceai):
            explained_by_splicing += 1
            explained = True
        if is_outlier(vienna_scores, delta_mfe_cutoffs, delta_efe_cutoffs, delta_cfe_cutoffs, delta_cde_cutoffs, delta_meafe_cutoffs, delta_end_cutoffs, mfeed_cutoffs, cfeed_cutoffs, mead_cutoffs, efeed_cutoffs):
            explained_by_rna_structure += 1
            explained = True
        if is_number(codon_efficiency) and abs(codon_efficiency) >= 0.5: # TODO arbitrary cutoff for now
            explained_by_codon_efficiency += 1
            explained = True
        if extramp in ["Yes", "ramp gain", "ramp loss", "ramp length change"]:
            explained_by_extramp += 1
            explained = True
        if not explained:
            not_explained += 1

    total_variants = len(variant_to_data)
    print(f"Total variants: {total_variants}")
    print(f"miRNA binding: {explained_by_miRNA} ({explained_by_miRNA/total_variants:.2%})")
    print(f"Splicing: {explained_by_splicing} ({explained_by_splicing/total_variants:.2%})")
    print(f"TFBS: {explained_by_tfbs} ({explained_by_tfbs/total_variants:.2%})")
    print(f"mRNA secondary structure: {explained_by_rna_structure} ({explained_by_rna_structure/total_variants:.2%})")
    print(f"Codon efficiency: {explained_by_codon_efficiency} ({explained_by_codon_efficiency/total_variants:.2%})")
    print(f"Ramps: {explained_by_extramp} ({explained_by_extramp/total_variants:.2%})")
    print(f"No explanation: {not_explained} ({not_explained/total_variants:.2%})")

    # plot the sorted counts, but no explanation last
    labels = ['Splicing', 'TFBS', 'microRNA', 'mRNA structure', 'Codon efficiency', 'Ramps', 'No explanation']
    counts = [explained_by_splicing, explained_by_tfbs, explained_by_miRNA, explained_by_rna_structure, explained_by_codon_efficiency, explained_by_extramp, not_explained]
    x = np.arange(len(labels))
    plt.bar(x, counts, color='#4285f4')
    # plot no explanation bar in gray
    plt.bar(x[-1], counts[-1], color='#888888')
    # plot counts and percents on top of bars
    for i, v in enumerate(counts):
        plt.text(i, v + max(counts)*0.01, f"{v} ({v/total_variants:.1%})", ha='center', va='bottom', fontsize=10)
    plt.xticks(x, labels, rotation=45, ha='right', fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylim(0, max(counts)*1.1)
    # no title
    plt.ylabel('Number of GWAS Variants', fontsize=14)
    plt.tight_layout()
    plt.savefig('variant_explanations.png', dpi=300)
    plt.close()

    # # get outlier cutoffs in delta mfe values (2 standard deviations from mean)
    # import numpy as np
    # mean_mfe = np.mean(delta_mfes)
    # std_mfe = np.std(delta_mfes)
    # cutoff_high = mean_mfe + 2*std_mfe
    # cutoff_low = mean_mfe - 2*std_mfe

    # print(f"Mean delta MFE: {mean_mfe:.2f}, Std Dev: {std_mfe:.2f}")
    # print(f"Outlier cutoffs: < {cutoff_low:.2f} or > {cutoff_high:.2f}")

    # histogram of delta mfe values
    # import matplotlib.pyplot as plt
    # plt.hist(delta_mfes, bins=50, color='blue', alpha=0.7)
    # plt.xlabel('Delta MFE (kcal/mol)')
    # plt.ylabel('Number of Variants')
    # plt.title('Distribution of Delta MFE Values from ViennaRNA')
    # plt.grid(axis='y', alpha=0.75)
    # fig_path = 'delta_mfe_distribution.png'
    # plt.savefig(fig_path)
    # print(f"Delta MFE distribution plot saved to {fig_path}")