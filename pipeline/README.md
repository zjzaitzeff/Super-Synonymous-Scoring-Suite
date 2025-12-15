# TODO
1. Clean up the folder
    1. Move unused files to the "old" folder
    2. Rename files so it's clearer which ones to run first, second, etc.
    3. Document the pipeline (I've written my process here, but we need a more generic process for the user, ideally with as few steps as possible)
2. Consider Enformer instead of Sei
3. Main pipeline
    1. Put everything into one file so users only have to run one command (vcf2fasta, prep_single_variant_vcfs.py, and bash get_data_pipeline_copy.sh)
    2. Reorganize the main pipeline's inputs to use options so order matters less, or consider if the order is good right now. Remove default options.
    3. Remove commented out code and add any necessary comments
    4. Test running without slurm option- might run each variant sequentially
    5. Better text at end (not just All done!)
    6. Make it more obvious to run the master summary script after the main pipeline finishes
4. Print the value determined by Sei instead of interpreting it before putting it in the summary.tsv
5. Consolidate required reference files
    1. duplicate ref FASTAs used by Sei, SpliceAI, and vcf2fasta- can we have them all use the vcf2fasta one?
    2. We don't need most of the reference files (CDS, RNA)
6. DOCKER :|
7. What options will the user want that we haven't made available?
8. Revisit the wij file- we tried using gTAI, but the values were REALLY weird. Is there a better file we should use? Do we want to do any tissue/cell-specific stuff?
9. Implement and test the features listed as TODO throughout the pipeline
    1. vcf2fasta
        1. Make sure all features/options are also implemented for GTFs and BED files
        2. Add an "ignore samples" option that creates a new VCF with all samples removed. vcf2fasta then can treat this file as REF ALT
        3. Handle multiple alts in VCF (ex: G,T). We are currently only testing the first one!
        4. If --mane is used, but the gff doesn't use the MANEselect tag, ignore that option
        5. Check that the VCF is indexed before proceeding
        6. Optimize:
            1. For the GFF, skip reading chromosomes not in the VCF (making sure to check different chromosome formats)
            2. Remove getSequences prep
            3. Remove duplicate code
        7. Less critical:
            1. If user intentionally uses "" for feat, select all features instead of defaulting to "gene"
            2. Fix print statements that write out lists that previously weren't lists
            3. Add ability to put all sequence pairs in the same file
    2. write_summary.py
        1. Write largest spliceAI score and make it clear if donor/accepter and loss/gain
    2. master_summary.py
        1. Use sys.argv instead
        2. Find a better fix for when GFF and input VCF genes are not the same
    3. analyze_master_summary.py
        1. Use sys.argv instead
        2. Need to explain my cutoffs:
            1. codon_efficiency (+-0.5)
            2. vienna (outliers)- need outliers for non-GWAS SNPs for this to be better
            3. Sei (>=0.18)- was this in their paper?
    4. write_summary.py
        1. Rename functions
        2. Print "delta" instead of the delta sign since the sign can break some text editors
        3. Reorganize functions
        4. Revisit the PLX mechanic
        5. "Hardcode" the paths to output files instead of using a for loop
10. Report effect size for GWAS synonymous variants tested
11. Can we do fine mapping on the GWAS synonymous variants before running them through out pipeline?
12. Adapt the pipeline to work for other species too
	1. SpliceAI needs to run on sequences when the species is not human
	2. Sei needs species-specific reference genome
	3. What do we do about the miRNA motif file?



# NOTES FOR THE PROGRAMMERS TO BE MOVED TO THE PAPER

# GWAS Catalog Synonymous Variants Test
## Get VCF
[get_syn.py](./get_syn.py) creates a [VCF of synonymous variants](./inputs/syn_variants.vcf) found in the [GWAS catalog](./inputs/gwas_catalog_v1.0.2-associations_e115_r2025-09-15.tsv). It requires pysam to be installed (see OpenSSL error below if pysam doesn't work at first). For each synonymous variant, the file's info column has the risk allele, mapped gene, and reported and mapped traits. The ref and alt alleles are determined using pysam using a dbSNP file since the ref allele is sometimes the risk allele. The [dbSNP file](./ref/GCF_000001405.40.gz) and [index file](./ref/GCF_000001405.40.gz.tbi) were copied from the ramps and snps project. They are build 157, downloaded from [here](https://ftp.ncbi.nih.gov/snp/latest_release/VCF/) on 25/07/02.
Note that more synonymous SNPs are present in the GWAS catalog than are in the output file, but many were filtered out because the risk allele was "?", or because they were part of multi-allelic associations.
2128 synonymous SNPs were written to the [output file](./inputs/syn_variants.vcf).
```
python get_syn.py
```
## BGzip and Index
The VCF was subsequently zipped and indexed using the following command, which requires bcftools to be installed.
```
bash bgzip_and_index.sh ./inputs/syn_variants.vcf
```
## VCF to FASTA
[vcf2fasta.py](./vcf2fasta.py) is a heavily modified script based on [vcf2fasta](https://github.com/santiagosnchez/vcf2fasta) by Santiago Sanchez. It enables us to convert the GWAS catalog VCF to the gene, mRNA, and CDS FASTAs necessary to investigate synonymous regulatory signals. The following command creates all the FASTA files, organized into folders by variant.
```
python vcf2fasta.py --vcf "inputs/syn_variants.vcf.gz" --fasta "ref/GCF_000001405.40_GRCh38.p14_genomic.fna" --gff "ref/GCF_000001405.40_genomic.gff" --out "outputs/inputs" --feat "mRNA,CDS" --addref --skip --force --skipscaffoldcontigs --blend --mane --byvariant
```
Notice that I have added quite a few custom options:

    --skipscaffoldcontigs ignores any regions that start with NW or NT, helpful when genes occur in chromosomes and scaffolds/contigs and you just want the chromosomal gene sequence.
    --mane writes out only MANE select sequences, good for analyzing only one isoform per gene to speed up investigations. Unfortunately, synonymous variants are sometimes functional in non MANE select isoforms, so this speed up might miss relevant information
    --byvariant ensures that each variant is written out into it's own folder instead of all variants in the VCF in the same gene being written out together.

Note: the pipeline works best when the VCF doesn't have any samples. I may resolve this in the future.

vcf2fasta resulted in 2417 folders to test separately, representing 2119 variants. Some variants are found in multiple genes and 9 weren't inserted in any gene because they were in pseudogenes.

## Get Single Variant VCFs
Some predictors use VCFs as input (SpliceAI and Sei), so we write single variant VCFs into each of the output folders created by vcf2fasta using the following command:
```
python prep_single_variant_vcfs.py ./inputs/syn_variants.vcf ./outputs/inputs
```

## Run Prediction Pipeline on All Variants
Finally, we can run the pipeline on all of the folders produced by vcf2fasta! Here is the command we used:
```
sbatch get_data_pipeline_copy.sh inputs/syn_variants.vcf outputs/inputs outputs ref/miR_Family_Info_Human.fa 500 True 00:30:00
```
A slurm process is spawned for each of the variants and, when each finishes, a summary.tsv files is generated within the variant's output folder.

## Summarize Results
[master_summary.py](./master_summary.py) combines all the summary.tsv files for each of the variants and puts the results in [master_summary.tsv](./outputs/master_summary.tsv)
```
python master_summary.py
```
The master summary can then be analyzed/interpreted by [analyze_master_summary.py](./analyze_master_summary.py). It plots vienna histograms for each vienna metric (currently commented out), prints out the number/percent of variants explained by each regulatory signal, and creates a [histogram](./variant_explanations.png) of the number of variants explained by each regulatory signal. Note that there is heavy overlap between each of the signals because many variants are explained by multiple signals. miRNAs predict 63% of the variants tested.
```
python analyze_master_summary.py
```

# pip install
## sei
pip install docopt selene-sdk torch torchvision
## SPOT-RNA
pip install pandas numpy tqdm argparse six

## Python Version
We are required to use Python version 3.12 or less for Tensorflow to work, so we used Python 3.12.

## OpenSSL Error
For SpliceAI, pysam triggers an OpenSSL error on the supercomputer due to its prebuilt wheel. The issue is documented [here](https://github.com/pysam-developers/pysam/issues/1276).
To get around this issue, pysam should be installed from source instead like so:
```
# uninstall any current versions of pysam
pip uninstall pysam spliceai
# install pysam from binary
pip install --no-binary=pysam pysam
pip install spliceai
```

## Reference Genome
The [hg38 reference genome](./reference/hg38.fa.gz) was downloaded on 25/06/20 [here](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz) as suggested by the SpliceAI GitHub. We may switch it to a newer NCBI release at some point if needed.

## Sei Cuttoff
We needed to select a cutoff score for Sei to know when to report that a variant caused a significant regulatory change.
We decided on a cutoff diff score of 0.18, based on the following methods:
In the Sei paper, they reported the scores of 854 deleterious variants from HGMD. The distribution of variant max diff scores was heavily right skewed, so we calculated Q1, median, and Q3 for those scores which were:

Q1	median	Q3
0.183151209	0.385322691	0.879895552

We decided on a cutoff score of 0.18, which is quartile 1 of Sei's dataset of 854 deleterious variants they reported on in their paper. This means that in our tool, we report "yes" if a variant's max diff score greater than or equal to 0.18 and "no" if the variant's score is less than 0.18.