# gets the overlap between our intput VCF and clinvar

import pysam

if __name__ == "__main__":
    input_vcf = "inputs/syn_variants.vcf.gz"
    clinvar_vcf = "ref/clinvar.vcf.gz"
    output_vcf = "outputs/clinvar_overlap.vcf"

    input_vcf = pysam.VariantFile(input_vcf)
    clinvar_vcf = pysam.VariantFile(clinvar_vcf)
    output_vcf = pysam.VariantFile(output_vcf, "w", header=clinvar_vcf.header)

    non_overlap_count = 0
    overlap_count = 0
    significance_counts = {}

    for record in input_vcf.fetch():
        chrom = record.chrom
        pos = record.pos
        ref = record.ref
        alt = record.alts[0]

        clinvar_records = clinvar_vcf.fetch(chrom, pos - 1, pos)
        clinvar_records = [clinvar_record for clinvar_record in clinvar_records]
        if len(clinvar_records) == 0:
            non_overlap_count += 1
        else:
            overlap_count += 1
            found = False
            for clinvar_record in clinvar_records:
                if clinvar_record.pos == pos and clinvar_record.ref == ref and clinvar_record.alts is not None and alt in clinvar_record.alts:
                    output_vcf.write(clinvar_record)
                    sig = clinvar_record.info["CLNSIG"]
                    sig = ",".join(sig) if isinstance(sig, (list, tuple)) else str(sig)
                    if sig not in significance_counts:
                        significance_counts[sig] = 0
                    significance_counts[sig] += 1
                    found = True
                    break
            if not found:
                non_overlap_count += 1
            
    print(f"Number of overlapping variants: {overlap_count}")
    print(f"Number of non-overlapping variants: {non_overlap_count}")

    print("Significance counts:")
    for sig, count in significance_counts.items():
        print(f"{sig}: {count}")