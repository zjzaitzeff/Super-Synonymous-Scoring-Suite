# the vcf file to bgzip and index
vcf_file=$1
fai_file=$2

if [ -z "$vcf_file" ]; then
  echo "Usage: $0 <vcf_file>"
  exit 1
fi

# bgzip the vcf file
bgzip -c $vcf_file > ${vcf_file}.gz
complete_code=$?
if [ $complete_code -ne 0 ]; then
  echo "Error: bgzip failed with exit code $complete_code"
  exit $complete_code
fi
# index the bgzipped vcf file
bcftools index -t ${vcf_file}.gz
complete_code=$?
if [ $complete_code -ne 0 ]; then
  echo "Error: bcftools index failed with exit code $complete_code"
  exit $complete_code
fi

echo "bgzipped and indexed '${vcf_file}' to '${vcf_file}.gz' and '${vcf_file}.gz.tbi'"
