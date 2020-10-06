# Download phased data from 1000 Genomes:
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

# Filter SNPs in CYP2C19 + 10 Mb up/down:
out="CYP2C19_10Mb"
bcftools filter -r 10:91522463-101612671 ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -o ${out}.vcf && bgzip ${out}.vcf && rm ${out}.vcf && bcftools index -t ${out}.vcf.gz

# Split multiallelic records:
input="CYP2C19_10Mb"
out1=${input}_split
bcftools norm -m -any ${input}.vcf.gz > ${out1}.vcf && bgzip ${out1}.vcf && rm ${out1}.vcf && bcftools index -t ${out1}.vcf.gz

# Calculate the INFO field:
out2=${out1}_info
bcftools +fill-tags ${out1}.vcf.gz -o ${out2}.vcf && bgzip ${out2}.vcf && rm ${out2}.vcf && bcftools index -t ${out2}.vcf.gz

# Convert VCF to HAP/LEGEND/SAMPLE format:
bcftools convert -h ${input}_ALL --vcf-ids ${out2}.vcf.gz

