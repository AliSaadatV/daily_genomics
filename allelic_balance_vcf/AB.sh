#Add allelic balance (AB) also known as variant allele fraction (VAF)
bcftools +fill-tags input.vcf -Oz -o out.vcf.gz -- -t FORMAT/VAF

#Filter the het-genotype of samples with AB < 0.2 or AB > 0.8
gatk VariantFiltration \
        -V out.vcf.gz \
        -R $REF \
        --genotype-filter-expression "isHet == 1 && VAF < 0.2" \
        --genotype-filter-name "HetLowAB" \
        --genotype-filter-expression "isHet == 1 && VAF > 0.8" \
        --genotype-filter-name "HetHighAB" \
        -O out_AB.vcf.gz
