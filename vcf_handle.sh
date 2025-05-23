#!/bin/bash/

WDIR="/opt/home/jong/oryza/asm/pggb_pangenomes/output_all_O_mer_p90/"
VCF_WDIR="$WDIR/vcf_clean"
vcf_txt="$WDIR/vcf_gz.txt"

###################
# Merge Chromosomes 
###################
# since I am dealing with whole-genome assemblies, no need for +fixploidy
# no need for merging biallelic records to multiallelic records as well as I need biallelic sites
# prepare for cleaning the decomposed files
# if dir not existing

if [ ! -d "$VCF_WDIR" ]; then
	mkdir -p "$VCF_WDIR"
fi

cd $VCF_WDIR

# list the paths of decomposed vcf files and save in vcf_gz.txt
# concatenate VCF files
bcftools concat -f $vcf_txt | bgzip --threads 16 -c > pggb_merged.vcf.gz

# sort
bcftools sort pggb_merged.vcf.gz | bcftools +fill-tags -- -t AN,AC,AF,F_MISSING | bcftools view -o pggb_merged_final.vcf.gz -Oz

# remove "AT" 
bcftools annotate --remove INFO/AT pggb_merged_final.vcf.gz -o pggb_merged_final_clean.vcf.gz -Oz --threads 30

# delete intermediate VCFs
rm -f pggb_merged_final.vcf.gz
rm -f pggb_merged.vcf.gz