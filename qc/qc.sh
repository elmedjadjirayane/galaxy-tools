#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

# Variables
vcf_file="$1"
pheno_file="$2"
maf_threshold="$3"
hwe_threshold="$4"
mgcf_threshold="$5"
max_heterozygosity="$6"
window_size="$7"
step_size="$8"
ld_threshold="$9"
pbed_output="${10}"
__tool_directory__="${11}"
perform_ld_pruning="${12}" # Boolean to control steps 10 and 11
fid="${13}"
pLD="&{14}"
echo "this what we have inside pLD" "$pLD"
# Temporary directory
TMPDIR=$(mktemp -d)

echo "Step 1: Prepare phenotype file"
if [ "$fid" == "false" ]; then
    echo "Adding a second column"
    awk 'BEGIN{OFS="\t"} {print $1, $1, substr($0, index($0, $2))}' "$pheno_file" > "$TMPDIR/temp_pheno"
else
    echo "No changes made to phenotype file"
    cp "$pheno_file" "$TMPDIR/temp_pheno"
fi
echo "check phenotype file"
cat "$TMPDIR/temp_pheno"
echo "Step 2: Create output directories"
mkdir -p ./plink2_output


# Step 3: Convert VCF to PLINK binary format
plink --allow-extra-chr --allow-no-sex --vcf "$vcf_file" --make-bed --double-id --pheno "$TMPDIR/temp_pheno" --out ./plink2_output/output || {
    echo "Error: Failed to convert VCF to PLINK format."; exit 1;
}

echo "Step 4: Standardize variant IDs"
plink2 --bfile ./plink2_output/output --allow-no-sex --set-all-var-ids @:# --make-bed --out ./plink2_output/output_renamed || {
    echo "Error: Failed to standardize variant IDs."; exit 1;
}

echo " Step 5: Filter by MAF and HWE thresholds"
plink2 --bfile ./plink2_output/output_renamed --allow-no-sex --maf "$maf_threshold" --hwe "$hwe_threshold" --make-bed --out ./plink2_output/output_filtered 2>&1 || {
    echo "Error: Failed to filter data by MAF and HWE."; exit 1;
}

echo " Step 6: Calculate allele frequencies"
plink --bfile ./plink2_output/output_filtered --allow-no-sex --freqx --out ./plink2_output/output_frequencies || {
    echo "Error: Failed to calculate allele frequencies."; exit 1;
}

echo " Step 7: Run quality control script"
python "$__tool_directory__/qc.py" "$mgcf_threshold" "$max_heterozygosity" || {
    echo "Error: QC script failed."; exit 1;
}

echo "Step 8: Filter SNPs based on MGCF"
plink2 --bfile ./plink2_output/output_filtered --allow-no-sex --extract ./plink2_output/mgcf_filtered_snps.txt --make-bed --out ./plink2_output/output_filtered_mgcf || {
    echo "Error: Failed to filter SNPs based on MGCF."; exit 1;
}

echo " Step 9: Filter SNPs based on heterozygosity"
plink2 --bfile ./plink2_output/output_filtered_mgcf --allow-no-sex --extract ./plink2_output/hz_filtered_snps.txt --make-bed --out ./plink2_output/output_filtered_hz || {
    echo "Error: Failed to filter SNPs based on heterozygosity."; exit 1;
}
echo "Creating extra files directory"
pbed_output_files="${pbed_output%.dat}"_files
mkdir $pbed_output_files
if [ "$perform_ld_pruning" = "true" ]; then
    echo "Perform LD pruning"
    plink2 --bfile ./plink2_output/output_filtered_hz --allow-no-sex --indep-pairwise "$window_size" "$step_size" "$ld_threshold" --out ./plink2_output/output_pruned || {
        echo "Error: LD pruning failed."; exit 1;
    }

    echo "Step 11: Create final dataset"
    plink2 --bfile ./plink2_output/output_filtered --allow-no-sex --extract ./plink2_output/output_pruned.prune.in --make-bed --out ./plink2_output/final || {
        echo "Error: Failed to create final dataset."; exit 1;
    }

    echo "Copy final files to output directory from pruned data"
    cp ./plink2_output/final.bim "$pbed_output_files/Composite Dataset.bim"
    cp ./plink2_output/final.bed "$pbed_output_files/Composite Dataset.bed"
    cp ./plink2_output/final.fam "$pbed_output_files/Composite Dataset.fam"
    cp ./plink2_output/final.bim "$pbed_output_files/RgeneticsData.bim"
    cp ./plink2_output/final.bed "$pbed_output_files/RgeneticsData.bed"
    cp ./plink2_output/final.fam "$pbed_output_files/RgeneticsData.fam"

else
    echo " LD pruning is not performed, use output_filtered_hz instead of final"
    cp ./plink2_output/output_filtered_hz.bim "$pbed_output_files/Composite Dataset.bim"
    cp ./plink2_output/output_filtered_hz.bed "$pbed_output_files/Composite Dataset.bed"
    cp ./plink2_output/output_filtered_hz.fam "$pbed_output_files/Composite Dataset.fam"
    cp ./plink2_output/output_filtered_hz.bim "$pbed_output_files/RgeneticsData.bim"
    cp ./plink2_output/output_filtered_hz.bed "$pbed_output_files/RgeneticsData.bed"
    cp ./plink2_output/output_filtered_hz.fam "$pbed_output_files/RgeneticsData.fam"
fi


# Cleanup
rm -rf "$TMPDIR"
echo "Pipeline completed successfully."
