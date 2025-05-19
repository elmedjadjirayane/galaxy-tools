import csv
import subprocess
import sys

##########################################

### Writes minor genotype class filtered SNPs ids in a txt file

##########################################

def filter_mgcf(mgcf_threshold):
    frqx_file = "./plink2_output/output_frequencies.frqx"
    output_file = "./plink2_output/mgcf_filtered_snps.txt"

    with open(frqx_file, 'r') as infile, open(output_file, 'w') as outfile:
        reader = csv.DictReader(infile, delimiter='\t')
        writer = csv.writer(outfile)

        for row in reader:
            total = int(row["C(HOM A1)"]) + int(row["C(HET)"]) + int(row["C(HOM A2)"])
            if total > 0:  # Avoid division by zero
                mgcf = min(int(row["C(HOM A1)"]), int(row["C(HET)"]), int(row["C(HOM A2)"])) / total
                if mgcf >= float(mgcf_threshold):
                    writer.writerow([row["SNP"]])

##########################################


##########################################

### filter on HZ

##########################################

def filter_hz(max_heterozygosity):
    frqx_file = "./plink2_output/output_frequencies.frqx"
    output_file = "./plink2_output/hz_filtered_snps.txt"

    with open(frqx_file, 'r') as infile, open(output_file, 'w') as outfile:
        reader = csv.DictReader(infile, delimiter='\t')
        writer = csv.writer(outfile)

        for row in reader:
            total_genotypes = int(row["C(HOM A1)"]) + int(row["C(HET)"]) + int(row["C(HOM A2)"])
            if total_genotypes > 0:  # Avoid division by zero
                heterozygosity = int(row["C(HET)"]) / total_genotypes
                if heterozygosity <= float(max_heterozygosity):
                    writer.writerow([row["SNP"]])
##########################################

def main():
    mgcf_threshold = sys.argv[1]
    max_heterozygosity = sys.argv[2]
    filter_mgcf(mgcf_threshold)
    filter_hz(max_heterozygosity)

    
if __name__ == "__main__":
    main()
