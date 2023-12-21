# Allele-specific and genetically-regulated gene expression prediction
We simulated log2 effect sizes for each eQTL from standard Normal distribution. Expected gene expression is estimated as the sum of the two haplotypic counts using a simulated e0 (expression for haplotype carrying reference allele), and the simulated sequencing read count for each gene is generated using a Poisson distribution. To generate simulated ASE data, we calculated the expected reference expression ratio for each individual, and simulated discrete read counts for each haplotype using Binomial distribution.

## Inputs

### input file
In the data folder there is a file named input_sample.txt. In this file for each variant in columnn "variant_id" we have the following information:

- gene_id : gene associated to that variant (required)


### VCF file

A tabix indexed gzip compressed VCF file containing genotypes.

Genotype, encoded as allele values separated by either of / or |. " /" means genotype unphased and "|" means genotype phased. The allele values are 0 for the reference allele (what is in the REF  field), 1 for the allele listed in ALT. For diploid calls examples could be 0/1, 1|0.

The assumption is that the build version is aligned to the GRCh38 assembly (The format is : chr_position_ref_alt_b38).

**The REF and ALT information should match the REF and ALT information in aFC file**

## Arguments

### Required

- **--input_path** : File containing gene_id and variant_id
- **--sep** : The seperator of the aFC file.
- **--vcf_path** : Tabix indexed and gzipped VCF file containing sample genotypes.
- **--variant_max** : Maximum number of variants per gene to be processed.
- **--geno** : Which field in VCF to use as the genotype. By default 'GT' = genotype
- **--output** : Output directory.
- **--phased** : If True only phased genotypes in VCF will be considered. If False both phased and unphased genotypes will be considered, in case of unphased genotypes the results are approximate. The default is True.
- **--mean** : The mean for normal distribution for simulating aFcs. By default 'mean' = 0
- **--sd** : The standard deviation for normal distribution for simulating aFcs. By default 'sd' = 1

### Optional
- **--chr** : Limit to a specific chromosome.

## Output file

- **gene_expression.txt** : Contains predicted gene expression for each gene_id and individual. The expression is relative to the expression of haplotype carrying reference allele(s) (e<sub>R</sub>) as described in the [manuscript](https://www.biorxiv.org/content/10.1101/2022.01.28.478116v1).
- **ASE.txt** : Contains relative estimated ASE for each gene_id and individual.
- **simulated_afc.txt** : Simulated aFC

# Resources

## gene_expr_pred.py
This script uses the lookup tables to predict expression for each haplotype, reading the individual genotypes from vcf file. To run the script use the following command:

```Shell
    python gene_expr_pred.py --input_path input_data --sep seperator --vcf_path vcf_path --variant_max max_number_variants --geno vcf_genotype_index --output/--o output_dir --mean mean --sd sd
``` 

## License
This project is covered under the Apache 2.0 License.
