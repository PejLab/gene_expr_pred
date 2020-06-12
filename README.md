# gene expression prediction
This script is designed to predict allele-specific expressian and total gene expression using allelic fold change (aFC). See the [manuscript](https://genome.cshlp.org/content/27/11/1872.short) for method description.

More information on estimating aFCs, could be found [here](https://github.com/wickdChromosome/aFC-n). 

## Inputs

### aFC file
In the data folder there is a file named aFC.txt. This file provids allelic fold change for each eQTL. In this file for each variant in columnn "variant_id" we have the following information:

- gene_id : gene associated to that variant (required)

- log2_aFC : The effect size of the variant (required)

- rank : rank of that variant in eQTL calling process (optional)

- CHROM : the chromosome of the variant (optional)

- POS : position of that variant (optional)

- REF : The Reference allele (optional)

- ALT : The Alternative allele (optional)

### VCF file

A tabix indexed gzip compressed VCF file containing genotypes.

Genotype, encoded as allele values separated by either of / or |. " /" means genotype unphased and "|" means genotype phased. The allele values are 0 for the reference allele (what is in the REF  field), 1 for the allele listed in ALT. For diploid calls examples could be 0/1, 1|0.

**The REF and ALT information should match the REF and ALT information in aFC.txt**

## Arguments

### Required

- **--aFC_path** : File containing allelic fold change for each eQTL. For each variant in columnn "variant_id" should contain 'gene_id', gene associated to that variant and 'log2_aFC', the effect size of the variant. 
- **--sep** : The seperator of the aFC file.
- **--vcf_path** : Tabix indexed and gzipped VCF file containing sample genotypes.
- **--variant_max** : Maximum number of variants per gene to be processed.
- **--geno** : Which field in VCF to use as the genotype. By default 'GT' = genotype
- **--output** : Output file.

## Output file

- **gene_expression.txt** : Contains predicted gene expression for each gene_id and individual
- **ASE.txt** : Contains predicted ASE for each gene_id and individual

# Resources

## gene_expr_pred.py
This script uses the lookup tables to predict expression for each haplotype, reading the individual genotypes from vcf file. To run the script use the following command:

```Shell
    python gene_expr_pred.py --aFC_path aFC_path --sep seperator --vcf_path vcf_path --variant_max max_no_of_eQTLs_per_gene --geno genotype_field_VCF  --output/--o output_file
``` 
## gene_expression_lookupTable.R

This R script, counts the number of variants for each gene and produces lookup tables representing log transformed expression values for all genotypes. The lookup tables are generated based on the number of the eQTLs for each gene, thus "haplotype_logExpression_var_2" includes genes with two eQTLS. To run the script use the following command:

```Shell
    Rscript gene_expression_lookupTable.R  aFC_path output_path sep variant_max
```    
