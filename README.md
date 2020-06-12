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

A phased vcf file is required to extract the genotypes of individuals.

Genotype, encoded as allele values separated by either of / or |. " /" means genotype unphased and "|" means genotype phased. The allele values are 0 for the reference allele (what is in the REF  field), 1 for the allele listed in ALT. For diploid calls examples could be 0/1, 1|0.

**The REF and ALT information should match the REF and ALT information in aFC.txt**

Here is an example to predict the expression of a specific gene for an individual using these inputs


# Resources

## R/gene_expression_lookupTable.R

This R script gets a sorted aFC file (**sorted based on gene_id**), counts the number of variants for each gene and produces lookup tables representing expression values for all genotypes. The lookup tables are generated based on the number of the eQTLs for each gene, thus "lookupTable_variantNo_2" includes genes with two eQTLS. To run the script use the following command:

```Shell
    Rscript gene_expression_lookupTable.R data\aFC_Whole_Blood.txt
```    
## python/ASE_prediction.ipynb
This script uses the lookup tables to predict expression for each haplotype, reading the individual genotypes from vcf file.





