# Allele-specific and genetically-regulated gene expression prediction
This script is designed to predict allele-specific expressian and total gene expression using allelic fold change (aFC). See the [manuscript](https://genome.cshlp.org/content/27/11/1872.short) for method description.

More information on estimating aFCs, could be found [here](https://github.com/wickdChromosome/aFC-n). 

## Inputs

### aFC file
In the data folder there is a file named aFC.txt. This file provids allelic fold change for each eQTL. In this file for each variant in columnn "variant_id" we have the following information:

- gene_id : gene associated to that variant (required)

- log2_aFC : The effect size of the variant in log<sub>2</sub> (required)

### VCF file

A tabix indexed gzip compressed VCF file containing genotypes.

Genotype, encoded as allele values separated by either of / or |. " /" means genotype unphased and "|" means genotype phased. The allele values are 0 for the reference allele (what is in the REF  field), 1 for the allele listed in ALT. For diploid calls examples could be 0/1, 1|0.

**The REF and ALT information should match the REF and ALT information in aFC file**

## Arguments

### Required

- **--aFC_path** : File containing allelic fold change for each eQTL. For each variant in columnn "variant_id" should contain 'gene_id', gene associated to that variant and 'log2_aFC', the effect size of the variant. 
- **--sep** : The seperator of the aFC file.
- **--vcf_path** : Tabix indexed and gzipped VCF file containing sample genotypes.
- **--variant_max** : Maximum number of variants per gene to be processed.
- **--geno** : Which field in VCF to use as the genotype. By default 'GT' = genotype
- **--output** : Output directory.
- **--phased** : If True only phased genotypes in VCF will be considered. If False both phased and unphased genotypes will be considered, in case of unphased genotypes the results are approximate. The default is True. 

### Optional
- **--chr** : Limit to a specific chromosome.

## Output file

- **gene_expression.txt** : Contains predicted gene expression for each gene_id and individual. The expression is relative to the expression of haplotype carrying reference allele(s) (e<sub>R</sub>) as described in the [manuscript](https://www.biorxiv.org/content/10.1101/2022.01.28.478116v1).
- **ASE.txt** : Contains relative estimated ASE for each gene_id and individual.

# Resources

## gene_expr_pred.py
This script uses the lookup tables to predict expression for each haplotype, reading the individual genotypes from vcf file. To run the script use the following command:

```Shell
    python gene_expr_pred.py --aFC_path aFC_path --sep seperator --vcf_path vcf_path --variant_max max_number_variants --geno vcf_genotype_index --output/--o output_dir
``` 
## gene_expression_lookupTable.R

This R script, counts the number of variants for each gene and produces lookup tables representing log transformed expression values for all genotypes. The lookup tables are generated based on the number of the eQTLs for each gene, thus "haplotype_logExpression_var_2" includes genes with two eQTLS. To run the script use the following command:

```Shell
    Rscript gene_expression_lookupTable.R  aFC_path output_path sep variant_max
```  

## Data
The folder contains a sample of aFC file and a tabix indexed VCF file. The VCF file is a subset of original file generated by permutation of genotypes. (It should take less than a minute to generate output on this sample size data).

## License
This project is covered under the Apache 2.0 License.
