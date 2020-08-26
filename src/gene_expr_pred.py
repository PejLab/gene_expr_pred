#!/usr/bin/env python
# coding: utf-8




import pandas as pd
import numpy as np
from itertools import islice
import os
import pysam
import gzip
import pdb
import time
import argparse
import warnings
import sys
from tqdm.auto import tqdm



def genotype_finder(path_vcf, path_lookup_table,max_var,output_path,vcf_geno):
#looking for genotypes from vcf file

    def genotypeOfHaplotypes(genotype_info):
        
        # set 22 for nonphased genotypes
        if not ("|" in genotype_info):
            return "22"
    
        genotype_info_Split = genotype_info.split("|")
    
        genotype_hap1 = 0 if genotype_info_Split[0] == "0" else 1
        genotype_hap2 = 0 if genotype_info_Split[1] == "0" else 1
        return (str(genotype_hap1) + str(genotype_hap2))
    
    def sample_column_map(path, start_col=9, line_key="#CHR"):
        stream_in = gzip.open(path, "r")

        out_map = {}
        for line in stream_in:
            if isinstance(line, bytes) and not isinstance(line, str):
                line = line.decode()
            if line_key in line:
                line = line.rstrip().split("\t")
                for i in range(start_col,len(line)):
                    out_map[line[i]] = i

                break

        stream_in.close()
        return(out_map)

    
    
    

    vcf_map = sample_column_map(os.path.expanduser(vcf_path))
    tabix_vcf = pysam.Tabixfile(os.path.expanduser(vcf_path),"r")

    for var_count in range(1,max_var+1):
        lookup_table= pd.read_table(path_lookup_table+'haplotype_logExpression_var_'+str(var_count)+'_sort.txt', sep=" ")
        lookup_table = lookup_table.reset_index(drop= True)
        
        prediction_matrix_haplotype1= pd.DataFrame();
        prediction_matrix_haplotype2= pd.DataFrame();
        
        hap1_individual_genotype_lst = [];
        hap2_individual_genotype_lst = [];
        temp = pd.DataFrame();
        individual_starting_column = 9
        pred_index =1 ;

        individuals = vcf_map.keys()
        individuals_lst = list(vcf_map.keys())

        genotype_info_hap1 = dict.fromkeys(np.arange(len(lookup_table)))
        genotype_info_hap2 = dict.fromkeys(np.arange(len(lookup_table)))
        
        count = 0


        for index, row in tqdm(lookup_table.iterrows(),total=lookup_table.shape[0],position=0, leave=True, desc="   "+str(var_count)+"-eQTL genes"):
            
#            tqdm.write(str(var_count)+"-eQTL genes", end = "\r")
            columns = pd.Series();
            columns = columns.append(pd.Series({'gene_id':row['gene_id']}))
            
            if ('chrX' in row['variant_id_1']):# chrX won't be found in vcf file
                for chrx_count in range(1,var_count+1):
                    var= 'variant_id_' + str(chrx_count)
                    columns = columns.append(pd.Series({var:row[var]}))
                    temp = columns.to_frame().T; 
                prediction_matrix_haplotype1= prediction_matrix_haplotype1.append(temp,ignore_index=True);
                prediction_matrix_haplotype2= prediction_matrix_haplotype2.append(temp,ignore_index=True);
                genotype_info_hap1[index] = ['2' for i in range(len(individuals_lst))]
                genotype_info_hap2[index] = ['2' for i in range(len(individuals_lst))]      
                continue;
            hap1_individual_genotype_lst.clear();
            hap2_individual_genotype_lst.clear();
            
            count = count + 1
#            print (count, end='\r')
            for i in range(1,var_count+1):
                var= 'variant_id_' + str(i)
                columns = columns.append(pd.Series({var:row[var]}))
                temp = columns.to_frame().T; 

                var_split = row[var].split('_')
                # reading vcf file

                
                records = tabix_vcf.fetch(var_split[0],int(var_split[1])-1,int(var_split[1]))
                snp_found = 0
                for record in records:
                    cols = record.rstrip().split("\t")
                    if cols[2] == row[var]:

                        gt_index = cols[8].split(":").index(vcf_geno)
                        snp_found = 1
                        break
                
                

                
                if (snp_found==0):     # chrX won't be found
                    gt_index = 1000
                    print("WARNING: eSNP %s not found in VCF"%(row[var]))
                
                hap_ind = 0;

    
                for ind in individuals:
            
                    if (gt_index ==1000):
                        genotype_both = './.'
                    else:  
                        sample_col = cols[vcf_map[ind]]
                        genotype_both = sample_col.split(":")[gt_index]


                    genotype_both_haps= genotypeOfHaplotypes(genotype_both);
            
                    if (i>1):
                       
                        hap1_individual_genotype_lst.insert(hap_ind, hap1_individual_genotype_lst[hap_ind]+str(genotype_both_haps[0]) + ", ");
                        del hap1_individual_genotype_lst[hap_ind + 1];
                        hap2_individual_genotype_lst.insert(hap_ind, hap2_individual_genotype_lst[hap_ind]+str(genotype_both_haps[1]) + ", ");
                        del hap2_individual_genotype_lst[hap_ind + 1];
                    else:
                        hap1_individual_genotype_lst.insert(hap_ind, str(genotype_both_haps[0]) + ", ");
                        hap2_individual_genotype_lst.insert(hap_ind, str(genotype_both_haps[1]) + ", ");
                    hap_ind = hap_ind+1;
                    
                    
            prediction_matrix_haplotype1= prediction_matrix_haplotype1.append(temp,ignore_index=True);
            prediction_matrix_haplotype2= prediction_matrix_haplotype2.append(temp,ignore_index=True);
    
    
    
    
            #Create genotype array to be inserted for each gene        
            clean_genotype_lst1 = [genotype[:-2] for genotype in hap1_individual_genotype_lst]
            clean_genotype_lst2 = [genotype[:-2] for genotype in hap2_individual_genotype_lst]
            
                
            #Do the insertion one gene at a time
            genotype_info_hap1[index] = clean_genotype_lst1
            genotype_info_hap2[index] = clean_genotype_lst2                       
        
        
        output_genotype_hap1 = pd.concat([prediction_matrix_haplotype1, pd.DataFrame.from_dict(genotype_info_hap1, columns=individuals, orient='index')],axis = 1, sort = False);
        output_genotype_hap2 = pd.concat([prediction_matrix_haplotype2, pd.DataFrame.from_dict(genotype_info_hap2, columns=individuals, orient='index')],axis = 1, sort = False);
        
        
        

        output_genotype_hap1.to_csv(os.path.expanduser(output_path+'/genotype_hap1_var_count_'+str(var_count)+".csv"));
        output_genotype_hap2.to_csv(os.path.expanduser(output_path+'/genotype_hap2_var_count_'+str(var_count)+".csv"));    



def expr_pred(lookup_Table_path, max_var, path_genotype, output_path):
# predicting ASE and total expression     


    def Allelic_imbalance_predicted(hap1_count,hap2_count):

        return hap1_count - hap2_count;
    
    def ASE_predicted_value(hap1_count,hap2_count):
        h1 = 2**hap1_count
        htotal = 2**hap1_count + 2**hap2_count
        p = h1 / htotal
        h1_percent = round(p,2)
        h2_percent  = round(1 - h1_percent,2)
        s = str(h1_percent)+"|"+str(h2_percent)
        
        return s

    def total_predicted_value(hap1_count,hap2_count):
        hap1_count = 2**hap1_count;        
        hap2_count = 2**hap2_count;           
    

        return np.log2(hap1_count + hap2_count);
    
    

    for var_count in range(1,max_var+1):
        lookup_Table = pd.read_table(lookup_Table_path+"haplotype_logExpression_var_"+str(var_count)+".txt", sep=" ");
        lookup_Table.columns = lookup_Table.columns.str.replace(',','').str.replace(" ", "");


        genotype_haplotype1=pd.read_csv(path_genotype+"/genotype_hap1_var_count_"+str(var_count)+".csv")
        genotype_haplotype2= pd.read_csv(path_genotype+"/genotype_hap2_var_count_"+str(var_count)+".csv")


        start_index_individuals = var_count + 2;

        

        individuals = genotype_haplotype1.columns.tolist()[start_index_individuals :]

        # make empty dict
        total_exp_output = dict.fromkeys(np.arange(len(genotype_haplotype1)))     
        ASE_exp_output =  dict.fromkeys(np.arange(len(genotype_haplotype1))) 

        # make gene_ids seperately
        genes = pd.DataFrame(np.arange(len(genotype_haplotype1)))
        count = 0
        
        total_count_lst = [];
        ASE_count_lst = [];
        
        for index, row in tqdm(genotype_haplotype1.iterrows(),total=genotype_haplotype1.shape[0],position=0, leave=True, desc="   "+str(var_count)+"-eQTL genes"):
            
            total_count_lst.clear();
            ASE_count_lst.clear();
            count = count +1
            

#            if ('chrX' in row['variant_id_1']):
#                continue;
#            print (var_count, "   ", count, end='\r')
            
            placeholder =  list(np.where(lookup_Table['gene_id'] == row['gene_id']))[0] + 1;  # index starts from 1 in lookup table 
            lookup_Table_index = placeholder [0];   



            genes.at[index,'gene_id'] = row['gene_id']


            ind_index = 0
            for ind in individuals:


                if (not '2' in str(row[ind]) and not '2' in str(genotype_haplotype2.loc[index,ind])):
               

                    x1 = str(genotype_haplotype1.loc[index,ind]).replace(",","").replace(" ", "")
                    x2 = str(genotype_haplotype2.loc[index,ind]).replace(",","").replace(" ", "")
                    
                    if (var_count == 1):
                        x1 = x1.split('.')[0]
                        x2 = x2.split('.')[0]
                    # read predicted value from lookup table
                    total_prediction = total_predicted_value(lookup_Table.loc[lookup_Table_index,x1],lookup_Table.loc[lookup_Table_index,x2])
                    total_count_lst.insert(ind_index,round(total_prediction,4))
                    
                    ASE_prediction = ASE_predicted_value(lookup_Table.loc[lookup_Table_index,x1],lookup_Table.loc[lookup_Table_index,x2])
                    ASE_count_lst.insert(ind_index,ASE_prediction)
                    
                else:
                    total_count_lst.insert(ind_index,np.nan)
                    ASE_count_lst.insert(ind_index,np.nan)
                ind_index = ind_index + 1
                
    

            total_count_lst1 = [expr for expr in total_count_lst]
            ASE_count_lst1 = [expr for expr in ASE_count_lst]
            total_exp_output[index] = total_count_lst1
            ASE_exp_output[index] = ASE_count_lst1
            

        
        # make different dicts for different number of eQTLs
        globals()['total_exp_output_df_%s'%var_count] =  pd.concat([genes,pd.DataFrame.from_dict(total_exp_output, columns=individuals, orient='index')],axis = 1, sort = False)
        globals()['ASE_exp_output_df_%s'%var_count] = pd.concat([genes,pd.DataFrame.from_dict(ASE_exp_output, columns=individuals, orient='index')], axis =1 , sort = False)
        


    frames_total = []
    frames_ASE = []
    for j in range(1,max_var+1):
        frames_total.append(globals()['total_exp_output_df_%s'%j])
        frames_ASE.append(globals()['ASE_exp_output_df_%s'%j])
        
    total_exp_output_df = pd.concat(frames_total, sort=False)
    ASE_exp_output_df = pd.concat(frames_ASE, sort = False)
    
    total_exp_output_df = total_exp_output_df.drop(total_exp_output_df.columns[0], axis = 1)
    ASE_exp_output_df = ASE_exp_output_df.drop(ASE_exp_output_df.columns[0], axis = 1)
    total_exp_output_df.to_csv(os.path.expanduser(output_path+'/gene_expression.txt'),index=False, sep='\t')
    ASE_exp_output_df.to_csv(os.path.expanduser(output_path+'/ASE.txt'),index = False, sep ='\t')





if __name__ == "__main__":
    
    
    parser = argparse.ArgumentParser()
    # REQUIRED
    parser.add_argument("--aFC_path", required=True, help="aFC path")
    parser.add_argument("--vcf_path", required=True, help="vcf path")
    parser.add_argument("--sep", required=True, help="aFC seperator")
    parser.add_argument("--variant_max",type = int, required=True, help="max no. of variants to process")
    parser.add_argument("--geno", required=False, default="GT", help="Which field in VCF to use as the genotype. By default 'GT'") 
   
    parser.add_argument("--output_path", "--o", required=True, help="Output file")
    
    # OPTIONAL
    parser.add_argument("--chr", type=str, help="Limit to a specific chromosome.")
    
    
    
    args = parser.parse_args() 
    
# sample input    
    
#    aFC_path = '~/genome_biology/data/aFC/aFC_v1.csv'
#    vcf_path = '~/data/GTEx_Analysis_v8_phASER/phASER_GTEx_v8_merged.vcf.gz'
#    output_path = '~/data/final_output/final_test/new/'
#    variant_max = 16
#    sep = ','    
#    geno = 'GT'

    
    aFC_path = args.aFC_path
    vcf_path = args.vcf_path
    output_path = args.output_path
    variant_max = int(args.variant_max)
    geno = args.geno
    sep = args.sep
    
    
    warnings.filterwarnings("ignore")
    script_path = sys.path[0]
    print(script_path)
    
    print("Run settings ...")
    print("     aFC File: %s"%(aFC_path))
    print("     Genotype VCF: %s"%(vcf_path))
    if args.chr != None:
        
        print("     Chromosome: %s"%(args.chr))

    
    
    
    
    # make the output directory
    bashCommand = 'mkdir -p ' + output_path
    os.system(bashCommand)
    
    bashCommand = 'mkdir ' + output_path + "/temp"
    os.system(bashCommand)
    
    
    
#    logging.basicConfig(filename=output_path + "/temp/"+tissue_name+ 'logs/info.log',level=logging.INFO)
#    logging.basicConfig(filename=output_path + "/temp/"+tissue_name+'logs/info.debug.log',level=logging.DEBUG)
        
#    logging.info("Run starting at " + str(datetime.datetime.now))
    
    
    # make lookup tables for all genotypes
    bashCommand = 'mkdir ' + output_path + "/temp/lookup_table/"
    os.system(bashCommand)
    
    path_lookup = output_path + "/temp/lookup_table/"
    
    print("Preparing lookup tables ...")
    if (args.chr == None):
        bashCommand='Rscript '+script_path+'/gene_expression_lookupTable.R '+ aFC_path + " "+path_lookup+" "+sep+" "+str(variant_max)
    else:
        bashCommand='Rscript '+script_path+'/gene_expression_lookupTable.R '+ aFC_path + " "+path_lookup+" "+sep+" "+str(variant_max)+" "+args.chr
    os.system(bashCommand)
    
    bashCommand= 'bash '+script_path+'/bash/group_sort_args.sh '+path_lookup+" "+str(variant_max)
    os.system(bashCommand)
    
    
    
    
    
    # read genotypes
    print("Reading genotypes ...")
    bashCommand = 'mkdir ' + output_path + "/temp/genotypes/"
    os.system(bashCommand)
    genotype_finder(vcf_path,path_lookup,variant_max,output_path+"/temp/genotypes/", geno)
    
    
    
    
    # predicting total and ASE expressions
    print('Expression prediction ...')
    
    
    genotype_path = output_path + "/temp/genotypes/"
    
    bashCommand = 'mkdir ' + output_path + "/expressions/"
    os.system(bashCommand)
    
    
    expr_pred(path_lookup, variant_max, genotype_path, output_path+"/expressions/")
    

