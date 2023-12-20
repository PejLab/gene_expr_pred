#lookup table for n variants

library(stringr)
library(data.table)

#Initialize
args <- commandArgs(trailingOnly = TRUE)
filename = args[1]
output = args[2]
seperator = args[3]
max_var = args[4]
chrom = args[5]
mean_val = args[6]
sd_val = args[7]


# The AFC file should be sorted on gene_id
AFC_data=read.table(file=filename, header=TRUE, sep=seperator)
AFC <- AFC_data[order(AFC_data$gene_id),]

# if the chr argument is not null
if (!is.na(chrom)){
    chr_id = paste0(chrom,'_')
    AFC<-AFC[AFC$variant_id %like% chr_id, ]
    if (dim(AFC)[1] == 0){
        stop("No gene for the specified chromosome")
    }
    }



# maximum number of variants
maximum_variants <- as.integer(max_var)
mean_val <- as.integer(mean_val)
sd_val <- as.integer(sd_val)

# h1: first haplotype
# h2: second haplotype
# s: effect size of variants
gene_expression_estimation<-function(h1,h2,s){
  log_expression_h1<-t(h1)%*%s
  log_expression_h2<-t(h2)%*%s
  total_expression<-2^(log_expression_h1)+2^(log_expression_h2)
  expression_ratio_h1<-2^(log_expression_h1)/total_expression
  expression_ratio_h2<-2^(log_expression_h2)/total_expression
  result<-c(log_expression_h1,log_expression_h2,total_expression,expression_ratio_h1,expression_ratio_h2)
  return(result)
  
}

# producing all genotypes
number2binary = function(number, noBits) {
  binary_vector = rev(as.numeric(intToBits(number)))
  if(missing(noBits)) {
    return(binary_vector)
  } else {
    binary_vector[-(1:(length(binary_vector) - noBits))]
  }
}

########################

afc_table<-data.frame()
afc_table_index <-1
for (var_count in 1:maximum_variants){    
  lookup_table<-data.frame()
  lookup_table_index<-1
  count<-1
  result<-c()
  count_result<-c()
  for (index in 1:nrow(AFC)){
    
    # check the number of variants 
    current_gene<-paste(AFC[index,'gene_id'])
    next_gene<-paste(AFC[index+1,'gene_id'])
    
    if((current_gene == next_gene)){
      count<-count+1
      
    }
    else{
      
      if (var_count==count){
        
        gene_id<-gsub("\\..*","",AFC[index,'gene_id'])   #remove version No.
        lookup_table[lookup_table_index,'gene_id']<-paste0(gene_id)
        AFC_vector<-c()
        var_index<-1
        for (temp_index in (index-count+1):index){
                
          AFC_value<-rnorm(1, mean = mean_val, sd = sd_val) #simulate AFC value
          
          # save zero for nan values
          if (AFC_value=='NaN'| AFC_value=='nan' | is.na(AFC_value))
            AFC_value<-0
          
          
          # Cut off afc values to remove outliers
          if (AFC_value < -6.64){
            AFC_value <- -log2(100)
          }
          if (AFC_value > 6.64){
            AFC_value <- log2(100)
          }
          # Cut off afc values
          afc_table[afc_table_index,'gene_id'] = gene_id
          afc_table[afc_table_index,'variant_id'] = AFC_data[temp_index,'variant_id']
          afc_table[afc_table_index,'afc_simulated'] = AFC_value
          afc_table_index <- afc_table_index +1 
          
            AFC_vector<- c(AFC_vector,AFC_value)
          
          # create lookup table
          
          lookup_table[lookup_table_index,paste0('variant_id_',var_index)]<-paste0(AFC[temp_index,'variant_id'])
          
#          lookup_table[lookup_table_index,paste0('rank_',var_index)]<-paste0(AFC[temp_index,'rank'])
          var_index<-var_index+1    
        }
        for (j in 0:(2^count-1)){
          
          ref_alt_vector<-c(number2binary(j,count))    #get all possible genotypes
          ref_vector<-c(replicate(count,0))
          output<-gene_expression_estimation(ref_alt_vector,ref_vector,AFC_vector)
          lookup_table[lookup_table_index,paste0("",toString(ref_alt_vector))]<-paste0(output[1])
          
          
          
        }
        
        lookup_table_index<-lookup_table_index+1   
        
        
        
      }
      count<-1
    }#end of else
  }
  # create output  
#  print(args[2])
#  write.csv(lookup_table,file=paste0(args[2],"/haplotype_logExpression_var_",var_count,".csv"))
  write.table(lookup_table,file=paste0(args[2],"/haplotype_logExpression_var_",var_count,".txt"))
  write.table(afc_table,file=paste0(args[2],"/simulated_afc.txt"))

  
  
}
