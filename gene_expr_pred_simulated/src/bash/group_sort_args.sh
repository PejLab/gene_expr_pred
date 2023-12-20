
dir=$1
max_var=$2
echo $dir
for i in `seq 1 ${max_var}`; do
    awk 'NR==1; NR > 1 {print $0 | "sort -k 3"}' $dir'haplotype_logExpression_var_'$i'.txt' > $dir'haplotype_logExpression_var_'$i'_sort.txt'
done
