#!/bin/bash


USAGE='
usage: report_stats [options]

-s		sites file of m6A or m6Am to be annotated, required
-m		methlylation type: m6A or m6Am, required
-g		genome used, hg38 (default) or mm10  
-l		length of sequence around sites to be extracted, default is 5
'

while getopts 'hs:m:g:l:' opt
do
  case $opt in
    h) echo "$USAGE"
       exit 1;;
    s) sites=$OPTARG;;
    m) methyl=$OPTARG;;
    g) species=$OPTARG;;
    l) motif_len=$OPTARG;;
    ?) echo "Unrecognized options $1!"; exit 1;;
  esac
done

species=${species:-'hg38'}
motif_len=${motif_len:-'5'}

if [[ ${species} = "hg38" ]]; then
    genome2="/lustre2/chengqiyi_pkuhpc/zhouzhe/genome/GLORI/hg38/genome/hg38_chr_only.fa"
elif [[ ${species} = "mm10" ]]; then
    genome2="/lustre2/chengqiyi_pkuhpc/zhouzhe/genome/GLORI/mm10/genome/GCF_000001635.26_GRCm38.p6.only.fa"
else
    echo "Please specifiy genome used (hg38 or mm10)!"
    exit 1
fi
metaPlotR_dir=/lustre2/chengqiyi_pkuhpc/zhouzhe/software/metaPlotR/
annot_prx=/lustre2/chengqiyi_pkuhpc/zhouzhe/genome/anno_metaPlotR/${species}/${species}

# 1. MetaPlotR
echo -e "\n[ $(date) ]: MetaPlotR --------"
source activate metaPlotR
if [[ ${methyl} == "m6A" ]]; then
  cat ${sites} | awk -v OFS='\t' 'NR>1 {print $1, $2-1, $2, $4, 0, $3}' > ${sites}.bed
else
  cat ${sites} | awk -v OFS='\t' 'NR>1 {print $1, $2-1, $2, $5, 0, $3}' > ${sites}.bed
fi
perl ${metaPlotR_dir}/annotate_bed_file.pl --bed ${sites}.bed \
    --bed2 ${annot_prx}_annot.sorted.bed_longest > ${sites}_annot.bed
perl ${metaPlotR_dir}/rel_and_abs_dist_calc.pl --bed ${sites}_annot.bed \
    --regions ${annot_prx}_region_sizes.txt_longest > ${sites}_dist.measures.txt
rm -rf ${sites}.bed ${sites}_annot.bed
source activate zzhou_bio2

# 02. get bases around m6A sites
echo -e "\n[ $(date) ]: Get ${motif_len} bases around sites --------"
if [[ ${methyl} == "m6A" ]]; then
  cat ${sites} | awk -v OFS='\t' -v l="$(( motif_len / 2 ))" 'NR!=1 && $2>=3 {print $1, $2-l-1, $2+l, $4, 0, $3}' \
    > ${sites}_len${motif_len}.bed
else
  cat ${sites} | awk -v OFS='\t' -v l="$(( motif_len / 2 ))" 'NR!=1 && $2>=3 {print $1, $2-l-1, $2+l, $5, 0, $3}' \
    > ${sites}_len${motif_len}.bed
fi
bedtools getfasta -s -fi ${genome2} -bed ${sites}_len${motif_len}.bed -bedOut > ${sites}_len${motif_len}_motif.bed
rm -rf ${sites}_len${motif_len}.bed
