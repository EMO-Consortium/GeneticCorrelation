
## reformat data as numeric
```R
.libPaths("/groups/umcg-griac/tmp01/projects/umcg-cqi/software/Rpackage/4.0/")
library(tidyverse)
aa<-read_table("without_UKB_TC_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results")
bb<-aa[,c("rsID","A2","A1","BETA","SE","pvalue")]
write.table(bb,file="TC_INV_EUR_HRC_1KGP3.txt",sep="\t",row.names=F,quote=F)

```

## summary statisitcs
```bash
cd /groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Tools/ldsc/
source activate /groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Tools/Conda/envs/ldsc
ml PythonPlus/2.7.16-foss-2018b-v19.08.1

cd /groups/umcg-griac/tmp01/projects/umcg-cqi/tmp.Yanni

/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Tools/ldsc/munge_sumstats.py \
--sumstats TC_INV_EUR_HRC_1KGP3.txt \
--N 500000 \
--out TC_INV_EUR_HRC_1KGP3 --merge-alleles w_hm3.snplist \
--chunksize 500000

```
## run ldsc
```bash
/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Tools/ldsc/ldsc.py \
--ref-ld-chr /groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Reference/eur_w_ld_chr/ \
--out /groups/umcg-griac/tmp01/projects/umcg-cqi/tmp.Yanni/out/FLI_vs_TC_INV_EUR_HRC_1KGP3 \
--rg /groups/umcg-griac/tmp01/projects/umcg-cqi/tmp.Yanni/FLI.sumstats.gz,/groups/umcg-griac/tmp01/projects/umcg-cqi/tmp.Yanni/TC_INV_EUR_HRC_1KGP3.sumstats.gz \
--w-ld-chr /groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Reference/eur_w_ld_chr/ 

```
