# GeneticCorrelation

This file is about the GWAS pre-process at Linux system.

For genetic correlation files, see sumstats.

For Oxidative Stress in CD patients, see Oxidative Stree Project.

Contact person: hushx9@mail.sysu.edu.cn


***Summary statistics***

```

cd /groups/umcg-weersma/tmp01/Shixian/GeneticCorrelation/MAFLD/Tools/ldsc/

source /groups/umcg-weersma/tmp01/Shixian/Tools/Conda/bin/activate ldsc

cd /groups/umcg-gastrocol/tmp01/Shixian/GeneticCorrelation/Mental_disease/filter/

/groups/umcg-weersma/tmp01/Shixian/GeneticCorrelation/MAFLD/Tools/ldsc/munge_sumstats.py \
--sumstats XXXX.sumstats.txt \
--N 382549 \
--out XXXX --merge-alleles w_hm3.snplist \
--chunksize 500000

```
Check allele direction.


***ldsc***

```
cd /groups/umcg-weersma/tmp01/Shixian/GeneticCorrelation/MAFLD/Tools/ldsc/

source /groups/umcg-weersma/tmp01/Shixian/Tools/Conda/bin/activate ldsc

for i in /groups/umcg-gastrocol/tmp01/Shixian/GeneticCorrelation/Mental_disease/sumstats/*gz
do

p1=$(basename $i)

for n in /groups/umcg-gastrocol/tmp01/Shixian/GeneticCorrelation/Mental_disease/sumstats/*gz
do

p2=$(basename $n)

./ldsc.py \
--ref-ld-chr /groups/umcg-gastrocol/tmp01/Shixian/GeneticCorrelation/Mental_disease/filter/eur_w_ld_chr/ \
--out /groups/umcg-gastrocol/tmp01/Shixian/GeneticCorrelation/Mental_disease/LDSC_out/$p1\_vs_$p2 \
--rg $i,$n \
--w-ld-chr /groups/umcg-gastrocol/tmp01/Shixian/GeneticCorrelation/Mental_disease/filter/eur_w_ld_chr/ 

done

done

```

***Merge results from ldsc***

```
for i in *sumstats.gz.log

do

intercept=$(grep Intercept $i | sed -n "3,3p"  |awk '{print $2}')
correlation=$(grep "Genetic Correlation" $i | sed -n "2,2p" | awk '{print $3}')
Pvalue=$(grep "P:" $i | awk '{print $2}')

echo -e "$i \t $intercept \t $correlation \t $Pvalue" >> Merged.results.txt

echo $i

done

```

***Merge LAVA results and annotate with genes***
```

library(data.table)

filenames <- list.files("LAVAoutput/", pattern="bivar.lava", full.names=TRUE)
lava.out <- rbindlist(lapply(filenames,fread))
lava.out=na.omit(lava.out)

# select pairs, between gastro diseases and respiratory diseases
disease_gastro=c("CD","UC","Celiac_disease","Diverticular_disease","Esophageal_cancer","Gastro_reflux","Gastric_cancer",
                 "IBS","Peptic_ulcer","Chronic_gastritis","CRC")
disease_respiratory=c("Allergy","MtoS_asthma","FEV1","FVC","Asthma","FEV1_to_FVC","Atopic_dermatitis","hayfever",
                      "AOA","COA")

tmp1=c()
for(i in disease_gastro){
  for(n in disease_respiratory){
    mm=paste(i,n,sep = "_")
    tmp1=append(tmp1,mm)
  }
}
tmp2=c()
for(i in disease_respiratory){
  for(n in disease_gastro){
    mm=paste(i,n,sep = "_")
    tmp2=append(tmp2,mm)
  }
}
uniq_pairs=unique(c(tmp1,tmp2))

lava.out$pairs=paste(lava.out$phen1,lava.out$phen2,sep = "_")
unique(lava.out$pairs[!lava.out$pairs %in% uniq_pairs])
lava.out=lava.out[lava.out$pairs %in% uniq_pairs,]
lava.out$FDR=p.adjust(lava.out$p,method = "BH")

lava.sig=lava.out[lava.out$FDR<0.05,]

library(biomaRt)
mart <- useMart("ensembl")
mart <- useMart( "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "www.ensembl.org" ) # could be disrupted, try more times

attributes <- c("ensembl_gene_id","start_position","end_position","strand","hgnc_symbol","chromosome_name","ucsc","band")
filters <- c("chromosome_name","start","end")

locus_summary=lava.sig[,c("locus","chr","start","stop")]
locus_summary=locus_summary[!duplicated(locus_summary),]
locus_summary$symbol=NA
locus_summary$ensemble=NA
for(i in 1:nrow(locus_summary)){
  tmp.values=list(chromosome=as.character(locus_summary$chr[i]),start=as.character(locus_summary$start[i]),end=as.character(locus_summary$stop[i]))
  tmp.genes=getBM(attributes=attributes, filters=filters, values=tmp.values, mart=mart)
  tmp.symbol=unique(tmp.genes$hgnc_symbol)
  tmp.ensemble=unique(tmp.genes$ensembl_gene_id)
  
  tmp.symbol=paste(tmp.symbol, collapse = ',')
  tmp.ensemble=paste(tmp.ensemble, collapse = ',')
  locus_summary$symbol[i]=(tmp.symbol)
  locus_summary$ensemble[i]=tmp.ensemble
  
  cat(i,"\n")
}
lava=merge(lava,locus_summary[,c("locus","symbol","ensemble")],by="locus",all=F)

```


***Split locus: to extract LAVA enriched SNPs or genes for next coloc or SMR***

```
cat Lava.significant.results.txt | while read line

do

start=$(echo $line | cut -f5 -d" ")
stop=$(echo $line | cut -f6 -d" ")
chr=$(echo $line | cut -f4 -d" ")
number=$(echo $line | cut -f1 -d" ")

zcat ../raw_summary/CD.txt.gz | awk -v STA="$start" -v STO="$stop" -v CHR="$chr" '{if ($3==CHR && $4>STA && $4<STO) print $2}' > Locus.$chr.$number.txt

echo -e "CHR $chr Locus $number done"

done
```
***Extract eQTL from GTEx***

```
tar -zxvf GTEx_Analysis_v7_eQTL_all_associations.tar.gz GTEx_Analysis_v7_eQTL_all_associations/Small_Intestine_Terminal_Ileum.allpairs.txt.gz
tar -zxvf GTEx_Analysis_v7_eQTL_all_associations.tar.gz GTEx_Analysis_v7_eQTL_all_associations/Colon_Sigmoid.allpairs.txt.gz
tar -zxvf GTEx_Analysis_v7_eQTL_all_associations.tar.gz GTEx_Analysis_v7_eQTL_all_associations/Colon_Transverse.allpairs.txt.gz


```

***Subset GTEx SMR summary***
```
ml PLINK/1.9-beta6-20190617
/groups/umcg-gastrocol/tmp01/Shixian/Tools/smr_v1.3.1_linux_x86_64_static/smr_v1.3.1_linux_x86_64_static \
--beqtl-summary /groups/umcg-gastrocol/tmp01/Shixian/DUOX2/eQTLTissue/Colon_Sigmoid \
--extract-probe All.LAVA.gene.txt \
--make-besd \
--out  Colon_Sigmoid.LAVA
```

***re-format mbQTL summaries***
```
file="GWAS_HSERMETANA.PWY..L.methionine.biosynthesis.III.rearrange.tsv.gz"
name=${file%.tsv*}

zcat $file | awk 'OFS="\t"{if($1!="NA")print $3,$1,$4,7899,$5,$6,$7,$8,$9,5,$2}' | awk '!seen[$2]++'> $name.tsv
sed -i "1s/.*/chr rs  ps  n_miss  allel1  allel0  af  beta  se  l_remle p_wald/" $name.tsv
cat Flist.file.all.txt | { head -1; grep $name; } > $name.flist

ml PLINK/1.9-beta6-20190617

/groups/umcg-gastrocol/tmp01/Shixian/Tools/smr_v1.3.1_linux_x86_64_static/smr_v1.3.1_linux_x86_64_static \
--eqtl-flist $name.flist \
--cis-wind 1 --trans-wind 1000 --peqtl-trans 1.0e-5 --peqtl-other 1.0e-5 \
--gemma-format --make-besd \
--out  $name.BESD

/groups/umcg-gastrocol/tmp01/Shixian/Tools/smr_v1.3.1_linux_x86_64_static/smr_v1.3.1_linux_x86_64_static \
--bfile /groups/umcg-gastrocol/tmp01/Shixian/GeneticCorrelation/Metabolic_disease/reference/1000G.Euro \
--gwas-summary /groups/umcg-gastrocol/tmp01/Shixian/GeneticCorrelation/Metabolic_disease/GWAS_beta/LDL.txt \
--beqtl-summary $name.BESD \
--peqtl-smr 1e-05 --thread-num 10 \
--diff-freq-prop 0.9 --diff-freq 1 --out SMR.$name \
--trans --trans-wind 1000 

```


***Coloc eQTL***

```
for i in Locus*txt

do

locus=${i%.txt}
awk ' FNR==NR { a[$1]=$1; next } $1 in a { print $2 }' $i ../All.pairs.txt > $locus.tmp
cat $locus.tmp | sort | uniq > $locus.gene
rm $locus.tmp

done

```



```
cd  Locus.11.1733
cat Locus.11.1733.gene | while read line

do
awk -v gene="$line" '{if($2==gene) print}' ../All.pairs.txt > $line.summary.txt
echo $line
done

```

extract eQTLGen
```
cd  Locus.11.1733

/groups/umcg-gastrocol/tmp01/Shixian/Tools/smr_v1.3.1_linux_x86_64_static/smr_v1.3.1_linux_x86_64_static \
--beqtl-summary /groups/umcg-gastrocol/tmp01/Shixian/GeneticCorrelation/Mental_disease/eQTL/eQTLGen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense \
--extract-snp Locus.11.1733.txt \
--query 1 \
--out Locus.11.1733.eQTLGen 

awk '{print $7}' Locus.11.1733.eQTLGen.txt  | sort | uniq > Locus.11.1733.eQTLGen.gene.txt

cat Locus.11.1733.eQTLGen.gene.txt | while read line

do
awk -v gene="$line" 'OFS="\t"{if($7==gene)print $1,$12,$13}' Locus.11.1733.eQTLGen.txt > Locus.11.1733.$line.trait
sed -i "1i\SNP\tbeta\tse" Locus.11.1733.$line.trait
echo $line

done


```

extract GTEx

```

for n in /groups/umcg-gastrocol/tmp01/Shixian/GeneticCorrelation/Mental_disease/eQTL/Coloc/Locus*

do

i=$(basename $n)
cd $n

/groups/umcg-gastrocol/tmp01/Shixian/Tools/smr_v1.3.1_linux_x86_64_static/smr_v1.3.1_linux_x86_64_static \
--beqtl-summary //groups/umcg-gastrocol/tmp01/Shixian/GeneticCorrelation/Mental_disease/eQTL/GTEx/Small_Intestine_Terminal_Ileum/Small_Intestine_Terminal_Ileum \
--extract-snp //groups/umcg-gastrocol/tmp01/Shixian/GeneticCorrelation/Mental_disease/LAVA_locus/$i.txt \
--query 1 \
--out $i.GTEx.Ileum

awk '{print $7}' $i.GTEx.Ileum.txt | sort | uniq > $i.Ileum.gene.txt

cat $i.Ileum.gene.txt | while read line

do
awk -v gene="$line" 'OFS="\t"{if($7==gene)print $1,$12,$13}' $i.GTEx.Ileum.txt > $i.$line.Ileum.trait

sed -i "1i\SNP\tbeta\tse" $i.$line.Ileum.trait

echo $line
done

echo -e "$i +++++"

done


```

extract GWAS

```
for n in /groups/umcg-gastrocol/tmp01/Shixian/GeneticCorrelation/Mental_disease/eQTL/Coloc/Locus*

do

locus=$(basename $n)

cd $n

for i in //groups/umcg-gastrocol/tmp01/Shixian/GeneticCorrelation/Mental_disease/GWAS_beta/*beta.se.txt

do

gwas=$(basename $i)
gwas=$(echo $gwas | cut -f1 -d".")
awk 'FNR==NR { a[$1]=$1; next } $1 in a { print $0 }' /groups/umcg-gastrocol/tmp01/Shixian/GeneticCorrelation/Mental_disease/LAVA_locus/$locus.txt $i > $locus.$gwas.trait
sed -i "1i\SNP\tbeta\tse" $locus.$gwas.trait
sed -i "s/ /\t/g" $locus.$gwas.trait

echo -e "$locus +++ $gwas"
done

done



```

harmonize all traits

```
cd  Locus.11.1733

Rscript ----------

options = commandArgs(trailingOnly = TRUE)
name1=options[1]
name2=options[2]

.libPaths( c( .libPaths(), "/groups/umcg-weersma/tmp01/Shixian/Rpackage") )

# args[1] is beta.file
# args[2] is se.file

require(plyr)
library(data.table)

temp = list.files(pattern="*.trait")
for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i],header = T,stringsAsFactors = F,sep = "\t",row.names = NULL))

trait_list=list()
for(i in 1:length(temp)){
  mm=as.data.frame(get(temp[i]))
  mm=mm[!duplicated(mm$SNP),]
  rownames(mm)=mm$SNP
  mm$SNP=NULL
  colnames(mm)=paste(colnames(mm),temp[i],sep = "_")
  mm$rn=rownames(mm)
  trait_list=append(trait_list,list(mm))
}

aa <- join_all(trait_list, by = 'rn', type = 'full')
rownames(aa)=aa$rn
aa$rn=NULL

beta=aa[,colnames(aa) %like% "beta_Locus"]
se=aa[,colnames(aa) %like% "se_Locus"]

write.table(beta,file = name1,row.names = T,sep = "\t",quote = F)
write.table(se,file = name2,row.names = T,sep = "\t",quote = F)

```

***HyprColoc***


```

library(foreach)
trait1="CD"
trait2="Anorexia_nervosa"
traitN=temp[temp %like% "ENSG"]

result=foreach(i=1:length(traitN),.combine = rbind) %do%  {
  tmp.trait=traitN[i]
  
  tmp.beta=beta[,c(colnames(beta)[colnames(beta) %like% trait1],
                   colnames(beta)[colnames(beta) %like% trait2],
                   colnames(beta)[colnames(beta) %like% tmp.trait]),drop=F]
  tmp.se=se[,c(colnames(se)[colnames(se) %like% trait1],
                   colnames(se)[colnames(se) %like% trait2],
                   colnames(se)[colnames(se) %like% tmp.trait]),drop=F]
  tmp.beta=na.omit(tmp.beta)
  tmp.se=na.omit(tmp.se)
  
  tmp.snp=intersect(rownames(tmp.beta),rownames(tmp.se))
  tmp.beta=tmp.beta[rownames(tmp.beta) %in% tmp.snp,]
  tmp.se=tmp.se[rownames(tmp.se) %in% tmp.snp,]
  
  colnames(tmp.beta)=c(trait1,trait2,tmp.trait)
  colnames(tmp.se)=c(trait1,trait2,tmp.trait)
  
  traits <- c(trait1,trait2,tmp.trait)
  rsid <- rownames(tmp.snp)
  
  mm=hyprcoloc(as.matrix(tmp.beta), as.matrix(tmp.se), trait.names=traits, snp.id=rsid)
  aa=as.data.frame(mm$results)
  aa$GeneTrait=tmp.trait
  
  return.string=aa
  
}


```


***MR***

```
ml RPlus/4.2.1-foss-2022a-v22.12.1
ml GCC/7.3.0-2.30

file="GWAS_HSERMETANA.PWY..L.methionine.biosynthesis.III.rearrange.tsv.gz"
name=$(basename $file)
name=${name%.tsv*}

zcat $file | awk 'OFS="\t"{if($1!="NA")print $3,$1,$4,7899,$5,$6,$7,$8,$9,5,$2}' | awk '!seen[$2]++'> $name.tsv
sed -i "1s/.*/chr rs  ps  n_miss  allel1  allel0  af  beta  se  l_remle p_wald/" $name.tsv

Rscript TwoSampleMR.R $name.tsv CD.txt

```

```
args <- commandArgs(trailingOnly = TRUE)
.libPaths("/groups/umcg-weersma/tmp01/Shixian/Rpackage/")

library(TwoSampleMR)
library(MRInstruments)
library(ieugwasr)

test_exposure=read.table(args[1],header = T,stringsAsFactors = F)
test_outcome=read.table(args[2],header = T,stringsAsFactors = F)
test_outcome$NA.=NULL

print("+++ Read files done")

exposure.name=tools::file_path_sans_ext(basename(args[1]))
outcome.name=tools::file_path_sans_ext(basename(args[2]))

colnames(test_exposure)=c("chr","SNP","ps","n_miss","effect_allele","other_allele","eaf","beta","se","l_remle","P")
test_exposure$Phenotype=exposure.name
dat_exp <- format_data(test_exposure, type="exposure")
colnames(test_outcome)=c("SNP","effect_allele","other_allele","eaf","beta","se","P")
test_outcome$Phenotype=outcome.name
dat_out <- format_data(test_outcome, type="outcome")

print("+++ Change format done")

# clumping
dat_clump <- ld_clump(
  dplyr::tibble(rsid=dat_exp$SNP, pval=dat_exp$pval.exposure, id=dat_exp$exposure),
  plink_bin = "/groups/umcg-weersma/tmp01/Shixian/Tools/plink",
  bfile = "/groups/umcg-gastrocol/tmp01/Shixian/GeneticCorrelation/Metabolic_disease/reference/EUR"
)
dat_exp=dat_exp[dat_exp$SNP %in% dat_clump$rsid,]
dat_exp=dat_exp[dat_exp$pval.exposure<5e-6,]

print("+++ Clumping done")

# harmonize
dat_harm <- harmonise_data(
  exposure_dat = dat_exp,
  outcome_dat = dat_out
)
dat_harm$mr_keep="TRUE"
dat_harm$mr_keep=as.logical(dat_harm$mr_keep)

# run analysis
test_MR <- mr(dat_harm,method_list=c("mr_egger_regression", "mr_ivw","mr_wald_ratio","mr_weighted_median"))
test_hetero=mr_heterogeneity(dat_harm)
test_pleiotropy=mr_pleiotropy_test(dat_harm)

print("+++ Analysis done")

write.table(test_MR,file = paste(outcome.name,exposure.name,"MR.txt",sep = "."),sep = "\t",quote = F,row.names = F)
write.table(test_hetero,file = paste(outcome.name,exposure.name,"heterogeneity.txt",sep = "."),sep = "\t",quote = F,row.names = F)
write.table(test_pleiotropy,file = paste(outcome.name,exposure.name,"pleiotropy.txt",sep = "."),sep = "\t",quote = F,row.names = F)
```

***SMR***
```

#!/bin/bash

helpFunction()
{
   echo ""
   echo "Usage: $0 -a parameterA -b parameterB -c parameterC"
   echo -e "\t-a parameterA absolute location of GWAS data"
   echo -e "\t-b parameterB absolute location of QTL data"
   echo -e "\t-c parameterC absolute location of output file"
   exit 1 # Exit script after printing help
}

while getopts "a:b:c:" opt
do
   case "$opt" in
      a ) parameterA="$OPTARG" ;;
      b ) parameterB="$OPTARG" ;;
      c ) parameterC="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$parameterA" ] || [ -z "$parameterB" ] || [ -z "$parameterC" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# Begin script in case all parameters are correct
echo -e "GWAS data ***** $parameterA"
echo -e "QTL data ***** $parameterB"
echo -e "Output file ***** $parameterC"


ml PLINK/1.9-beta6-20190617
/groups/umcg-gastrocol/tmp01/Shixian/Tools/smr_v1.3.1_linux_x86_64_static/smr_v1.3.1_linux_x86_64_static \
--bfile /groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Reference/g1000_eur \
--gwas-summary $parameterA \
--beqtl-summary $parameterB \
--peqtl-smr 1e-05 --thread-num 10 \
--diff-freq-prop 0.9 --diff-freq 1 \
--out $parameterC \
--smr-multi
```

***CELLECT***

```
cd /groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Tools/CELLEX/cellex/
python
import sys
sys.path.append("/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Tools/CELLEX/cellex")

import numpy as np # needed for formatting data for this tutorial
import pandas as pd # needed for formatting data for this tutorial
import cellex
```

```
Default python version Python 3.11.3 

cd /groups/umcg-griac/tmp01/projects/umcg-cqi/geneticcorrelation/tools/tmp.conda/bin
./conda info --envs
conda activate snakemake

# Note: pybedtools sometimes needs update which could lead to "no module named pybedtools", to avoid this:
cd pybedtools
python setup.py develop
python setup.py cythonize develop
see https://daler.github.io/pybedtools/main.html

Then follow the tutorial https://github.com/perslab/CELLECT/wiki/CELLECT-LDSC-Tutorial and https://github.com/perslab/CELLECT/wiki/CELLECT-MAGMA-Tutorial
```



