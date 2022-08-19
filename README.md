# GeneticCorrelation

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

***Split locus***

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
cd  Locus.11.1733

/groups/umcg-gastrocol/tmp01/Shixian/Tools/smr_v1.3.1_linux_x86_64_static/smr_v1.3.1_linux_x86_64_static \
--beqtl-summary /groups/umcg-gastrocol/tmp01/Shixian/GeneticCorrelation/Mental_disease/eQTL/eQTLGen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense \
--extract-snp Locus.11.1733.txt \
--query 1 \
--out Locus.11.1733.eQTLGen 

```

extract GWAS

```
cd  Locus.11.1733

awk 'FNR==NR { a[$1]=$1; next } $1 in a { print $0 }'  Locus.11.1733.txt ../../../GWAS_beta/CD.beta.se.txt > Locus.11.1733.CD.trait

sed -i "1i\SNP\tbeta\tse" Locus.11.1733.CD.trait
sed -i "s/ /\t/g" Locus.11.1733.CD.trait

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
for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i],header = T,stringsAsFactors = F,sep = "\t",row.names = 1))

trait_list=list()
for(i in 1:length(temp)){
  mm=as.data.frame(get(temp[i]))
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






