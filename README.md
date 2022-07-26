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









