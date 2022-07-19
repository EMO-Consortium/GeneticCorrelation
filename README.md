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




