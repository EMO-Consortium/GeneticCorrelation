

```bash
ml PLINK/1.9-beta6-20190617

/groups/umcg-griac/tmp01/projects/umcg-cqi/software/smr-1.3.1-linux-x86_64/smr-1.3.1 \
--bfile /groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Reference/g1000_eur \
--gwas-summary /groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/GWAS_clean/CD.txt \
--beqtl-summary /groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/eQTLGen/Candidate.SMR.eQTLGen \
--out WDR18.GWAS.eQTL --plot --probe ENSG00000065268 \
--probe-wind 500 --gene-list /groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Reference/glist-hg19 \
--peqtl-smr 1e-05 --thread-num 10 --diff-freq-prop 0.9 --diff-freq 1

smr --bfile mydata --gwas-summary mygwas.ma --beqtl-summary myeqtl --out myplot --plot --probe rs2365709 --probe-wind 500

```
