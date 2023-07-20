## Export data (examples) for SMR plot 
* load plink
```bash
ml PLINK/1.9-beta6-20190617
```
* Example1 ENSG00000090621
```bash
### extract eQTL ENSG00000090621 blood
/groups/umcg-griac/tmp01/projects/umcg-cqi/software/smr-1.3.1-linux-x86_64/smr-1.3.1 \
--beqtl-summary /groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/eQTLGen/Candidate.SMR.eQTLGen \
--query 0.99 \
--probe ENSG00000090621 \
--out eQTL.blood.PABPC4

### extract GWAS IBS FEV1_to_FVC
awk 'NR==FNR { data[$1]++; next } $1 in data' eQTL.PABPC4.txt /groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/GWAS_clean/IBS.txt > IBS.blood.PABPC4.txt
awk 'NR==FNR { data[$1]++; next } $1 in data' eQTL.PABPC4.txt /groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/GWAS_clean/FEV1_to_FVC.txt > FEV1_to_FVC.blood.PABPC4.txt

### extract eQTL ENSG00000090621 Colon_Sigmoid
/groups/umcg-griac/tmp01/projects/umcg-cqi/software/smr-1.3.1-linux-x86_64/smr-1.3.1 \
--beqtl-summary /groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/eQTLTissue/Candidate.SMR.Colon_Sigmoid \
--query 0.99 \
--probe ENSG00000090621 \
--out eQTL.Colon_Sigmoid.PABPC4

### extract GWAS IBS FEV1_to_FVC
awk 'NR==FNR { data[$1]++; next } $1 in data' eQTL.Colon_Sigmoid.PABPC4.txt /groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/GWAS_clean/IBS.txt > IBS.Colon_Sigmoid.PABPC4.txt
awk 'NR==FNR { data[$1]++; next } $1 in data' eQTL.Colon_Sigmoid.PABPC4.txt /groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/GWAS_clean/FEV1_to_FVC.txt > FEV1_to_FVC.Colon_Sigmoid.PABPC4.txt

```
* Example 2 ENSG00000214944
```bash
### extract eQTL ENSG00000214944
/groups/umcg-griac/tmp01/projects/umcg-cqi/software/smr-1.3.1-linux-x86_64/smr-1.3.1 \
--beqtl-summary /groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/eQTLTissue/Candidate.SMR.Esophagus_Muscularis \
--query 0.99 \
--probe ENSG00000214944 \
--out eQTL.Esophagus_Muscularis.ARHGEF28

### extract GWAS asthma Diverticular_disease
awk 'NR==FNR { data[$1]++; next } $1 in data' eQTL.Esophagus_Muscularis.ARHGEF28.txt /groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/GWAS_clean/Asthma.txt > Asthma.Esophagus_Muscularis.ARHGEF28.txt
awk 'NR==FNR { data[$1]++; next } $1 in data' eQTL.Esophagus_Muscularis.ARHGEF28.txt /groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/GWAS_clean/Diverticular_disease.txt > Diverticular_disease.Esophagus_Muscularis.ARHGEF28.txt

```
* Example 3 ENSG00000197375
```bash
### extract eQTL ENSG00000197375
/groups/umcg-griac/tmp01/projects/umcg-cqi/software/smr-1.3.1-linux-x86_64/smr-1.3.1 \
--beqtl-summary /groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/eQTLTissue/Candidate.SMR.Esophagus_Muscularis \
--query 0.99 \
--probe ENSG00000197375 \
--out eQTL.Esophagus_Muscularis.SLC22A5

### extract GWAS asthma Diverticular_disease
awk 'NR==FNR { data[$1]++; next } $1 in data' eQTL.Esophagus_Muscularis.SLC22A5.txt /groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/GWAS_clean/Asthma.txt > Asthma.Esophagus_Muscularis.SLC22A5.txt
awk 'NR==FNR { data[$1]++; next } $1 in data' eQTL.Esophagus_Muscularis.SLC22A5.txt /groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/GWAS_clean/Gastro_reflux.txt > Gastro_reflux.Esophagus_Muscularis.SLC22A5.txt

```




