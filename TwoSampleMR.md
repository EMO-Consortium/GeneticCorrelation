## Setting up
* Insatll packages
```R
## R version: R/4.2.2-foss-2022a-bare
.libPaths("/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Tools/Rpackage")

```
## Start MR
```R
.libPaths("/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Tools/Rpackage")
library(TwoSampleMR)
#library(MRInstruments)
library(ieugwasr)

### test one pair, load data
test_exposure=read.table("/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/GWAS_clean/FEV1_to_FVC.txt",header = T,stringsAsFactors = F)
test_outcome=read.table("/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/GWAS_clean/CD.txt",header = T,stringsAsFactors = F)
test_outcome$NA.=NULL

### re-format
exposure.name=tools::file_path_sans_ext(basename("/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/GWAS_clean/FEV1_to_FVC.txt"))
outcome.name=tools::file_path_sans_ext(basename("/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/GWAS_clean/CD.txt"))

colnames(test_exposure)=c("SNP","effect_allele","other_allele","eaf","beta","se","P","Phenotype")
test_exposure$Phenotype=exposure.name
dat_exp <- format_data(test_exposure, type="exposure")
colnames(test_outcome)=c("SNP","effect_allele","other_allele","eaf","beta","se","P","Phenotype")
test_outcome$Phenotype=outcome.name
dat_out <- format_data(test_outcome, type="outcome")

### clumping
dat_clump <- ld_clump(
  dplyr::tibble(rsid=dat_exp$SNP, pval=dat_exp$pval.exposure, id=dat_exp$exposure),
  plink_bin = "/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Tools/plink",
  bfile = "/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Reference/EUR"
)
dat_exp=dat_exp[dat_exp$SNP %in% dat_clump$rsid,]
dat_exp=dat_exp[dat_exp$pval.exposure<5e-6,]

### harmonize
dat_harm <- harmonise_data(
  exposure_dat = dat_exp,
  outcome_dat = dat_out
)
dat_harm$mr_keep="TRUE"
dat_harm$mr_keep=as.logical(dat_harm$mr_keep)

### run analysis
test_MR <- mr(dat_harm,method_list=c("mr_egger_regression", "mr_ivw","mr_wald_ratio","mr_weighted_median"))
test_hetero=mr_heterogeneity(dat_harm)
test_pleiotropy=mr_pleiotropy_test(dat_harm)

print("+++ Analysis done")

write.table(test_MR,file = paste(outcome.name,exposure.name,"MR.txt",sep = "."),sep = "\t",quote = F,row.names = F)
write.table(test_hetero,file = paste(outcome.name,exposure.name,"heterogeneity.txt",sep = "."),sep = "\t",quote = F,row.names = F)
write.table(test_pleiotropy,file = paste(outcome.name,exposure.name,"pleiotropy.txt",sep = "."),sep = "\t",quote = F,row.names = F)


```
