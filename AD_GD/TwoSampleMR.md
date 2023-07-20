## Setting up
* Insatll packages
```R
## R version: R/4.2.2-foss-2022a-bare
.libPaths("/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Tools/Rpackage")

```
## Start MR
* Buil R script
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

colnames(test_exposure)=c("SNP","effect_allele","other_allele","eaf","beta","se","P")
test_exposure$Phenotype=exposure.name
dat_exp <- format_data(test_exposure, type="exposure")
colnames(test_outcome)=c("SNP","effect_allele","other_allele","eaf","beta","se","P")
test_outcome$Phenotype=outcome.name
dat_out <- format_data(test_outcome, type="outcome")

### clumping
dat_clump <- ld_clump(
  dplyr::tibble(rsid=dat_exp$SNP, pval=dat_exp$pval.exposure, id=dat_exp$exposure),
  plink_bin = "/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Tools/plink",
  bfile = "/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Reference/EUR"
)
dat_exp=dat_exp[dat_exp$SNP %in% dat_clump$rsid,]
dat_exp=dat_exp[dat_exp$pval.exposure<5e-8,]

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

*Generate bash jobs
```bash
while IFS= read -r file1; do
    while IFS= read -r file2; do
        # 在这里执行你的操作，使用 $file1 和 $file2
        echo "处理文件：$file1 和 $file2"
    done < gwas_sum_gut.txt
done < gwas_sum_lung.txt
```

*Merge results
```R
## load results
file_names <- list.files("./Out_D1", pattern = "*.MR.txt", full.names = TRUE)
combined_data <- data.frame()
combined_data <- do.call(rbind, lapply(file_names, function(file) {
  read.table(file, sep = "\t", header = TRUE)}))

file_names <- list.files("./Out_D1", pattern = "*.heterogeneity.txt", full.names = TRUE)
heter <- data.frame()
heter <- do.call(rbind, lapply(file_names, function(file) {
  read.table(file, sep = "\t", header = TRUE)}))

file_names <- list.files("./Out_D1", pattern = "*.pleiotropy.txt", full.names = TRUE)
plei <- data.frame()
plei <- do.call(rbind, lapply(file_names, function(file) {
  read.table(file, sep = "\t", header = TRUE)}))

## Merge results
combined_data$matchid<-paste0(combined_data$outcome,"_",combined_data$exposure)

heter1<-heter[which(heter$method=="MR Egger"),]
heter1$matchid<-paste0(heter1$outcome,"_",heter1$exposure)

heter2<-heter1[match(combined_data$matchid,heter1$matchid),]
df.res1<-cbind(combined_data,heter2[,c("Q","Q_df","Q_pval")])

combined_data$matchid2<-paste0(combined_data$outcome,"_",combined_data$exposure)
plei$matchid2<-paste0(plei$outcome,"_",plei$exposure)
plei1<-plei[match(combined_data$matchid2,plei$matchid2),]
df.res<-cbind(df.res1,plei1[,c("egger_intercept","se","pval")])

write.csv(df.res,file="MR_results_merge1_D1.csv")

df.egger<-df.res[which(df.res$method=="MR Egger"),]

```

