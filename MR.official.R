
library(TwoSampleMR)
library(ieugwasr)
library(RadialMR)
library(MRPRESSO)

# test_exposure is the exposure data, columns follow c("SNP","effect_allele","other_allele","eaf","beta","se","P")
# test_outcome is the outcome data, columns follow c("SNP","effect_allele","other_allele","eaf","beta","se","P")
# exposure.name is the exposure trait
# outcome.name is the outocme trait

#  *********************************** format prepare ***********************************

test_exposure$Phenotype=exposure.name
dat_exp <- format_data(test_exposure, type="exposure")
test_outcome$Phenotype=outcome.name
dat_out <- format_data(test_outcome, type="outcome")

print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++ Change format done")

#  *********************************** IVs selection ***********************************
# there are different methods and assumptions to check IVs, all the following steps are suggested, but not limited in these

# (1) independent: clumpling with r2 <0.02, genomic distance >1000Kb, MAF >0.05
# (2) significant with exposure, P <5e-08
# these all default setting in dat_clump function

dat_clump <- ld_clump(
  dplyr::tibble(rsid=dat_exp$SNP, pval=dat_exp$pval.exposure, id=dat_exp$exposure),
  plink_bin = "/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Tools/plink",
  bfile = "/groups/umcg-griac/tmp01/projects/umcg-cqi/GeneticCorrelation/Reference/EUR"
)
dat_exp=dat_exp[dat_exp$SNP %in% dat_clump$rsid,]
dat_exp=dat_exp[dat_exp$pval.exposure<5e-8,]

print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++ Clumping done")

# (3) harmonize and choose proxy, or use extract_outcome_data (with default proxy parameter)
dat_harm <- harmonise_data(
  exposure_dat = dat_exp,
  outcome_dat = dat_out
)
dat_harm$mr_keep="TRUE"
dat_harm$mr_keep=as.logical(dat_harm$mr_keep)

print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++ harmonize done")

# (4) outlier removing, heterogenious SNPs are removed by ivw (meta-based) and egger (regression-based) methods
mod_outlier=format_radial(dat_harm$beta.exposure,dat_harm$beta.outcome,dat_harm$se.exposure,dat_harm$se.outcome,RSID = dat_harm$SNP)
mod_radial_ivw=ivw_radial(mod_outlier,alpha = 0.05,weights = 1,tol = 0.0001)
mod_radial_egger=egger_radial(mod_outlier,alpha = 0.05,weights = 1)
outlier=c(as.character(mod_radial_ivw$outliers$SNP),as.character(mod_radial_egger$outliers$SNP))
dat_harm=dat_harm[!dat_harm$SNP %in% outlier,]
write.table(outlier,file = paste("OutToExp",outcome.name,exposure.name,"Outliers.test.txt",sep = "."),sep = "\t")

print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++ outlier done")

# (5) confounder removing
# remove BMI, alcohol, smoking related SNPs, check phenoscan database
write.table(dat_harm,file = paste("OutToExp",outcome.name,exposure.name,"Final.used.SNP.test.txt",sep = "."))

# (6) Steiger test, check whether IVs explain more variation in exposure comapred with outocmes 
tryCatch({
  mod_steiger=mr_steiger(p_exp=dat_harm$pval.exposure, p_out=dat_harm$pval.outcome, 
                         n_exp=50000, n_out=50000, r_exp=dat_harm$beta.exposure, r_out=dat_harm$beta.outcome)
  output_steiger=data.frame(Exposure=exposure.name,Outcome=outcome.name,
                            R2_exposure=mod_steiger$r2_exp,R2_outcome=mod_steiger$r2_out,
                            P_steiger=mod_steiger$steiger_test,HypothesisNULL=mod_steiger$correct_causal_direction)
  write.table(output_steiger,file = paste("OutToExp",outcome.name,exposure.name,"Steiger.test.txt",sep = "."),sep = "\t")
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++ Steiger done")

# (7) F-statistic (final validation)

tryCatch({
  BXG=abs(dat_harm$beta.exposure)
  seBetaXG=dat_harm$se.exposure
  F   = BXG^2/seBetaXG^2
  mF  = mean(F)
  
  output_F<-cbind(exposure.name, mF)
  colnames(output_F) <- c("Exposure", "F value")
  write.table(output_F,file = paste("OutToExp",outcome.name,exposure.name,"F.test.txt",sep = "."),sep = "\t")
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++ F test done")

#  *********************************** run MR ***********************************
test_MR <- mr(dat_harm,method_list=c("mr_egger_regression", "mr_ivw","mr_wald_ratio","mr_weighted_median"))
write.table(test_MR,file = paste("OutToExp",outcome.name,exposure.name,"MR.txt",sep = "."),sep = "\t")

print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++ MR run done")

#  *********************************** sensitivity ***********************************

tryCatch({
  test_hetero=mr_heterogeneity(dat_harm)
  write.table(test_hetero,file = paste("OutToExp",outcome.name,exposure.name,"heterogeneity.test.txt",sep = "."),sep = "\t")
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++ Qtest done")

tryCatch({
  test_pleiotropy=mr_pleiotropy_test(dat_harm)
  write.table(test_pleiotropy,file = paste("OutToExp",outcome.name,exposure.name,"pleiotropy.test.txt",sep = "."),sep = "\t")
  
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++ Pleiotropy ME Egger  done")

tryCatch({
  test_MRpresso=mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", 
                          SdExposure = "se.exposure", OUTLIERtest = T, DISTORTIONtest = T, data = dat_harm, NbDistribution = 1000,  
                          SignifThreshold = 0.05)
  output_MRPresso_global=test_MRpresso$`MR-PRESSO results`$`Global Test`$Pvalue
  write.table(output_MRPresso_global,file = paste("OutToExp",outcome.name,exposure.name,"MR_Presso.global.test.txt",sep = "."),sep = "\t")
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++ Pleiotropy MR PPRESSO done")

tryCatch({
  res_loo <- mr_leaveoneout(dat_harm,method=mr_ivw)
  output_leave_one_out=data.frame(SNP=res_loo$SNP,Beta=res_loo$b,SE=res_loo$se,Pvalue=res_loo$p)
  write.table(output_leave_one_out,file = paste("OutToExp",outcome.name,exposure.name,"Leave_one_out.test.txt",sep = "."),sep = "\t")
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++ Leave one out done")
