
options = commandArgs(trailingOnly = TRUE)
eqtl_file=options[1]
mbqtl_path=options[2]
name=options[3]

.libPaths( c( .libPaths(), "/groups/umcg-weersma/tmp01/Shixian/Rpackage") )

library(coloc)
library(data.table)
library(foreach)

eqtl=read.table(eqtl_file,header=T,stringsAsFactors = F)

summary.coloc=matrix(ncol = 8)
colnames(summary.coloc)=c("nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf","taxa","Probe")

temp = list.files(path = mbqtl_path,pattern="*rearrange.tsv.gz")
coloc.result = foreach(i=1:length(temp),.combine = rbind) %do%  {
  tmp.taxa=temp[i]
  mbqtl= fread(paste(mbqtl_path,tmp.taxa,sep="/"))

  overlap=intersect(eqtl$SNP,mbqtl$variant_id)

  if(length(overlap)>50){
    tmp.data1=eqtl[eqtl$SNP %in% overlap,]
    tmp.data2=mbqtl[mbqtl$variant_id %in% overlap,]

    tmp.data1=tmp.data1[order(tmp.data1$SNP),]
    tmp.data2=tmp.data2[order(tmp.data2$variant_id),]

    mm=coloc.abf(dataset1=list(pvalues=tmp.data1$p,N=1000,type="quant"),
                 dataset2=list(pvalues=tmp.data2$p_value,N=7749,type="quant"),
                 MAF=tmp.data2$effect_allele_frequency)

    tmp.result=as.data.frame(t(mm$summary))
    tmp.result$taxa=tmp.taxa
    tmp.result$probe=eqtl_file

    return.string=tmp.result
  }
}

write.table(coloc.result,file = paste(eqtl_file,name,"coloc.txt",sep = "_"),quote=F,row.names=F)


