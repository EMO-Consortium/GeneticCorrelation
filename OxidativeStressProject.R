# 
library(edgeR)
library(limma)
library(foreach)
library(ggplot2)
library(foreach)
library(lme4)
library(nlme)
library(factoextra)
library(ggsci)
library(vegan)
require(cowplot)
library(ggpubr)
library(annotate)
library(biomaRt)
library(crayon)
library(metafor)

individual_cal=function(tmp.batch){
  tmp.expression=expression_all[rownames(expression_all) %in% pheno_all$Sample.ID[pheno_all$Data==tmp.batch],]
  tmp.covar=pheno_all[pheno_all$Data==tmp.batch,]
  
  individual_test = foreach(i=1:ncol(expression_all),.combine = rbind) %do%  {
    tmp.gene=colnames(expression_all)[i]
    tmp.data=merge(tmp.expression[,tmp.gene,drop=F],tmp.covar,by.x="row.names",by.y="Sample.ID",all=F)
    tmp.data=tmp.data[,apply(tmp.data,2,function(x){
      length(unique(x))
    }) >1,drop = FALSE]
    
    tmp.data$Inflammation.status[tmp.data$Inflammation.status=="inflamed"]=2
    tmp.data$Inflammation.status[tmp.data$Inflammation.status=="uninflamed"]=1
    tmp.data$Inflammation.status[tmp.data$Inflammation.status=="normal"]=0
    tmp.data$Inflammation.status=as.numeric(tmp.data$Inflammation.status)
    colnames(tmp.data)[2]="Gene"
    
    cat(yellow(tmp.batch,"===",tmp.gene,"\n"))
    tmp.result = as.data.frame(summary(lm(Gene ~ Inflammation.status,data = tmp.data))$coef)
    return.string=data.frame(Data=tmp.batch,Gene=tmp.gene,
                             beta=tmp.result$Estimate[rownames(tmp.result)=="Inflammation.status"],
                             sd=tmp.result$`Std. Error`[rownames(tmp.result)=="Inflammation.status"],
                             pvalue=tmp.result$`Pr(>|t|)`[rownames(tmp.result)=="Inflammation.status"])
    
  }
  return(individual_test)
}

######################################

# data check

######################################
# array data: GSE112366,GSE75214,GSE179285
# rna data: GSE165512,GSE137344,GSE111889

array1=read.table("input/filtered/matrix/GSE112366.os_gene.filtered.revised.txt",row.names = 1,header = T)
array2=read.table("input/filtered/matrix/GSE75214.os_gene.filtered.revised.txt",row.names = 1,header = T)
array3=read.table("input/filtered/matrix/GSE179285.os_gene.filtered.revised.txt",row.names = 1,header = T)
seq1=read.table("input/filtered/matrix/GSE165512.os_gene.filtered revised.txt",row.names = 1,header = T)
seq2=read.table("input/filtered/matrix/GSE137344.os_gene.filtered revised.txt",row.names = 1,header = T)
seq3=read.table("input/filtered/matrix/GSE111889.os_gene.filtered.revised.txt",row.names = 1,header = T)

pheno_array1=read.table("../样本信息/GSE112366.txt",header = T,stringsAsFactors = F,sep = "\t")
pheno_array2=read.table("../样本信息/GSE75214.txt",header = T,stringsAsFactors = F,sep = "\t")
pheno_array3=read.table("../样本信息/GSE179285.txt",header = T,stringsAsFactors = F,sep = "\t")
pheno_seq1=read.table("../样本信息/GSE165512.txt",header = T,stringsAsFactors = F,sep = "\t")
pheno_seq2=read.table("../样本信息/GSE137344.txt",header = T,stringsAsFactors = F,sep = "\t")
pheno_seq3=read.table("../样本信息/GSE111889.txt",header = T,stringsAsFactors = F,sep = "\t")

pheno_array3$Gender=NULL

pheno_array1$Inflammation.status="NA"
pheno_array1$Data="GSE112366"
pheno_array2$Data="GSE75214"
pheno_array3$Data="GSE179285"
pheno_seq1$Data="GSE165512"
pheno_seq2$Data="GSE137344"
pheno_seq3$Data="GSE111889"
pheno_seq1$Inflammation.status="NA"
pheno_array1$Gender="NA"
pheno_array2$Gender="NA"
pheno_array3$Gender="NA"
pheno_seq1$Gender="NA"
pheno_seq2$Inflammation.status="NA"
pheno_seq3$Inflammation.status="NA"

pheno_array1$Inflammation.status[pheno_array1$Samples=="CD"]="inflamed"
pheno_array1$Inflammation.status[pheno_array1$Samples=="HC"]="normal"
pheno_seq1$Inflammation.status[pheno_seq1$Samples=="CD"]="inflamed"
pheno_seq1$Inflammation.status[pheno_seq1$Samples=="HC"]="normal"
pheno_seq2$Inflammation.status[pheno_seq2$Samples=="HC"]="normal"
pheno_seq2$Inflammation.status[pheno_seq2$Samples=="CD"]="inflamed"
pheno_seq3$Inflammation.status[pheno_seq3$Samples=="HC"]="normal"
pheno_seq3$Inflammation.status[pheno_seq3$Samples=="CD"]="inflamed"
pheno_seq3$Tissues[pheno_seq3$Tissues=="rectum"]="colon"

pheno_all=rbind(pheno_array1,pheno_array2,pheno_array3,pheno_seq1,pheno_seq2,pheno_seq3)
expression_all=Reduce(intersect,list(colnames(array1),colnames(array2),colnames(array3),colnames(seq1),colnames(seq2),colnames(seq3)))
array1=array1[,colnames(array1) %in% expression_all]
array2=array2[,colnames(array2) %in% expression_all]
array3=array3[,colnames(array3) %in% expression_all]
seq1=seq1[,colnames(seq1) %in% expression_all]
seq2=seq2[,colnames(seq2) %in% expression_all]
seq3=seq3[,colnames(seq3) %in% expression_all]
expression_all=rbind(array1,array2,array3,seq1,seq2,seq3)

# PCA
pca=prcomp(expression_all,scale = TRUE)
fviz_eig(pca)
eigenvalue=get_eig(pca)
ind <- get_pca_ind(pca)
pca_matrix=as.data.frame(ind$coord)
pca_matrix=merge(pca_matrix[,1:10],pheno_all,by.x="row.names",by.y="Sample.ID",all=F)

ggplot (pca_matrix, aes(Dim.1,Dim.2,color=Data)) + geom_point(size=2) + 
  theme_bw()+scale_color_npg()+
  theme(legend.position = 'top')
ggplot (pca_matrix, aes(Dim.1,Dim.2,color=Tissues)) + geom_point(size=2) + 
  theme_bw()+scale_color_npg()+
  theme(legend.position = 'top')
ggplot (pca_matrix, aes(Dim.1,Dim.2,color=Inflammation.status)) + geom_point(size=2) + 
  theme_bw()+scale_color_npg()+
  theme(legend.position = 'top')

#remove outliers
expression_all=expression_all[!rownames(expression_all) %in% c("CSMDRVXQ","HSM5FZAZ"),]
pheno_all=pheno_all[!pheno_all$Sample.ID %in% c("CSMDRVXQ","HSM5FZAZ"),]

# individual analysis
result_array1=individual_cal("GSE112366")
result_array2=individual_cal("GSE75214")
result_array3=individual_cal("GSE179285")
result_seq1=individual_cal("GSE165512")
result_seq2=individual_cal("GSE137344")
result_seq3=individual_cal("GSE111889")

# heatmap visulization
library(pheatmap)
result_array1$Zscore=result_array1$beta/result_array1$sd
result_array2$Zscore=result_array2$beta/result_array2$sd
result_array3$Zscore=result_array3$beta/result_array3$sd
result_seq1$Zscore=result_seq1$beta/result_seq1$sd
result_seq2$Zscore=result_seq2$beta/result_seq2$sd
result_seq3$Zscore=result_seq3$beta/result_seq3$sd

matrix_all=cbind(result_array1[,"Zscore",drop=F],result_array2[,"Zscore",drop=F],result_array3[,"Zscore",drop=F],
                 result_seq1[,"Zscore",drop=F],result_seq2[,"Zscore",drop=F],result_seq3[,"Zscore",drop=F])
rownames(matrix_all)=(result_array1$Gene)
matrix_all=matrix_all[rownames(matrix_all) %in% meta_cal$Gene[meta_cal$FDR<0.05],]
heatmap_plot=as.data.frame(t(matrix_all))

pdf("output//heatmap.variation.pdf",height = 3,width = 10)
pheatmap(heatmap_plot,cluster_cols = T, cluster_rows = F,scale = "row",
         show_rownames=T, show_colnames=F, 
        border_color=F,
        color = colorRampPalette(c("#FDB366", "white", "#762A83"))(150))
dev.off()


# meta analysis
result_array1=result_array1[order(result_array1$Gene),]
result_array2=result_array2[order(result_array2$Gene),]
result_array3=result_array3[order(result_array3$Gene),]
result_seq1=result_seq1[order(result_seq1$Gene),]
result_seq2=result_seq2[order(result_seq2$Gene),]
result_seq3=result_seq3[order(result_seq3$Gene),]


meta_cal = foreach(i=1:ncol(expression_all),.combine = rbind) %do%  {
  tmp.gene=colnames(expression_all)[i]
  vec_beta=c(result_array1$beta[result_array1$Gene==tmp.gene],
             result_array2$beta[result_array2$Gene==tmp.gene],
             result_array3$beta[result_array3$Gene==tmp.gene],
             result_seq1$beta[result_seq1$Gene==tmp.gene],
             result_seq2$beta[result_seq2$Gene==tmp.gene],
             result_seq3$beta[result_seq3$Gene==tmp.gene])
  vec_sd=c(result_array1$sd[result_array1$Gene==tmp.gene],
             result_array2$sd[result_array2$Gene==tmp.gene],
             result_array3$sd[result_array3$Gene==tmp.gene],
             result_seq1$sd[result_seq1$Gene==tmp.gene],
             result_seq2$sd[result_seq2$Gene==tmp.gene],
             result_seq3$sd[result_seq3$Gene==tmp.gene])
  tmp.meta=summary(rma(yi = vec_beta,  sei = vec_sd,method = "FE"))
  
  cat(green(i,"===",tmp.gene,"\n"))
  
  return.string=data.frame(Gene=tmp.gene,Estimate=tmp.meta$beta,se=tmp.meta$se,
                           zvalue=tmp.meta$zval,pvalue=tmp.meta$pval,ci.lb=tmp.meta$ci.lb,ci.ub=tmp.meta$ci.ub,
                           I2=tmp.meta$I2,H2=tmp.meta$H2,Q=tmp.meta$QEp)
  
}
meta_cal$FDR=p.adjust(meta_cal$pvalue,method = "BH")
nrow(meta_cal[meta_cal$FDR<0.05,])
#write.table(meta_cal,file = "output/meta.results.txt",row.names = F,quote = F,sep = "\t")
meta_cal=read.table("output/meta.results.txt",header = T,stringsAsFactors = F)

meta_z = foreach(i=1:ncol(expression_all),.combine = rbind) %do%  {
  tmp.gene=colnames(expression_all)[i]
  vec_p=c(result_array1$pvalue[result_array1$Gene==tmp.gene],
             result_array2$pvalue[result_array2$Gene==tmp.gene],
             result_array3$pvalue[result_array3$Gene==tmp.gene],
             result_seq1$pvalue[result_seq1$Gene==tmp.gene],
             result_seq2$pvalue[result_seq2$Gene==tmp.gene],
             result_seq3$pvalue[result_seq3$Gene==tmp.gene])
  vec_w=c(160,100,200,120,150,60)
  tmp.meta=sumz(vec_p,weights = vec_w)
  
  cat(green(i,"===",tmp.gene,"\n"))
  
  return.string=data.frame(Gene=tmp.gene,Zvalue=tmp.meta$z,Pvalue=tmp.meta$p)
  
}
meta_z$FDR=p.adjust(meta_z$Pvalue,method = 'BH')
nrow(meta_z[meta_z$FDR<0.05,])

write.table(result_array1,"output/result.arrray1.txt",sep = "\t",row.names = F,quote = F)
write.table(result_array2,"output/result.arrray2.txt",sep = "\t",row.names = F,quote = F)
write.table(result_array3,"output/result.arrray3.txt",sep = "\t",row.names = F,quote = F)
write.table(result_seq1,"output/result.seq1.txt",sep = "\t",row.names = F,quote = F)
write.table(result_seq2,"output/result.seq2.txt",sep = "\t",row.names = F,quote = F)
write.table(result_seq3,"output/result.seq3.txt",sep = "\t",row.names = F,quote = F)
write.table(meta_z,"output/result.meta.z.txt",sep = "\t",row.names = F,quote = F)

library(biomaRt)
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position')
genes=c(as.character(meta_cal$Gene))
G_list<- getBM(attributes=attributes, filters="hgnc_symbol",values=genes,
               mart=mart, uniqueRows=T)
write.table(G_list,file = "output/Gene.list.Ensemble.txt",row.names = F,quote = F,sep = "\t")

# plot
tmp.gene="PRKAB1"
tmp.data=merge(expression_all[,tmp.gene,drop=F],pheno_all,by.x="row.names",by.y="Sample.ID",all=F)
tmp.data$Inflammation.status=factor(tmp.data$Inflammation.status,levels = c("inflamed","uninflamed","normal"))
ggplot(tmp.data, aes(x=Inflammation.status,y=PRKAB1)) + 
  stat_boxplot(geom = "errorbar",linetype=1)+
  geom_boxplot(aes(fill=Inflammation.status),position=position_dodge(2),outlier.shape = NA)+
  theme_bw()+
  geom_jitter(width = 0.2,aes(fill=Inflammation.status),size=2,shape=21)+
  scale_fill_manual(values=c("#BB5566","#DDAA33","#009988"))+
  scale_color_manual(values=c("#BB5566","#DDAA33","#009988"))+
  guides(color=FALSE)+guides(fill=FALSE)+
  facet_wrap(. ~ Data, scales="free_y")+
  theme(strip.background =element_rect(fill="white"))

#=====================
# volcano  plot
#=====================
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggtext)

meta_plot=merge(meta_cal,G_list,by.x="Gene",by.y="hgnc_symbol",all=F)
meta_plot$position=(meta_plot$end_position-meta_plot$start_position)/2+meta_plot$start_position

meta_plot$threshold=NA
meta_plot$threshold[meta_plot$FDR<0.05]="Significant"
meta_plot$threshold[meta_plot$FDR>0.05]="Non-Significant"
ggplot(meta_plot) +
  geom_point(aes(x=Estimate, y=-log10(pvalue), fill=threshold),
             shape = 21, color = "black", size = 4) +
  xlab("Effect Size") +
  ylab("-log10 p-value") +
  geom_hline(yintercept=-log10(2.9e-02),linetype="dashed", color = "black")+
  geom_vline(xintercept=0,linetype="dashed", color = "grey")+
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) +xlim(-1,1)+
  theme_bw()+scale_fill_manual(values = c("grey90","darkred"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave("output//Volcano..pdf",width = 6,height = 4)

meta_plot$logP=-log10(meta_plot$Pvalue)
meta_plot=meta_plot[!meta_plot$chromosome_name %like% "CHR",]
ggplot(meta_plot,aes(position, y = logP)) +
  geom_point(aes(color=chromosome_name)) + theme_classic()+
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+xlab("")+
  geom_hline(yintercept=0)+
  geom_hline(yintercept=1.5,linetype="dashed", 
             color = "red", size=0.5)+
  geom_segment(aes(x = max(position[chromosome_name==1]), xend = max(position[chromosome_name==1]),y = -0.2, yend = 0.2))+
  geom_segment(aes(x = max(position[chromosome_name==2]), xend = max(position[chromosome_name==2]),y = -0.2, yend = 0.2))+
  geom_segment(aes(x = max(position[chromosome_name==3]), xend = max(position[chromosome_name==3]),y = -0.2, yend = 0.2))+
  geom_segment(aes(x = max(position[chromosome_name==4]), xend = max(position[chromosome_name==4]),y = -0.2, yend = 0.2))+
  geom_segment(aes(x = max(position[chromosome_name==5]), xend = max(position[chromosome_name==5]),y = -0.2, yend = 0.2))+
  geom_segment(aes(x = max(position[chromosome_name==6]), xend = max(position[chromosome_name==6]),y = -0.2, yend = 0.2))+
  geom_segment(aes(x = max(position[chromosome_name==7]), xend = max(position[chromosome_name==7]),y = -0.2, yend = 0.2))+
  geom_segment(aes(x = max(position[chromosome_name==8]), xend = max(position[chromosome_name==8]),y = -0.2, yend = 0.2))+
  geom_segment(aes(x = max(position[chromosome_name==9]), xend = max(position[chromosome_name==9]),y = -0.2, yend = 0.2))+
  geom_segment(aes(x = max(position[chromosome_name==10]), xend = max(position[chromosome_name==10]),y = -0.2, yend = 0.2))+
  geom_segment(aes(x = max(position[chromosome_name==11]), xend = max(position[chromosome_name==11]),y = -0.2, yend = 0.2))+
  geom_segment(aes(x = max(position[chromosome_name==12]), xend = max(position[chromosome_name==12]),y = -0.2, yend = 0.2))+
  geom_segment(aes(x = max(position[chromosome_name==13]), xend = max(position[chromosome_name==13]),y = -0.2, yend = 0.2))+
  geom_segment(aes(x = max(position[chromosome_name==14]), xend = max(position[chromosome_name==14]),y = -0.2, yend = 0.2))+
  geom_segment(aes(x = max(position[chromosome_name==15]), xend = max(position[chromosome_name==15]),y = -0.2, yend = 0.2))+
  geom_segment(aes(x = max(position[chromosome_name==16]), xend = max(position[chromosome_name==16]),y = -0.2, yend = 0.2))+
  geom_segment(aes(x = max(position[chromosome_name==17]), xend = max(position[chromosome_name==17]),y = -0.2, yend = 0.2))+
  geom_segment(aes(x = max(position[chromosome_name==18]), xend = max(position[chromosome_name==18]),y = -0.2, yend = 0.2))+
  geom_segment(aes(x = max(position[chromosome_name==19]), xend = max(position[chromosome_name==19]),y = -0.2, yend = 0.2))+
  geom_segment(aes(x = max(position[chromosome_name==20]), xend = max(position[chromosome_name==20]),y = -0.2, yend = 0.2))+
  geom_segment(aes(x = max(position[chromosome_name==21]), xend = max(position[chromosome_name==21]),y = -0.2, yend = 0.2))+
  geom_segment(aes(x = max(position[chromosome_name==22]), xend = max(position[chromosome_name==22]),y = -0.2, yend = 0.2))

#=====================
# cell enrichment  plot
#=====================
library(data.table)

cell_enrich=read.table("output/WebCSEA_All_tissue_cell_type.txt",header = T,sep = "\t")
cell_enrich$FDR=p.adjust(cell_enrich$input_list_combined_p,method = "BH")
select=data.frame(table(cell_enrich$General_cell_type))
select=select[select$Freq>=10,]
cell_enrich=cell_enrich[cell_enrich$General_cell_type %in% select$Var1,]
cell_enrich=cell_enrich[cell_enrich$Tissue %like% "Intestine" | cell_enrich$Tissue %like% "Colon",]
orders=as.character(cell_enrich$General_cell_typ[!duplicated(cell_enrich$General_cell_typ)])
cell_enrich$General_cell_type=factor(cell_enrich$General_cell_type,levels = orders)

library(randomcoloR)
n <- 23
palette <- distinctColorPalette(n)
ggplot(cell_enrich,aes(General_cell_type, input_list_combined_log10p)) +
  geom_jitter(aes(fill=General_cell_type), shape =21,size=4)+
  guides(fill=F)+theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept=1.84,linetype="dashed", color = "black")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+xlab("")+ylab("-log(10) P value of CSEA")+scale_fill_manual(values = palette)
ggsave("output/Cell_specific_enrichment.pdf",width = 6,height = 4)

#=====================
# locus  plot
#=====================

file_ld=read.table("output/locuszoom.JAK2.r2.ld.txt",header = T)
file_eqtl=read.table("output/Candidate.SMR.eQTLGen.txt.txt",header = T)
file_gwas=read.table("output/CD.gwas.txt",header = T)
file_meqtl=read.table("output/Candidate.SMR.meQTL.txt",header = T,stringsAsFactors = F)
file_mbqtl=read.table("output/Candidate.SMR.mbQTL.txt",header = T,stringsAsFactors = F)

file_plot=file_eqtl
file_plot=file_plot[file_plot$Probe=="ENSG00000111725",]
file_plot$shape=21
file_plot$shape=as.factor(file_plot$shape)
ggplot(file_plot, aes(BP, -log10(p))) +
  geom_point(aes(shape=shape), color = "black", size = 2,fill="#DDAA33")+theme_bw()+
  #xlim(4500000,6000000)+
  scale_fill_gradient2(low="navy", mid = "lightblue",high="red2",midpoint = 0.4)+
  scale_shape_manual(values = c(21,23))+guides(shape=F)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave("output/PRKAB1.eqtl.zoom.pdf",width = 10,height = 2)

file_plot=file_gwas
file_plot=file_plot[file_plot$SNP %in% file_eqtl$SNP[file_eqtl$Probe=="ENSG00000211445"],]
file_plot=merge(file_plot,file_eqtl[,c("SNP","BP")],by="SNP",all=F)
write.table(file_plot,file = "output/tmp.zoom.file.txt",quote = F,row.names = F)
file_plot$shape=NA
file_plot$shape[file_plot$SNP=="rs744166"]=23
file_plot$shape[file_plot$SNP!="rs744166"]=21
file_plot$shape=as.factor(file_plot$shape)
ggplot(file_plot, aes(BP, -log10(p))) +
  geom_point(aes(shape=shape), color = "black", size = 2,fill="#4477AA")+theme_bw()+
  #xlim(4500000,6000000)+
  scale_shape_manual(values = c(21,23))+guides(shape=F)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave("output/PRKAB1.gwas.zoom.pdf",width = 10,height = 2)

file_plot=file_mbqtl
file_plot=file_plot[file_plot$SNP %in% file_eqtl$SNP[file_eqtl$Probe=="ENSG00000111725"],]
file_plot$shape=21
file_plot$shape=as.factor(file_plot$shape)
file_plot1=file_plot[file_plot$trait=="P164.PWY..purine.nucleobases.degradation.I..anaerobic.",]
ggplot(file_plot1, aes(pos, -log10(pvalue))) +
  geom_point(aes(shape=shape), color = "black", size = 2,fill="#225555")+theme_bw()+
  #xlim(4500000,6000000)+
  scale_shape_manual(values = c(21,23))+guides(shape=F)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave("output/PRKAB1.mbqtl.1.zoom.pdf",width = 10,height = 2)
file_plot2=file_plot[file_plot$trait=="PWY.7237..myo...chiro..and.scillo.inositol.degradation",]
ggplot(file_plot2, aes(pos, -log10(pvalue))) +
  geom_point(aes(shape=shape), color = "black", size = 2,fill="#225555")+theme_bw()+
  #xlim(4500000,6000000)+
  scale_shape_manual(values = c(21,23))+guides(shape=F)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave("output/MUC1.mbqtl.2.zoom.pdf",width = 10,height = 2)

file_plot=file_meqtl
file_plot=file_plot[file_plot$Probe=="cg11195002",]
file_plot=file_plot[file_plot$SNP %in% file_eqtl$SNP[file_eqtl$Probe=="ENSG00000079999"],]
file_plot$shape=21
file_plot$shape=as.factor(file_plot$shape)
file_plot$p=as.numeric(as.character(file_plot$p))
file_plot$BP=as.numeric(as.character(file_plot$BP))
ggplot(file_plot, aes(BP, -log10(p))) +
  geom_point(aes(shape=shape), color = "black", size = 2,fill="#009988")+theme_bw()+
  scale_fill_gradient2(low="navy", mid = "lightblue",high="red2",midpoint = 0.4)+
  scale_shape_manual(values = c(21,23))+guides(shape=F)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave("output/KEAP1.cg11195002.zoom.pdf",width = 10,height = 2)

#=====================
# SMR
#=====================

source("plot/plot_SMR.r")
library(devtools)
file_plot=ReadSMRData("output/PRKAB1.GWAS.eQTL.ENSG00000111725.txt")
file_plot$eQTL$V2=-file_plot$eQTL$V2
pdf("output/test.pdf",width = 5,height = 4)
SMREffectPlot(data=file_plot, trait_name="cg22863889") 
dev.off()
#SMRLocusPlot(data=file_plot, smr_thresh=8.4e-6, heidi_thresh=0.05, plotWindow=1000, max_anno_probe=16)



file_plot1=file_meqtl
file_plot1=file_plot1[file_plot1$Probe=="cg18810026",]
file_plot1=file_plot1[,c("SNP","Chr","BP","A1","A2","b","SE")]
file_plot1$prbname="cg18810026"

file_plot2=file_eqtl
file_plot2=file_plot2[file_plot2$Probe=="ENSG00000168610",]
file_plot2=file_plot2[,c("SNP","Chr","BP","A1","A2","b","SE")]
file_plot2$prbname="ENSG00000168610"

SNP=rbind(file_plot1[,c("SNP","Chr","BP","A1","A2")],file_plot2[,c("SNP","Chr","BP","A1","A2")])
SNP=SNP[!duplicated(SNP$SNP),]

file_plot1=file_plot1[,c("prbname","SNP","b","SE")]
file_plot1$V4=NA
file_plot2=file_plot2[,c("SNP","b","SE")]

tmp.data=list(probeID="cg18810026",SMR="NA","SNP"=SNP,GWAS=file_plot2,eQTL=file_plot1)
SMREffectPlot(data=tmp.data, trait_name="cg18810026") 
ggplot(data, aes(x=type, y= est, color = type)) + 
  theme_bw()+
  geom_point(size = 3.5,position=pd) + 
  geom_errorbar(aes(ymin=l2.5, ymax=est+u97.5), width=0.1,position = pd)+
  geom_point(data=data2, aes(x=type,y=yval),size=1)


#=====================
# coloc analysis
#=====================

library(coloc)
library(snpStats)
library(locuscomparer)
gene_file=read.table("output/ENSG00000101017.Tissue.txt",header = T,stringsAsFactors = F)
pwy_file=read.table("output/ENSG00000101017_GWAS_PYRIDNUCSYN.PWY.txt",header = F,stringsAsFactors = F)
colnames(pwy_file)=c("SNP","p","chr","bp","A1","A2","maf","beta","se","trait")

overlap=intersect(gene_file$SNP,pwy_file$SNP)
tmp.data1=gene_file[gene_file$SNP %in% overlap,]
tmp.data2=pwy_file[pwy_file$SNP %in% overlap,]

tmp.data1=tmp.data1[order(tmp.data1$SNP),]
tmp.data2=tmp.data2[order(tmp.data2$SNP),]

mm=coloc.abf(dataset1=list(pvalues=tmp.data1$p,N=1000,type="quant"),
             dataset2=list(pvalues=tmp.data2$p,N=7749,type="quant"),
             MAF=tmp.data2$maf)

tmp.data1=tmp.data1[,c("SNP","p")]
tmp.data2=tmp.data2[,c("SNP","p")]
colnames(tmp.data1)=c("rsid","pval")
colnames(tmp.data2)=c("rsid","pval")

marker_col="SNP"
pval_col="p"
pdf("output/test.pdf",width = 2.6,height = 2.6)
pp=locuscompare(in_fn1=tmp.data1, in_fn2=tmp.data2, title1="GWAS", title2="eQTL",
             marker_col1= marker_col, pval_col1=pval_col, marker_col2=marker_col, pval_col2=pval_col,combine = F)
pp$locuscompare
dev.off()


