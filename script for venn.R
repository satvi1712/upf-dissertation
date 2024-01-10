library(tidyverse)
library(readr)
library(VennDiagram)
library(ggplot2)
library(gplots)  
library(clusterProfiler)
C_rv_poly<- read_csv("data/ctrl_rv_poly_comparison_transcripts.csv")
C_rv_total<- read_csv("data/ctrl_rv_total_comparison_transcripts.csv")
upf_med_total<- read_csv("data/upf_med_total_comparison_transcripts.csv")
upf_med_poly<- read_csv("data/upf1_med_poly_comparison_transcripts.csv")
upf_rv_total<- read_csv("data/upf1_rv_total_comparison_transcripts.csv")
upf_rv_poly<- read_csv("data/upf1_rv_poly_comparison_transcripts.csv")
#extract genes differentially expressed in control in RV 
ctrl_genes_total<- C_rv_total%>%
  filter(!is.na(padj))%>%
  filter(padj<=0.05)%>%
  select(...1)
ctrl_genes_poly<- C_rv_poly%>%
  filter(!is.na(padj))%>%
  filter(padj<=0.05)%>%
  select(...1)
ctrl_totalvspoly_genes<- inner_join(ctrl_genes_total, ctrl_genes_poly)
#draw venndiagram
grid.newpage()
draw.pairwise.venn(area1= 121+7, area2 = 207+7, cross.area = 7, category = c("total", "poly"), fill="yellow", filename= "ctrl_totalvspoly_venn.png")
#extract genes differentially expressed after upf knockdown in medium 
upf_genes_med_total<- upf_med_total%>%
  filter(!is.na(padj))%>%
  filter(padj<=0.05)%>%
  select(...1)
upf_genes_med_poly<- upf_med_poly%>%
  filter(!is.na(padj))%>%
  filter(padj<=0.05)%>%
  select(...1)
upf_med_totalvspoly<- inner_join(upf_genes_med_total, upf_genes_med_poly)
#draw venndiagram
grid.newpage()
draw.pairwise.venn(area1= 332+14, area2 = 137+14, cross.area = 14, category = c("total", "poly"), fill = "lightblue", filename= "upf_med_totalvspoly.png")
#extract genes differentially expressed after upf knockdown in rv
upf1_genes_rv_total<- upf_rv_total%>%
  filter(!is.na(padj))%>%
  filter(padj<=0.05)%>%
  select(...1)
upf1_genes_rv_poly<- upf_rv_poly%>%
  filter(!is.na(padj))%>%
  filter(padj<=0.05)%>%
  select(...1)
upf1_rv_polyvstotal<- inner_join(upf1_genes_rv_total, upf1_genes_rv_poly)
#draw venndiagram
grid.newpage()
draw.pairwise.venn(area1= 224+2, area2 = 121+2, cross.area = 2, category = c("total", "poly"), fill= "lightpink", filename= "upf_rv_totalvspoly.png")



#import genome data 
genome_data<- read_csv("data/genome_data.csv")
#count biotype of all genes
biotype_counts<-genome_data%>%
  group_by(biotype)%>%
  summarise(n=n())




#extracting NMD genes differentially expressed in poly in upf med 
upf_med_poly_padj<- upf_med_poly%>%
  filter(!is.na(padj))%>%
  filter(padj<=0.05)%>%
  select(transcripts, log2FoldChange, baseMean, pvalue, padj, stat, lfcSE)
 





#extracting genes differentially expressed in total in upfmed total 
upf_med_total_padj<- upf_med_total%>%
  filter(!is.na(padj))%>%
  filter(padj<=0.05)%>%
  select(transcripts, log2FoldChange, baseMean, pvalue, padj, stat, lfcSE)
  
#extracting genes of all biotypes differentially expressed in upf med poly and total from genome_data 
upf_effects_poly_med_allbiotypes<- inner_join(genome_data, upf_med_poly_padj, by= "transcripts")
upf_effects_total_med_allbiotypes<- inner_join(genome_data, upf_med_total_padj, by= "transcripts")
#biotype counts for upf_med data fro both poly and total 
upf_med_poly_biotypecounts<- upf_effects_poly_med_allbiotypes%>%
  group_by(biotype)%>%
  summarise(n=n())
upf_med_total_biotypecounts<- upf_effects_total_med_allbiotypes%>%
  group_by(biotype)%>%
  summarise(n=n())
  
#extracting genes of all biotypes differentially expressed in ctrl rv poly and total from genome data
ctrl_rv_poly_padj<- C_rv_poly%>%
  filter(!is.na(padj))%>%
  filter(padj<=0.05)%>%
  select(transcripts, log2FoldChange, baseMean, pvalue, padj, stat, lfcSE)
  

ctrl_rv_total_padj<- C_rv_total%>%
  filter(!is.na(padj))%>%
  filter(padj<=0.05)%>%
  select(transcripts, log2FoldChange, baseMean, pvalue, padj, stat, lfcSE)
  
rv_effects_ctrl_poly_allbiotypes<- inner_join(genome_data, ctrl_rv_poly_padj, by="transcripts")  
rv_effects_ctrl_total_allbiotypes<- inner_join(genome_data, ctrl_rv_total_padj, by="transcripts")
#biotype counts for rv_crtrl data for total and poly
rv_ctrl_poly_biotypeecounts<- rv_effects_ctrl_poly_allbiotypes%>%
  group_by(biotype)%>%
  summarise(n=n())
rv_ctrl_total_biotypecounts<- rv_effects_ctrl_total_allbiotypes%>%
  group_by(biotype)%>%
  summarise(n=n())

#extracting genes of all biotypes differentially expressed in upf rv poly and total from genome data.
upf_rv_poly_padj<- upf_rv_poly%>%
  filter(!is.na(padj))%>%
  filter(padj<=0.05)%>%
  select(transcripts, log2FoldChange, baseMean, pvalue, padj, stat, lfcSE)
  
upf_rv_total_padj<- upf_rv_total%>%
  filter(!is.na(padj))%>%
  filter(padj<=0.05)%>%
  select(transcripts, log2FoldChange, baseMean, pvalue, padj, stat, lfcSE)
  
upf_rv_effects_poly_allbiotypes<- inner_join(genome_data, upf_rv_poly_padj, by="transcripts")
upf_rv_effects_total_allbiotypes<- inner_join(genome_data, upf_rv_total_padj, by="transcripts")
#biotype counts for upf rv extracted data for both total and poly data
upf_rv_poly_biotypecounts<- upf_rv_effects_poly_allbiotypes%>%
  group_by(biotype)%>%
  summarise(n=n())
upf_rv_total_biotypecounts<- upf_rv_effects_total_allbiotypes%>%
  group_by(biotype)%>%
  summarise(n=n())

#deg and dbg modulated in both upffrv and ctrlrv in poly and total
upfrvctrlrvgenes_total<- inner_join(upf_rv_total_padj, ctrl_rv_total_padj, by = "transcripts")
upfrvctrlrvgenes_poly<- inner_join(upf_rv_poly_padj, ctrl_rv_poly_padj, by = "transcripts")
#deg and dbg modulated in upfrv and upf med
upfmedupfrvgenes_total<- inner_join(upf_med_total_padj, upf_rv_total_padj, by="transcripts")
upfmedupfrvgenes_poly<- inner_join(upf_med_poly_padj, upf_rv_poly_padj, by="transcripts")


#NMD candidates directionality change- postive or negative (log2FC)
#upf med poly and total
positiveFC_NMD_upf_med_poly<- NMD_upf_med_genes_poly%>%
  filter(log2FoldChange>0)
negativeFC_NMD_upf_med_poly<- NMD_upf_med_genes_poly%>%
  filter(log2FoldChange<0)
positiveFC_NMD_upf_med_total<- NMD_upf_med_genes_total%>%
  filter(log2FoldChange>0)
negativeFC_NMD_upf_med_total<- NMD_upf_med_genes_total%>%
  filter(log2FoldChange<0)

#ctrl rv poly vs total 
positiveFC_NMD_ctrl_rv_poly<- NMD_ctrl_rv_genes_poly%>%
  filter(log2FoldChange>0)
negativeFC_NMD_ctrl_rv_poly<- NMD_ctrl_rv_genes_poly%>%
  filter(log2FoldChange<0)
positiveFC_NMD_ctrl_rv_total<- NMD_ctrl_rv_genes_total%>%
  filter(log2FoldChange>0)
negativeFC_NMD_ctrl_rv_total<-NMD_ctrl_rv_genes_total%>%
  filter(log2FoldChange<0)

#upf rv poly vs total
positiveFC_NMD_upf_rv_poly<- NMD_upf_rv_genes_poly %>%
  filter(log2FoldChange>0)
negativeFC_NMD_upf_rv_poly<- NMD_upf_rv_genes_poly%>%
  filter(log2FoldChange<0)
postiveFC_NMD_upf_rv_total<- NMD_upf_rv_genes_total%>%
  filter(log2FoldChange>0)
negativeFC_NMD_upf_rv_total<- NMD_upf_rv_genes_total%>%
  filter(log2FoldChange<0)

#export reactome data for upd med, rv ctrl and upf rv- total and poly- used data from allbiotype data 
upf_med_poly_pathwaydata<- read_csv("data/upf_med_pathwaymapped_poly.csv")
upf_med_total_pathwaydata<- read_csv("data/upf_med_pathwaymapped_total.csv")
rv_ctrl_poly_pathwaydata<- read_csv("data/rv_ctrl_pathwaymapped_poly.csv")
rv_ctrl_total_pathwaydata<- read_csv("data/rv_ctrl_pathwaymapped_total.csv")
upf_rv_total_pathwaydata<- read_csv("data/upf_rv_total_pathwaymapped.csv")
upf_rv_poly_pathwaydata<-read_csv("data/upf_rv_poly_pathwaymapped.csv")

#extracting entities fdr<0.05 for all the data sets 
upf_med_poly_fdr<- upf_med_poly_pathwaydata%>%
  filter(`Entities_FDR`<=0.05)%>%
  dplyr::select(`Pathway_name`, `Entities_found`, `Entities_FDR`) 
upf_med_total_fdr<- upf_med_total_pathwaydata%>%
  filter(`Entities_FDR`<=0.05)%>%
  dplyr::select(`Pathway_name`, `Entities_found`, `Entities_FDR`)
rv_ctrl_poly_fdr<- rv_ctrl_poly_pathwaydata%>%
  filter(`Entities_FDR`<=0.05)%>%
  dplyr::select(`Pathway_name`, `Entities_found`, `Entities_FDR`)
rv_ctrl_total_fdr<- rv_ctrl_total_pathwaydata%>%
  filter(`Entities_FDR`<=0.05)%>%
  dplyr::select(`Pathway_name`, `Entities_found`, `Entities_FDR`)
rv_upf_poly_fdr<- upf_rv_poly_pathwaydata%>%
  filter(`Entities_FDR`<=0.05)%>%
  dplyr::select(`Pathway_name`, `Entities_found`, `Entities_FDR`)
rv_upf_total_fdr<- upf_rv_total_pathwaydata%>%
  filter(`Entities_FDR`<=0.05)%>%
  dplyr::select(`Pathway_name`, `Entities_found`, `Entities_FDR`)



#extracting NMD genes and non NMD from upf rv
NMD_upf_rv_genes_poly<- upf_rv_effects_poly_allbiotypes%>%
  filter(biotype=="nonsense_mediated_decay")
NMD_upf_rv_genes_total<- upf_rv_effects_total_allbiotypes%>%
  filter(biotype=="nonsense_mediated_decay")

nonNMD_upf_rv_genes_total<- upf_rv_effects_total_allbiotypes%>%
  filter(!biotype=="nonsense_mediated_decay")
nonNMD_upf_rv_genes_poly<- upf_rv_effects_poly_allbiotypes%>%
  filter(!biotype=="nonsense_mediated_decay")

#extracting NMD and non-NMD genes from ctrl rv
NMD_ctrl_rv_genes_poly<- rv_effects_ctrl_poly_allbiotypes%>%
  filter(biotype=="nonsense_mediated_decay")
NMD_ctrl_rv_genes_total<- rv_effects_ctrl_total_allbiotypes%>%
  filter(biotype=="nonsense_mediated_decay")

nonNMD_ctrl_rv_genes_poly<- rv_effects_ctrl_poly_allbiotypes%>%
  filter(!biotype=="nonsense_mediated_decay")
nonNMD_ctrl_rv_genes_total<- rv_effects_ctrl_total_allbiotypes%>%
  filter(!biotype=="nonsense_mediated_decay")

#extracting NMD and non-NMD genes from upf med
NMD_upf_med_genes_poly<- upf_effects_poly_med_allbiotypes%>%
  filter(biotype=="nonsense_mediated_decay")
NMD_upf_med_genes_total<- upf_effects_total_med_allbiotypes%>%
  filter(biotype=="nonsense_mediated_decay")
nonNMD_upf_med_genes_poly<- upf_effects_poly_med_allbiotypes%>%
  filter(!biotype=="nonsense_mediated_decay")
nonNMD_upf_med_genes_total<- upf_effects_total_med_allbiotypes%>%
  filter(!biotype=="nonsense_mediated_decay")

#exporting data for sever asthmatics vs healthy controls- supplental data 5 and 7
savshc_poly_dei<- read_csv("data/poly_isoforms_p0.05.csv")
savshc_total_dei<- read_csv("data/total_isoforms_p0.05.csv")

#merging data with upf med data and savshc
DEI_upfmed_svhc_poly<- inner_join(savshc_poly_dei, upf_effects_poly_med_allbiotypes, by="cdinggene")
DEI_upfmed_svhc_total<- inner_join(savshc_total_dei, upf_effects_total_med_allbiotypes, by="cdinggene")

#make dei poly dataframe 
DEI_poly<- data.frame(
  Isoform = c("uc003jvs.4","uc003jhg.2","uc021pix.1","uc003aui.3","uc001gkj.1","uc002zdc.1","uc004cgv.4","uc021ult.1","uc002zoq.1","uc002nkf.3","uc003dwv.1","uc002grz.4","uc001ojk.1","uc001yqh.3","uc003pvd.3","uc002onn.3","uc001oht.3", "uc003uew.3", "uc003xlc.3") ,
  siUPF1 = c(-9.058685844,-8.302756795,-8.049344339,-1.083580294,6.923971601,4.479729742,9.044828988,9.024802653,-7.992347707,-9.395251194,-0.99239979,9.472065182,1.396008736,-6.871351715,8.937087915,8.844278368,-4.869702417,-0.965516007,-8.271628198), 
  SA = c(0.657,4.611,0.674,-0.605,0.672,2.053,-0.572,-0.588,-0.586,-0.536,0.51,-0.684,-0.488,-0.774,3.401,-0.719,-3.938,-0.726,2.639)
)
rownames(	DEI_poly) <- DEI_poly$Isoform 
DEI_poly <- DEI_poly[, -1]

DEI_poly_datamatrix<- data.matrix(DEI_poly)# convert into data matrix
#colour pallete fro heatmap 
cs<- colorRampPalette(c("blue", "white", "red"))
#create heatmap using log2fold change and log ratio 
heatmap.2(DEI_poly_datamatrix, main = "DEI in poly",col = cs, trace = "none", margins= c(10,12))
#make dei total dataframe 
DEI_total<- data.frame(
  Isoform = c("uc003ptt.2","uc003tns.3","uc002dpo.3","uc003tnq.3","uc011lbv.2","uc021paz.1","uc001rqe.3","uc010fzo.3","uc010unh.1","uc002igm.2","uc002exm.2","uc004enz.1","uc001jtb.2","uc002vyw.4","uc010znp.2","uc021wcr.1","uc001agm.3","uc001mjp.3","uc010uxs.2","uc002njd.2","uc003zax.3","uc011mkm.2","uc001hvg.3","uc003xsk.4","uc003aui.3","uc003zpk.3"),
  siUPF1 =c (8.857637,1.974022253,-7.662387,1.974022253,-5.7323508,-6.32743,-1.30513,-8.957886,-8.122807,-0.7334933,3.707793,-7.785235105,-8.357679876,-0.5997059,8.140293,-8.145928206,-2.529800072,-8.458302113,8.440247426,-1.453329,-0.4600329,-8.037110648,8.549458,-9.502695932,-0.860556151,-0.338120637),
  SA =c(-0.295,-4.086,-0.507,-1.405,-3.56,-0.327,0.651,3.345,-2.187,0.225,3.843,-0.493,-0.193,-0.377,2.53,-0.491,-3.061,1.858,-2.814,-2.688,-0.246,-2.981,-0.655,0.435,-0.329,1.063)										
)
rownames(DEI_total) <- DEI_total$Isoform
DEI_total<- DEI_total[,-1]

DEI_total_datamatrix<- data.matrix(DEI_total)
#create heatmap
heatmap.2(DEI_total_datamatrix, main = "DEI in total",col = cs, trace = "none", margins = c(10,12), cexRow = 0.5)


#genomewide NMD and non-NMD proportions
NMD_transcripts_genomedata<- genome_data%>%
  filter(biotype=="nonsense_mediated_decay")
nonNMD_trancripts_genomedata<- genome_data%>%
  filter(!biotype=="nonsense_mediated_decay")

#finding downregulated and upregulated genes in ctrl rv 
upregulated_ctrl_rv_total<- rv_effects_ctrl_total_allbiotypes%>%
  filter(log2FoldChange>0)%>%
  select(cdinggene)
downregulated_ctrl_rv_total<- rv_effects_ctrl_total_allbiotypes%>%
  filter(log2FoldChange<0)%>%
  select(cdinggene)
upregulatedctrlrvpoly<- rv_effects_ctrl_poly_allbiotypes%>%
  filter(log2FoldChange>0)%>%
  select(cdinggene)
downregulatedctrlrvpoly<- rv_effects_ctrl_poly_allbiotypes%>%
  filter(log2FoldChange<0)%>%
  select(cdinggene)
  



#finding downregulated and upregulated genes in upf rv
downregulated_total_upfrv<- upf_rv_effects_total_allbiotypes%>%
  filter(log2FoldChange<0)
 
upregulated_poly_upfrv<- upf_rv_effects_poly_allbiotypes%>%
  filter(log2FoldChange>0)
  

upregulated_total_upfrv<- upf_rv_effects_total_allbiotypes%>%
  filter(log2FoldChange>0)
  
downregulated_poly_upfrv<- upf_rv_effects_poly_allbiotypes%>%
  filter(log2FoldChange<0)
 


#finding downregulated and upregulated genes in upf med 
downregulated_total<- upf_effects_total_med_allbiotypes%>%
  filter(log2FoldChange<0)
upregulated_poly<- upf_effects_poly_med_allbiotypes%>%
  filter(log2FoldChange>0)
upregulated_total<- upf_effects_total_med_allbiotypes%>%
  filter(log2FoldChange>0)
downregulated_poly<- upf_effects_poly_med_allbiotypes%>%
  filter(log2FoldChange<0)

# exporting upregulated and downregulated upf med data into excel 
write.csv(upregulated_poly, file = "upfmedpolyupregulated.csv")
write.csv(upregulated_total, file = "upfmedtotalupregulated.csv")
write.csv(downregulated_poly, file = "upfmedpolydownregulated.csv")
write.csv(downregulated_total, file = "upfmedtotaldownregulated.csv")

#importing upregulation and downregulation upf med pathway enrichment data
upfmed_downregulated_poly_reactome<- read_csv("data/upfmed_poly_downregulated_reactome.csv")
upfmed_downregulated_total_reactome<- read_csv("data/upfmed_total_downregulated_reactome.csv")
upfmed_upregulated_poly_reactome<- read_csv("data/upfmed_upregulated_reactome_poly.csv")
upfmed_upregulated_total_reactome<- read_csv("data/upfmed_total_upregulated_reactome.csv")

#filter FDR<= 0.05
upfmedpoly_dr_fdr<- upfmed_downregulated_poly_reactome%>%
filter(`Entities_FDR`<=0.05)%>%
  dplyr::select(`Pathway_name`, `Entities_found`, `Entities_FDR`) 
upfmedtotal_dr_fdr<- upfmed_downregulated_total_reactome%>%
filter(`Entities_FDR`<=0.05)%>%
  dplyr::select(`Pathway_name`, `Entities_found`, `Entities_FDR`) 
upfmedpoly_ur_fdr<- upfmed_upregulated_poly_reactome%>%
filter(`Entities_FDR`<=0.05)%>%
  dplyr::select(`Pathway_name`, `Entities_found`, `Entities_FDR`) 
upfmedtotal_ur_fdr<-upfmed_upregulated_total_reactome%>%
filter(`Entities_FDR`<=0.05)%>%
  dplyr::select(`Pathway_name`, `Entities_found`, `Entities_FDR`) 

#make plots
#poly fractions had no significant changes(FDR<=0.05) in deg after upf knockdown 
#total fraction
ggplot(data = upfmedtotal_dr_fdr, aes(x = Entities_found, y = Pathway_name, color= Entities_FDR))+ geom_bar(stat = 'identity')
plot_ur_upfmedtotal<-ggplot(data = upfmedtotal_ur_fdr, aes(x = Entities_found, y = Pathway_name, color= Entities_FDR))+ geom_bar(stat = 'identity')
png("upfmedtotal_ur.png", res=600, width = 6000, height = 6000)
print(plot_ur_upfmedtotal)
dev.off()

# volcano plot datasets 
#filtering padj=NA
#ctrl rv 
ctrlrv_total_NA<- C_rv_total%>%
  filter(!is.na(padj))
  
ctrlrv_poly_NA<- C_rv_poly%>%
  filter(!is.na(padj))
  
#upf med
upfmed_total_NA<- upf_med_total%>%
  filter(!is.na(padj))
 
upfmed_poly_NA<- upf_med_poly%>%
  filter(!is.na(padj))
  
#upf rv 
upfrv_total_NA<- upf_rv_total%>%
  filter(!is.na(padj))
  
upfrv_poly_NA<- upf_rv_poly%>%
  filter(!is.na(padj))
  
#merging datsets with genome data to find gene symbols
#ctrl rv
ctrlrvtotal_NA_geneid<- inner_join(ctrlrv_total_NA, genome_data, by = "transcripts")
ctrlrvpoly_NA_geneid<- inner_join(ctrlrv_poly_NA, genome_data, by= "transcripts")
#upf med
upfmedtotal_NA_geneid<- inner_join(upfmed_total_NA, genome_data, by= "transcripts")
upfmedpoly_NA_geneid<- inner_join(upfmed_poly_NA, genome_data, by= "transcripts")
#upf rv
upfrvtotal_NA_geneid<- inner_join(upfrv_total_NA, genome_data, by= "transcripts")
upfrvpoly_NA_geneid<- inner_join(upfrv_poly_NA, genome_data, by= "transcripts")

#exporting to excel to find -log10(padj)
write.csv(ctrlrvtotal_NA_geneid, file = "ctrlrvtotalNAgeneid.csv")
write.csv(ctrlrvpoly_NA_geneid, file = "ctrlrvpolyNAgeneid.csv")
write.csv(upfmedtotal_NA_geneid, file = "upfmedtotalNAgeneid.csv")
write.csv(upfmedpoly_NA_geneid, file = "upfmedpolyNAgeneid.csv")
write.csv(upfrvtotal_NA_geneid, file = "upfrvtotalNAgeneid.csv")
write.csv(upfrvpoly_NA_geneid, file = "upfrvpolyNAgeneid.csv")




#ctrl rv biotype counts
ctrlrv_biotypecounts_poly<-rv_effects_ctrl_poly_allbiotypes%>%
  group_by(biotype)%>%
  summarise(n=n())
ctrlrv_biotypecounts_total<-rv_effects_ctrl_total_allbiotypes%>%
  group_by(biotype)%>%
  summarise(n=n())


