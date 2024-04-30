library(DMRcate)
library(lme4)
library(lmerTest)
library(data.table)
library(DescTools)
library(ggplot2)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)

#######################Genomic Feature Enrichment Analysis for 850K Adult SES CpGs#############################
#Reading in Significant CpGs from Adult 850K Biomarker and annotation files 
aSESEPIC_sig=read.csv('/net/orion/home/opsasnic/HRS Biomarker/aSESEPIC.csv', header=T)
genomicfeature=read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/genomeenrich.csv', header=T)
sig_anno=merge(aSESEPIC_sig,genomicfeature,by="cpg", all.x=T)

#Reading in all CpGs on 850K EPIC chip
EPIC=read.csv('/net/orion/home/opsasnic/HRS Biomarker/epic.csv', header=T)
#Cross Reactive Probes to be removed
crossreactive = read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/QC_processing_Smith/EPICcrossreactiveprobelist.csv', header=T)
merge=merge(EPIC,crossreactive,by.x="cpg", by.y="ProbeIDs", all.x=T)
merge2 <- subset(merge, is.na(merge$remove))
merge3 = subset(merge2, select = -c(remove,match47,match48,match49,match50,Total))
all_anno=merge(merge3,genomicfeature,by="cpg", all.x=T)

#Tabulating the significant and total CpGs that falls within each genomic region
table(sig_anno$enhancer)
table(all_anno$enhancer)

table(sig_anno$promoter)
table(all_anno$promoter)

table(sig_anno$shelfshore)
table(all_anno$shelfshore)

table(sig_anno$cpgisland)
table(all_anno$cpgisland)

table(sig_anno$DNAse_site)
table(all_anno$DNAse_site)

#Functional Gene Enrichment Analysis - need to remove the significant CpGs from the total CpGs to get final numbers 
dat <- data.frame(
  "promoter" = c(7,184402),
  "notpromoter" = c(18, 605043),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("Promoter", "Not promoter")
fisher.test(dat)

dat2 <- data.frame(
  "island" = c(8,153466),
  "not island" = c(17, 635979),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat2) <- c("island", "Not island")
fisher.test(dat2)

dat4 <- data.frame(
  "shelf" = c(8,201062),
  "not" = c(17, 588383),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat4) <- c("Shelf", "Not shelf")
fisher.test(dat4)


dat5 <- data.frame(
  "DNAse" = c(22,470194),
  "not DNAse" = c(3, 319251),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat5) <- c("DNAse", "Not DNAse")
fisher.test(dat5)


dat6 <- data.frame(
  "enhancer" = c(3,26037),
  "notenhancer" = c(22, 763408),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat6) <- c("Enhancer", "Not enhancer")
fisher.test(dat6)


#EQTM enrichment analysis 
#Reading in eQTM data file 
eqtm=read.csv("/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/eqtm/eqtm_genes.csv", header=T)
mergedsig=merge(eqtm,aSESEPIC_sig, all.y=T) #25 sign cpgs
noeqtm_sig=mergedsig[is.na(mergedsig$gene),] #17 biomarker CpGs have no transcript pairs 
eqtm_sig=mergedsig[!is.na(mergedsig$gene),] 
eqtm_sig2=eqtm_sig[!duplicated(eqtm_sig$cpg), ] # 8 biomarker CpGs have  transcript pairs

mergedall=merge(eqtm,all_anno, all.y=T) 
noeqtm_all=mergedall[is.na(mergedall$gene),] #752799 non-biomarker CpGs have no transcript pair (excluding significant CpGs)
eqtm_all=mergedall[!is.na(mergedall$gene),] 
eqtm_all2=eqtm_all[!duplicated(eqtm_all$cpg),] #36646 non-biomarker CpGs have transcript pair (excluding significant CpGs)

dat7 <- data.frame(
  "eqtm" = c(8,36646),
  "noteqtm" = c(17, 752799),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat7) <- c("eQTM", "not eQTM")
fisher.test(dat7)

#GO and KEGG ENRICHMENT Analysis Adult SES Biomarker 850K using GOmeth method 
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
sigcpgs=aSESEPIC_sig[,1]
allcpgs <- all_anno[,1]

gst <- gometh(sig.cpg = sigcpgs, all.cpg = allcpgs, 
              collection = "GO", array.type="EPIC")

kegg <- gometh(sig.cpg = sigcpgs, all.cpg = allcpgs, collection = "KEGG", array.type="EPIC", prior.prob=TRUE, anno=ann)

#GO and KEGG ENRICHMENT Analysis Adult SES Biomarker 850K using ClusterProfiler method with eQTM data
#Deduplicating all datasets
background_final=mergedall[!duplicated(mergedall$gene), ]
mergedsig_final=mergedsig[!duplicated(mergedsig$gene), ]

universe=background_final$gene
sig_id=mergedsig_final$gene

eg = bitr(universe, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
egsig = bitr(sig_id, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

eg2=eg$ENTREZID
egsig2=egsig$ENTREZID

ego <- enrichGO(gene          = sig_id,
                universe      = universe,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = .01,
                qvalueCutoff = .01,
                keyType = "SYMBOL")
summary_sig=data.frame(ego)

s_ego<-clusterProfiler::simplify(ego)
sum_GO=as.data.frame(s_ego)
dotplot(s_ego, font.size = 8)

kk <- enrichKEGG(gene         = egsig2,
                 universe =eg2,
                 organism     = 'hsa',
                 keyType = "kegg",
                 pvalueCutoff = .05,
                 qvalueCutoff=.05)
summary_keg=data.frame(kk)


################################Genomic Feature Enrichment Analysis for 450K Adult SES CpGs####################
#Reading in Significant CpGs from Adult 450K Biomarker and annotation files 
aSES450k_sig=read.csv('/net/orion/home/opsasnic/HRS Biomarker/aSES450K.csv', header=T)
genomicfeature=read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/genomeenrich.csv', header=T)
sig_anno=merge(aSES450k_sig,genomicfeature,by="cpg", all.x=T)

#Reading in all CpGs on 850K EPIC chip
EPIC=read.csv('/net/orion/home/opsasnic/HRS Biomarker/epic.csv', header=T)
four50k=read.csv('/net/orion/skardia_lab/clubhouse/research/projects/MESA/Methylation_081413/ValidProbe_DMRcate_MESA.csv', header=T)
cpgs=merge(EPIC,four50k,by.x="cpg",by.y="x")

#Removing cross reactive probes 
crossreactive = read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/QC_processing_Smith/EPICcrossreactiveprobelist.csv', header=T)
merge=merge(cpgs,crossreactive,by.x="cpg", by.y="ProbeIDs", all.x=T)
merge2 <- subset(merge, is.na(merge$remove))
merge3 = subset(merge2, select = -c(remove,match47,match48,match49,match50,Total))
all_anno=merge(merge3,genomicfeature,by="cpg", all.x=T)

#Tabulating the significant and total CpGs that falls within each genomic region
table(sig_anno$enhancer)
table(all_anno$enhancer)

table(sig_anno$promoter)
table(all_anno$promoter)

table(sig_anno$shelfshore)
table(all_anno$shelfshore)

table(sig_anno$cpgisland)
table(all_anno$cpgisland)

table(sig_anno$DNAse_site)
table(all_anno$DNAse_site)

#Functional Gene Enrichment Analysis - need to remove the significant CpGs from the total CpGs to get final numbers 
dat <- data.frame(
  "promoter" = c(10,115222),
  "notpromoter" = c(19, 263617),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("Promoter", "Not promoter")
fisher.test(dat)

dat2 <- data.frame(
  "island" = c(11,122339),
  "not island" = c(18, 256500),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat2) <- c("island", "Not island")
fisher.test(dat2)

dat4 <- data.frame(
  "shelf" = c(11,126066),
  "not" = c(18, 252773),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat4) <- c("Shelf", "Not shelf")
fisher.test(dat4)


dat5 <- data.frame(
  "DNAse" = c(23,242493),
  "not DNAse" = c(6, 136346),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat5) <- c("DNAse", "Not DNAse")
fisher.test(dat5)

dat6 <- data.frame(
  "enhancer" = c(2,6461),
  "notenhancer" = c(27, 372378),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat6) <- c("Enhancer", "Not enhancer")
fisher.test(dat6)

#EQTM enrichment analysis 
eqtm=read.csv("/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/eqtm.csv", header=T)
sep=separate(eqtm, col=GeneSymbol, into=c('gene1', 'gene2', 'gene3', 'gene4', 'gene5', 'gene6','gene7', 'gene8', 'gene9', 'gene10', 'gene11', 'gene12', 'gene13', 'gene14'), sep='\\|')
gene1=sep[,c(1,2)]
gene2=sep[,c(1,3)]
gene3=sep[,c(1,4)]
gene4=sep[,c(1,5)]
gene5=sep[,c(1,6)]
gene6=sep[,c(1,7)]
gene7=sep[,c(1,8)]
gene8=sep[,c(1,9)]
gene9=sep[,c(1,10)]
gene10=sep[,c(1,11)]
gene11=sep[,c(1,12)]
gene12=sep[,c(1,13)]
gene13=sep[,c(1,14)]
gene14=sep[,c(1,15)]


names(gene1)[names(gene1) == "gene1"] <- "gene"
names(gene2)[names(gene2) == "gene2"] <- "gene"
names(gene3)[names(gene3) == "gene3"] <- "gene"
names(gene4)[names(gene4) == "gene4"] <- "gene"
names(gene5)[names(gene5) == "gene5"] <- "gene"
names(gene6)[names(gene6) == "gene6"] <- "gene"
names(gene7)[names(gene7) == "gene7"] <- "gene"
names(gene8)[names(gene8) == "gene8"] <- "gene"
names(gene9)[names(gene9) == "gene9"] <- "gene"
names(gene10)[names(gene10) == "gene10"] <- "gene"
names(gene11)[names(gene11) == "gene11"] <- "gene"
names(gene12)[names(gene12) == "gene12"] <- "gene"
names(gene13)[names(gene13) == "gene13"] <- "gene"
names(gene14)[names(gene14) == "gene14"] <- "gene"

allgenes=rbind(gene1,gene2,gene3,gene4,gene5,gene6,gene7,gene8,gene9,gene10,gene11,gene12,gene13,gene14)
allgenes2=allgenes[!is.na(allgenes$gene),] 
mergedsig=merge(allgenes2,aSES450k_sig,all.y=T)

noeqtm_sig=mergedsig[is.na(mergedsig$gene),] 
noeqtm_sig2=noeqtm_sig[!duplicated(noeqtm_sig$cpg), ] #11 biomarker CpGs have no transcript pairs 

eqtm_sig=mergedsig[!is.na(mergedsig$gene),] 
eqtm_sig2=eqtm_sig[!duplicated(eqtm_sig$cpg), ] #18 biomarker CpGs have transcript pairs 

mergedall=merge(allgenes2,all_anno, all.y=T) 
noeqtm_all=mergedall[is.na(mergedall$gene),] #352290 non-biomarker CpGs have no transcript pairs (excluding significant CpGs)
eqtm_all=mergedall[!is.na(mergedall$gene),] 
eqtm_all2=eqtm_all[!duplicated(eqtm_all$cpg),] #26549 non-biomarker CpGs have transcript pairs (excluding significant CpGs)

dat7 <- data.frame(
  "eqtm" = c(18,26549),
  "noteqtm" = c(11, 352290),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat7) <- c("eQTM", "not eQTM")
fisher.test(dat7)

#GO and KEGG ENRICHMENT Analysis Adult SES Biomarker 450K using GOmeth method 
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
sigcpgs=aSES450k_sig[,1]
allcpgs <- all_anno[,1]

gst <- gometh(sig.cpg = sigcpgs, all.cpg = allcpgs, 
              collection = "GO", array.type="450K", anno = ann)

kegg <- gometh(sig.cpg = sigcpgs, all.cpg = allcpgs, collection = "KEGG", array.type="450K", prior.prob=TRUE, anno=ann)


#GO and KEGG ENRICHMENT Analysis Adult SES Biomarker 450K using ClusterProfiler method with eQTM data
#Deduplicating all datasets
background_final=mergedall[!duplicated(mergedall$gene), ]
mergedsig_final=mergedsig[!duplicated(mergedsig$gene), ]
mergedsig_final=na.omit(mergedsig_final)

universe=background_final$gene
sig_id=mergedsig_final$gene

eg = bitr(universe, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
egsig = bitr(sig_id, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

eg2=eg$ENTREZID
egsig2=egsig$ENTREZID

ego <- enrichGO(gene          = sig_id,
                universe      = universe,
                OrgDb         = org.Hs.eg.db,
                ont           = "all",
                pAdjustMethod = "BH",
                pvalueCutoff = .01,
                qvalueCutoff = .01,
                keyType = "SYMBOL")
summary_sig=data.frame(ego)
s_ego<-clusterProfiler::simplify(ego)
sum_GO=as.data.frame(s_ego)
dotplot(s_ego)

kk <- enrichKEGG(gene         = egsig2,
                 universe =eg2,
                 organism     = 'hsa',
                 keyType = "kegg",
                 pvalueCutoff = 0.01,
                 qvalueCutoff=0.01)
summary_keg=data.frame(kk)

#######################Genomic Feature Enrichment Analysis for 850K Childhood SES CpGs#############################
#Reading in Significant CpGs from Childhood 850K Biomarker and annotation files 
cSESEPIC_sig=read.csv('/net/orion/home/opsasnic/HRS Biomarker/cSESEPIC.csv', header=T)
genomicfeature=read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/genomeenrich.csv', header=T)
sig_anno=merge(cSESEPIC_sig,genomicfeature,by="cpg", all.x=T)

#Reading in all CpGs on 850K EPIC chip
EPIC=read.csv('/net/orion/home/opsasnic/HRS Biomarker/epic.csv', header=T)
#Cross Reactive Probes to be removed
crossreactive = read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/QC_processing_Smith/EPICcrossreactiveprobelist.csv', header=T)
merge=merge(EPIC,crossreactive,by.x="cpg", by.y="ProbeIDs", all.x=T)
merge2 <- subset(merge, is.na(merge$remove))
merge3 = subset(merge2, select = -c(remove,match47,match48,match49,match50,Total))
all_anno=merge(merge3,genomicfeature,by="cpg", all.x=T)

#Tabulating the significant and total CpGs that falls within each genomic region
table(sig_anno$enhancer)
table(all_anno$enhancer)

table(sig_anno$promoter)
table(all_anno$promoter)

table(sig_anno$shelfshore)
table(all_anno$shelfshore)

table(sig_anno$cpgisland)
table(all_anno$cpgisland)

table(sig_anno$DNAse_site)
table(all_anno$DNAse_site)


#Functional Gene Enrichment Analysis - need to remove the significant CpGs from the total CpGs to get final numbers 
dat <- data.frame(
  "promoter" = c(8,184401),
  "notpromoter" = c(8, 605053),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("Promoter", "Not promoter")
fisher.test(dat)

dat2 <- data.frame(
  "island" = c(6,153468),
  "not island" = c(10, 635986),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat2) <- c("island", "Not island")
fisher.test(dat2)


dat4 <- data.frame(
  "shelf" = c(6,201064),
  "not" = c(10, 588390),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat4) <- c("Shelf", "Not shelf")
fisher.test(dat4)

dat5 <- data.frame(
  "DNAse" = c(12,470204),
  "not DNAse" = c(4, 319250),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat5) <- c("DNAse", "Not DNAse")
fisher.test(dat5)

dat6 <- data.frame(
  "enhancer" = c(2,26038),
  "notenhancer" = c(14, 763416),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat6) <- c("Enhancer", "Not enhancer")
fisher.test(dat6)

#EQTM enrichment
eqtm=read.csv("/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/eqtm/eqtm_genes.csv", header=T)
mergedsig=merge(eqtm,cSESEPIC_sig, all.y=T) #17 sign cpgs
noeqtm_sig=mergedsig[is.na(mergedsig$gene),] #12 biomarker CpGs have no transcript pairs
eqtm_sig=mergedsig[!is.na(mergedsig$gene),] 
eqtm_sig2=eqtm_sig[!duplicated(eqtm_sig$cpg), ] # 5 biomarker CpGs have transcript pairs

mergedall=merge(eqtm,all_anno, all.y=T) 
noeqtm_all=mergedall[is.na(mergedall$gene),] #75804 non-biomarker CpGs have no transcript pairs (excluding significant CpGs)
eqtm_all=mergedall[!is.na(mergedall$gene),] 
eqtm_all2=eqtm_all[!duplicated(eqtm_all$cpg),] #36649 non-biomarker CpGs have transcript pairs (excluding significant CpGs)


dat7 <- data.frame(
  "eqtm" = c(5,36649),
  "noteqtm" = c(12, 752804),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat7) <- c("eQTM", "not eQTM")
fisher.test(dat7)

#GO and KEGG ENRICHMENT Analysis Childhood SES Biomarker 850K using GOmeth method 
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
sigcpgs=cSESEPIC_sig[,1]
allcpgs <- all_anno[,1]

gst <- gometh(sig.cpg = sigcpgs, all.cpg = allcpgs, 
              collection = "GO", array.type="EPIC", anno = ann)

kegg <- gometh(sig.cpg = sigcpgs, all.cpg = allcpgs, collection = "KEGG", array.type="EPIC", prior.prob=TRUE, anno=ann)



#GO and KEGG ENRICHMENT Analysis Childhood SES Biomarker 850K using ClusterProfiler method with eQTM data
#Deduplicating all datasets
background_final=mergedall[!duplicated(mergedall$gene), ]
mergedsig_final=mergedsig[!duplicated(mergedsig$gene), ]

universe=background_final$gene
sig_id=mergedsig_final$gene

eg = bitr(universe, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
egsig = bitr(sig_id, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

eg2=eg$ENTREZID
egsig2=egsig$ENTREZID


ego <- enrichGO(gene          = sig_id,
                universe      = universe,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = .01,
                qvalueCutoff = .01,
                keyType = "SYMBOL")
summary_sig=data.frame(ego)

s_ego<-clusterProfiler::simplify(ego)
sum_GO=as.data.frame(s_ego)
dotplot(s_ego, font.size = 8)

kk <- enrichKEGG(gene         = egsig2,
                 universe =eg2,
                 organism     = 'hsa',
                 keyType = "kegg",
                 pvalueCutoff = .01,
                 qvalueCutoff=.01)
summary_keg=data.frame(kk)


################################Genomic Feature Enrichment Analysis for 450K Childhood SES CpGs####################
#Reading in Significant CpGs from Childhood 450K Biomarker and annotation files 
cSES450k_sig=read.csv('/net/orion/home/opsasnic/HRS Biomarker/cSES450K.csv', header=T)
genomicfeature=read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/genomeenrich.csv', header=T)
sig_anno=merge(cSES450k_sig,genomicfeature,by="cpg", all.x=T)

#Reading in all CpGs on 850K EPIC chip
EPIC=read.csv('/net/orion/home/opsasnic/HRS Biomarker/epic.csv', header=T)
four50k=read.csv('/net/orion/skardia_lab/clubhouse/research/projects/MESA/Methylation_081413/ValidProbe_DMRcate_MESA.csv', header=T)
cpgs=merge(EPIC,four50k,by.x="cpg",by.y="x")

#Removing cross reactive probes 
crossreactive = read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/QC_processing_Smith/EPICcrossreactiveprobelist.csv', header=T)
merge=merge(cpgs,crossreactive,by.x="cpg", by.y="ProbeIDs", all.x=T)
merge2 <- subset(merge, is.na(merge$remove))
merge3 = subset(merge2, select = -c(remove,match47,match48,match49,match50,Total))
all_anno=merge(merge3,genomicfeature,by="cpg", all.x=T)

#Tabulating the significant and total CpGs that falls within each genomic region
table(sig_anno$enhancer)
table(all_anno$enhancer)

table(sig_anno$promoter)
table(all_anno$promoter)

table(sig_anno$shelfshore)
table(all_anno$shelfshore)

table(sig_anno$cpgisland)
table(all_anno$cpgisland)

table(sig_anno$DNAse_site)
table(all_anno$DNAse_site)

#Functional Gene Enrichment Analysis - need to remove the significant CpGs from the total CpGs to get final numbers 
dat <- data.frame(
  "promoter" = c(9,115223),
  "notpromoter" = c(6, 263630),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("Promoter", "Not promoter")
fisher.test(dat)


dat2 <- data.frame(
  "island" = c(7,122343),
  "not island" = c(8, 256510),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat2) <- c("island", "Not island")
fisher.test(dat2)


dat4 <- data.frame(
  "shelf" = c(7,126070),
  "not" = c(8, 252783),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat4) <- c("Shelf", "Not shelf")
fisher.test(dat4)


dat5 <- data.frame(
  "DNAse" = c(13,242503),
  "not DNAse" = c(2, 136350),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat5) <- c("DNAse", "Not DNAse")
fisher.test(dat5)

dat6 <- data.frame(
  "enhancer" = c(1,6462),
  "notenhancer" = c(14, 372391),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat6) <- c("Enhancer", "Not enhancer")
fisher.test(dat6)

#EQTM enrichment analysis 
eqtm=read.csv("/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/eqtm.csv", header=T)
sep=separate(eqtm, col=GeneSymbol, into=c('gene1', 'gene2', 'gene3', 'gene4', 'gene5', 'gene6','gene7', 'gene8', 'gene9', 'gene10', 'gene11', 'gene12', 'gene13', 'gene14'), sep='\\|')
gene1=sep[,c(1,2)]
gene2=sep[,c(1,3)]
gene3=sep[,c(1,4)]
gene4=sep[,c(1,5)]
gene5=sep[,c(1,6)]
gene6=sep[,c(1,7)]
gene7=sep[,c(1,8)]
gene8=sep[,c(1,9)]
gene9=sep[,c(1,10)]
gene10=sep[,c(1,11)]
gene11=sep[,c(1,12)]
gene12=sep[,c(1,13)]
gene13=sep[,c(1,14)]
gene14=sep[,c(1,15)]


names(gene1)[names(gene1) == "gene1"] <- "gene"
names(gene2)[names(gene2) == "gene2"] <- "gene"
names(gene3)[names(gene3) == "gene3"] <- "gene"
names(gene4)[names(gene4) == "gene4"] <- "gene"
names(gene5)[names(gene5) == "gene5"] <- "gene"
names(gene6)[names(gene6) == "gene6"] <- "gene"
names(gene7)[names(gene7) == "gene7"] <- "gene"
names(gene8)[names(gene8) == "gene8"] <- "gene"
names(gene9)[names(gene9) == "gene9"] <- "gene"
names(gene10)[names(gene10) == "gene10"] <- "gene"
names(gene11)[names(gene11) == "gene11"] <- "gene"
names(gene12)[names(gene12) == "gene12"] <- "gene"
names(gene13)[names(gene13) == "gene13"] <- "gene"
names(gene14)[names(gene14) == "gene14"] <- "gene"

allgenes=rbind(gene1,gene2,gene3,gene4,gene5,gene6,gene7,gene8,gene9,gene10,gene11,gene12,gene13,gene14)
allgenes2=allgenes[!is.na(allgenes$gene),] 
mergedsig=merge(allgenes2,aSES450k_sig,all.y=T)

noeqtm_sig=mergedsig[is.na(mergedsig$gene),] 
noeqtm_sig2=noeqtm_sig[!duplicated(noeqtm_sig$cpg), ] #6 biomarker CpGs have no transcript pairs 

eqtm_sig=mergedsig[!is.na(mergedsig$gene),] 
eqtm_sig2=eqtm_sig[!duplicated(eqtm_sig$cpg), ] #9 biomarker CpGs have transcript pairs 

mergedall=merge(allgenes2,all_anno, all.y=T) 
noeqtm_all=mergedall[is.na(mergedall$gene),] #352295 non-biomarker CpGs have no transcript pairs (excluding significant CpGs)
eqtm_all=mergedall[!is.na(mergedall$gene),] 
eqtm_all2=eqtm_all[!duplicated(eqtm_all$cpg),] #26558 non-biomarker CpGs have transcript pairs (excluding significant CpGs)

dat7 <- data.frame(
  "eqtm" = c(9,26558),
  "noteqtm" = c(6, 352295),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat7) <- c("eQTM", "not eQTM")
fisher.test(dat7)


#GO and KEGG ENRICHMENT Analysis Childhood SES Biomarker 450K using GOmeth method 
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
sigcpgs=cSES450k_sig[,1]
allcpgs <- all_anno[,1]

gst <- gometh(sig.cpg = sigcpgs, all.cpg = allcpgs, 
              collection = "GO", array.type="450K", anno = ann)

kegg <- gometh(sig.cpg = sigcpgs, all.cpg = allcpgs, collection = "KEGG", array.type="450K", prior.prob=TRUE, anno=ann)

#GO and KEGG ENRICHMENT Analysis Adult SES Biomarker 450K using ClusterProfiler method with eQTM data
#Deduplicating all datasets
background_final=mergedall[!duplicated(mergedall$gene), ]
mergedsig_final=mergedsig[!duplicated(mergedsig$gene), ]
mergedsig_final=na.omit(mergedsig_final)
universe=background_final$gene
sig_id=mergedsig_final$gene

eg = bitr(universe, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
egsig = bitr(sig_id, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

eg2=eg$ENTREZID
egsig2=egsig$ENTREZID

ego <- enrichGO(gene          = sig_id,
                universe      = universe,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = .01,
                qvalueCutoff = .01,
                keyType = "SYMBOL")
summary_sig=data.frame(ego)
sum_GO=as.data.frame(s_ego)
dotplot(s_ego, font.size = 8)


kk <- enrichKEGG(gene         = egsig2,
                 universe =eg2,
                 organism     = 'hsa',
                 keyType = "kegg",
                 pvalueCutoff = 0.01,
                 qvalueCutoff=0.01)
summary_keg=data.frame(kk)
