library(lme4)
library(lmerTest)
library(data.table)
library(zoo)




#input bash arguments
args = commandArgs(trailingOnly=TRUE)
start = as.numeric(args[1])
filenum = start/100000 + 1



#import methylation data
header = read.csv("/home/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer04.13.2021/allbeta.csv", header=F, nrows=1)
meth = read.csv("/home/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer04.13.2021/allbeta.csv", header=F, nrows=100000, skip=(start+1))

#transpose methylation data
methy = data.frame(sample=unlist(header)[-1], t(meth[,-c(1)]), row.names = NULL)
colnames(methy)[-1] = as.vector(meth$V1)
rm(meth)

#remove sites with more than 200 missing values
methy2 = methy[,colSums(is.na(methy))<201]
rm(methy)


#impute missing methylation values to the mean for each site
methyimp = na.aggregate(as.matrix(methy2[-1]))
methy3 = data.frame(methy2[1], methyimp)
rm(methy2)
rm(methyimp)


n = dim(methy3)[2]-1
cpg = colnames(methy3)[-1]




#import phenotype data created in HRS phenos.sas code and factorize appropriate variables
phenos = read.csv('/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/data/HRS_methyset_phenos.csv', header=T)
phenos$smoke = as.factor(phenos$smoke)
phenos$peduc = as.factor(phenos$peduc)
phenos$drinkscat = as.factor(phenos$drinkscat)
phenos$RAGENDER = as.factor(phenos$RAGENDER)
phenos$race = as.factor(phenos$race)
phenos$Slide = as.factor(phenos$Slide)
phenos$row = as.factor(substr(phenos$BCD_Well, 1, 1))
phenos$col = as.factor(substr(phenos$BCD_Well, 2, 3))



#PCs
PC.comb = read.csv('/home/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Genotype_Jan2015_Processing/PCA/Phase1234/Top_PCs.csv')[,1:5]
colnames(PC.comb) = c('local_id', 'PC1', 'PC2', 'PC3', 'PC4')





#combine phenotypes, Pcs, and methylation
merge = merge(merge(phenos, PC.comb, by='local_id', all.x=T), methy3, by='sample')
#mergex = subset(merge, !is.na(merge$race))


#######USE ONE OF THE FOLLOWING CHUNKS OF CODE CORRESPONDING TO MODEL YOU WANT###############

#mergex1 = data.frame(mergex[1])

#for (i in 1:n){
#  M = cpg[i]
#  temp=mergex[,c(M,'race', 'PAGE', 'RAGENDER', 'NK', 'B', 'CD4', 'CD8', 'MO', 'PC1', 'PC2', 'PC3', 'PC4', 'row', 'col', 'Sample_Plate')]
#  MODEL = as.formula(parse(text=paste(M,' ~ race + PAGE + RAGENDER + NK + B + CD4 + CD8 + MO + PC1 + PC2 + PC3 + PC4 + (1|row) + (1|col) + (1|Sample_Plate)', sep=''))[[1]])
#  mergex1[i+1] = residuals(lmer(MODEL, data=temp, na.action=na.exclude))
#  colnames(mergex1)[i+1] = M
#  print(i)
#}

#fileout = paste0('/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/data/HRS_methyation_resids_', filenum, '.csv')
#write.csv(mergex1, fileout, row.names=F, quote=F)







#mergex = subset(merge, !is.na(merge$race) & !is.na(merge$RAEDYRS) & !is.na(merge$smoke))

#mergex2 = data.frame(mergex[1])

#for (i in 1:n){
#  M = cpg[i]
#  temp=mergex[,c(M,'race', 'PAGE', 'RAGENDER', 'RAEDYRS', 'smoke', 'NK', 'B', 'CD4', 'CD8', 'MO', 'PC1', 'PC2', 'PC3', 'PC4', 'row', 'col', 'Sample_Plate')]
#  MODEL = as.formula(parse(text=paste(M,' ~ race + PAGE + RAGENDER + RAEDYRS + smoke + NK + B + CD4 + CD8 + MO + PC1 + PC2 + PC3 + PC4 + (1|row) + (1|col) + (1|Sample_Plate)', sep=''))[[1]])
#  mergex2[i+1] = residuals(lmer(MODEL, data=temp, na.action=na.exclude))
#  colnames(mergex2)[i+1] = M
#  print(i)
#}

#fileout = paste0('/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/data/HRS_methyation_resids2_', filenum, '.csv')
#write.csv(mergex2, fileout, row.names=F, quote=F)








#mergex = subset(merge, !is.na(merge$race) & !is.na(merge$smoke))

#mergex3 = data.frame(mergex[1])

#for (i in 1:n){
#  M = cpg[i]
#  temp=mergex[,c(M,'race', 'PAGE', 'RAGENDER', 'smoke', 'NK', 'B', 'CD4', 'CD8', 'MO', 'PC1', 'PC2', 'PC3', 'PC4', 'row', 'col', 'Sample_Plate')]
#  MODEL = as.formula(parse(text=paste(M,' ~ race + PAGE + RAGENDER + smoke + NK + B + CD4 + CD8 + MO + PC1 + PC2 + PC3 + PC4 + (1|row) + (1|col) + (1|Sample_Plate)', sep=''))[[1]])
#  mergex3[i+1] = residuals(lmer(MODEL, data=temp, na.action=na.exclude))
#  colnames(mergex3)[i+1] = M
#  print(i)
#}

#fileout = paste0('/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/data/HRS_methyation_resids3_', filenum, '.csv')
#write.csv(mergex3, fileout, row.names=F, quote=F)







#mergex = subset(merge, !is.na(merge$race) & !is.na(merge$smoke))

#mergex4 = data.frame(mergex[1])

#for (i in 1:n){
#  M = cpg[i]
#  temp=mergex[,c(M,'race', 'PAGE', 'RAGENDER', 'smoke', 'NK', 'B', 'CD4', 'CD8', 'MO', 'row', 'col', 'Sample_Plate')]
#  MODEL = as.formula(parse(text=paste(M,' ~ race + PAGE + RAGENDER + smoke + NK + B + CD4 + CD8 + MO + (1|row) + (1|col) + (1|Sample_Plate)', sep=''))[[1]])
#  mergex4[i+1] = residuals(lmer(MODEL, data=temp, na.action=na.exclude))
#  colnames(mergex4)[i+1] = M
#  print(i)
#}

#fileout = paste0('/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/data/HRS_methyation_resids4_', filenum, '.csv')
#write.csv(mergex4, fileout, row.names=F, quote=F)






mergex = subset(merge, !is.na(merge$race) & !is.na(merge$smoke))

mergex3 = data.frame(mergex[1])

for (i in 1:n){
  M = cpg[i]
  temp=mergex[,c(M,'race', 'PAGE', 'RAGENDER', 'smoke', 'NK', 'B', 'CD4', 'CD8', 'MO', 'PC1', 'PC2', 'PC3', 'PC4', 'row', 'col', 'Sample_Plate')]
  MODEL = as.formula(parse(text=paste(M,' ~ race + PAGE + RAGENDER + smoke + NK + B + CD4 + CD8 + MO + PC1 + PC2 + PC3 + PC4 + (1|row) + (1|col) + (1|Sample_Plate)', sep=''))[[1]])
  mergex3[i+1] = residuals(lmer(MODEL, data=temp, na.action=na.exclude)) + mean(temp[,M])
  colnames(mergex3)[i+1] = M
  print(i)
}

fileout = paste0('/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/data/HRS_methyation_resids5_', filenum, '.csv')
write.csv(mergex3, fileout, row.names=F, quote=F)


