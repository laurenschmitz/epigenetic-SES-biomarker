library(lme4)
library(lmerTest)
library(data.table)




#input bash arguments
args = commandArgs(trailingOnly=TRUE)
start = as.numeric(args[1])
filenum = start/100000 + 1
ses = args[2]
strata = args[3]
model = args[4]



#import methylation data
header = read.csv("/home/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer04.13.2021/allbeta.csv", header=F, nrows=1)
meth = read.csv("/home/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer04.13.2021/allbeta.csv", header=F, nrows=100000, skip=(start+1))

#transpose methylation data
methy = data.frame(id=unlist(header)[-1], t(meth[,-c(1)]), row.names = NULL)
colnames(methy)[-1] = as.vector(meth$V1)
rm(meth)

#remove sites with more than 200 missing values
methy2 = methy[,colSums(is.na(methy))<201]
rm(methy)


n = dim(methy2)[2]-1
cpg = colnames(methy2)[-1]



#import phenotype data created in HRS phenos.sas code and factorize appropriate variables
ptype = read.csv('/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/data/HRS_methyset_phenos.csv', header=T)

ptype$race=factor(ptype$race)
ptype$smoke=factor(ptype$smoke)
ptype$sex = ifelse(ptype$RAGENDER==2, 'F', ifelse(ptype$RAGENDER==1, 'M', NA))
ptype$SMOKE2=ifelse(ptype$smoke==2, 1, ifelse(ptype$smoke %in% c(0,1), 0, NA))
ptype$row = as.factor(substr(ptype$BCD_Well, 1, 1))
ptype$col = as.factor(substr(ptype$BCD_Well, 2, 3))



#PCs
PC.comb = read.csv('/home/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Genotype_Jan2015_Processing/PCA/Phase1234/Top_PCs.csv')[,1:5]
colnames(PC.comb) = c('local_id', 'PC1', 'PC2', 'PC3', 'PC4')

PC.EA = read.table('/home/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Genotype_Jan2015_Processing/Analysis/EA/PCA/EA_PCA_top_10.eigenvec', header=F)[,2:6]
colnames(PC.EA) = c('local_id', 'PC1', 'PC2', 'PC3', 'PC4')

PC.AA = read.table('/home/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Genotype_Jan2015_Processing/Analysis/AA/PCA/AA_PCA_top_10.eigenvec', header=F)[,2:6]
colnames(PC.AA) = c('local_id', 'PC1', 'PC2', 'PC3', 'PC4')

PC.HA = read.table('/home/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Genotype_Jan2015_Processing/Analysis/Hispanic/PCA/Hispanic_PCA_top_10.eigenvec', header=F)[,2:6]
colnames(PC.HA) = c('local_id', 'PC1', 'PC2', 'PC3', 'PC4')





#combine phenotypes with appropriate PCs and statification
if (strata=='ALL'){
  data = merge(ptype, PC.comb, by='local_id')

} else if (strata %in% c('M', 'F')){
  data = subset(merge(ptype, PC.comb, by='local_id'), sex==strata)  

} else if (strata=='1'){
  data = subset(merge(ptype, PC.EA, by='local_id'), ETHNICITY==strata)

} else if (strata=='2'){
  data = subset(merge(ptype, PC.AA, by='local_id'), ETHNICITY==strata)
 
} else if (strata=='3'){
  data = subset(merge(ptype, PC.HA, by='local_id'), ETHNICITY==strata)

} else if (strata=='M1'){
  data = subset(merge(ptype, PC.EA, by='local_id'), sex=='M' & ETHNICITY=='1')

} else if (strata=='M2'){
  data = subset(merge(ptype, PC.AA, by='local_id'), sex=='M' & ETHNICITY=='2')
  
} else if (strata=='M3'){
  data = subset(merge(ptype, PC.HA, by='local_id'), sex=='M' & ETHNICITY=='3')
  
} else if (strata=='F1'){
  data = subset(merge(ptype, PC.EA, by='local_id'), sex=='F' & ETHNICITY=='1')
  
} else if (strata=='F2'){
  data = subset(merge(ptype, PC.AA, by='local_id'), sex=='F' & ETHNICITY=='2')
  
} else if (strata=='F3'){
  data = subset(merge(ptype, PC.HA, by='local_id'), sex=='F' & ETHNICITY=='3')
  
}






#combine phenotypes with methylation
hrs = merge(data, methy2, by.x='sample', by.y='id')







#############EWAS MODELS##########################
snpmod = function(out) {
  if (ses %in% c('factor_ana_1216', 'F1_PC2_1216', 'socenv', 'occ_prest', 'lnhhincx', 'lnwealth2x')){
    if (model=='1' & strata=='ALL'){
      temp = na.omit(hrs[,c(out, ses, 'PAGE', 'sex', 'race', 'NK', 'B', 'CD4', 'CD8', 'MO', 'RAEDYRS', 'row', 'col', 'Sample_Plate', 'PC1', 'PC2', 'PC3', 'PC4')])
      eq = as.formula(parse(text=paste0(out, ' ~ ', ses, ' + PAGE + sex + race + NK + B + CD4 + CD8 + MO + RAEDYRS + PC1 + PC2 + PC3 + PC4 + (1|row) + (1|col) + (1|Sample_Plate)'))[[1]])
    } 
    else if (model=='2' & strata=='ALL'){
      temp = na.omit(hrs[,c(out, ses, 'PAGE', 'sex', 'race', 'NK', 'B', 'CD4', 'CD8', 'MO', 'RAEDYRS', 'smoke', 'row', 'col', 'Sample_Plate', 'PC1', 'PC2', 'PC3', 'PC4')])
      eq = as.formula(parse(text=paste0(out, ' ~ ', ses, ' + PAGE + sex + race + NK + B + CD4 + CD8 + MO + RAEDYRS + smoke + PC1 + PC2 + PC3 + PC4 + (1|row) + (1|col) + (1|Sample_Plate)'))[[1]])
    }
    else if (model=='1' & strata %in% c('1', '2', '3')){
      temp = na.omit(hrs[,c(out, ses, 'AGE', 'SEX', 'NK', 'B', 'CD4', 'CD8', 'MO', 'EDYEARS', 'row', 'col', 'Sample_Plate', 'PC1', 'PC2', 'PC3', 'PC4')])
      eq = as.formula(parse(text=paste0(out, ' ~ ', ses, ' + AGE + SEX + NK + B + CD4 + CD8 + MO + EDYEARS + PC1 + PC2 + PC3 + PC4 + (1|row) + (1|col) + (1|Sample_Plate)'))[[1]])
    }
    else if (model=='2' & strata %in% c('1', '2', '3')){
      temp = na.omit(hrs[,c(out, ses, 'AGE', 'SEX', 'NK', 'B', 'CD4', 'CD8', 'MO', 'EDYEARS', 'SMOKE', 'row', 'col', 'Sample_Plate', 'PC1', 'PC2', 'PC3', 'PC4')])
      eq = as.formula(parse(text=paste0(out, ' ~ ', ses, ' + AGE + SEX + NK + B + CD4 + CD8 + MO + EDYEARS + SMOKE + PC1 + PC2 + PC3 + PC4 + (1|row) + (1|col) + (1|Sample_Plate)'))[[1]])
    }
    else if (model=='1' & strata %in% c('M', 'F')){
      temp = na.omit(hrs[,c(out, ses, 'AGE', 'ETHNICITY', 'NK', 'B', 'CD4', 'CD8', 'MO', 'EDYEARS', 'row', 'col', 'Sample_Plate', 'PC1', 'PC2', 'PC3', 'PC4')])
      eq = as.formula(parse(text=paste0(out, ' ~ ', ses, ' + AGE + ETHNICITY + NK + B + CD4 + CD8 + MO + EDYEARS + PC1 + PC2 + PC3 + PC4 + (1|row) + (1|col) + (1|Sample_Plate)'))[[1]])
    }
    else if (model=='2' & strata %in% c('M', 'F')){
      temp = na.omit(hrs[,c(out, ses, 'AGE', 'ETHNICITY', 'NK', 'B', 'CD4', 'CD8', 'MO', 'EDYEARS', 'SMOKE', 'row', 'col', 'Sample_Plate', 'PC1', 'PC2', 'PC3', 'PC4')])
      eq = as.formula(parse(text=paste0(out, ' ~ ', ses, ' + AGE + ETHNICITY + NK + B + CD4 + CD8 + MO + EDYEARS + SMOKE + PC1 + PC2 + PC3 + PC4 + (1|row) + (1|col) + (1|Sample_Plate)'))[[1]])
    }
    else if (model=='2' & strata %in% c('M1', 'F1', 'M2', 'F2', 'M3', 'F3')){
      temp = na.omit(hrs[,c(out, ses, 'AGE', 'NK', 'B', 'CD4', 'CD8', 'MO', 'EDYEARS', 'SMOKE', 'row', 'col', 'Sample_Plate', 'PC1', 'PC2', 'PC3', 'PC4')])
      eq = as.formula(parse(text=paste0(out, ' ~ ', ses, ' + AGE + NK + B + CD4 + CD8 + MO + EDYEARS + SMOKE + PC1 + PC2 + PC3 + PC4 + (1|row) + (1|col) + (1|Sample_Plate)'))[[1]])
    }
  }
  else {
    if (model=='1' & strata=='ALL'){
      temp = na.omit(hrs[,c(out, ses, 'AGE', 'SEX', 'ETHNICITY', 'NK', 'B', 'CD4', 'CD8', 'MO', 'row', 'col', 'Sample_Plate', 'PC1', 'PC2', 'PC3', 'PC4')])
      eq = as.formula(parse(text=paste0(out, ' ~ ', ses, ' + AGE + SEX + ETHNICITY + NK + B + CD4 + CD8 + MO + PC1 + PC2 + PC3 + PC4 + (1|row) + (1|col) + (1|Sample_Plate)'))[[1]])
    } 
    else if (model=='2' & strata=='ALL'){
      temp = na.omit(hrs[,c(out, ses, 'PAGE', 'sex', 'race', 'NK', 'B', 'CD4', 'CD8', 'MO', 'smoke', 'row', 'col', 'Sample_Plate', 'PC1', 'PC2', 'PC3', 'PC4')])
      eq = as.formula(parse(text=paste0(out, ' ~ ', ses, ' + PAGE + sex + race + NK + B + CD4 + CD8 + MO + smoke + PC1 + PC2 + PC3 + PC4 + (1|row) + (1|col) + (1|Sample_Plate)'))[[1]])
    }
    else if (model=='1' & strata %in% c('1', '2', '3')){
      temp = na.omit(hrs[,c(out, ses, 'AGE', 'SEX', 'NK', 'B', 'CD4', 'CD8', 'MO', 'row', 'col', 'Sample_Plate', 'PC1', 'PC2', 'PC3', 'PC4')])
      eq = as.formula(parse(text=paste0(out, ' ~ ', ses, ' + AGE + SEX + NK + B + CD4 + CD8 + MO + PC1 + PC2 + PC3 + PC4 + (1|row) + (1|col) + (1|Sample_Plate)'))[[1]])
    }
    else if (model=='2' & strata %in% c('1', '2', '3')){
      temp = na.omit(hrs[,c(out, ses, 'AGE', 'SEX', 'NK', 'B', 'CD4', 'CD8', 'MO', 'SMOKE', 'row', 'col', 'Sample_Plate', 'PC1', 'PC2', 'PC3', 'PC4')])
      eq = as.formula(parse(text=paste0(out, ' ~ ', ses, ' + AGE + SEX + NK + B + CD4 + CD8 + MO + SMOKE + PC1 + PC2 + PC3 + PC4 + (1|row) + (1|col) + (1|Sample_Plate)'))[[1]])
    }
    else if (model=='1' & strata %in% c('M', 'F')){
      temp = na.omit(hrs[,c(out, ses, 'AGE', 'ETHNICITY', 'NK', 'B', 'CD4', 'CD8', 'MO', 'row', 'col', 'Sample_Plate', 'PC1', 'PC2', 'PC3', 'PC4')])
      eq = as.formula(parse(text=paste0(out, ' ~ ', ses, ' + AGE + ETHNICITY + NK + B + CD4 + CD8 + MO + PC1 + PC2 + PC3 + PC4 + (1|row) + (1|col) + (1|Sample_Plate)'))[[1]])
    }
    else if (model=='2' & strata %in% c('M', 'F')){
      temp = na.omit(hrs[,c(out, ses, 'AGE', 'ETHNICITY', 'NK', 'B', 'CD4', 'CD8', 'MO', 'SMOKE', 'row', 'col', 'Sample_Plate', 'PC1', 'PC2', 'PC3', 'PC4')])
      eq = as.formula(parse(text=paste0(out, ' ~ ', ses, ' + AGE + ETHNICITY + NK + B + CD4 + CD8 + MO + SMOKE + PC1 + PC2 + PC3 + PC4 + (1|row) + (1|col) + (1|Sample_Plate)'))[[1]])
    }
    else if (model=='2' & strata %in% c('M1', 'F1', 'M2', 'F2', 'M3', 'F3')){
      temp = na.omit(hrs[,c(out, ses, 'AGE', 'NK', 'B', 'CD4', 'CD8', 'MO', 'SMOKE', 'row', 'col', 'Sample_Plate', 'PC1', 'PC2', 'PC3', 'PC4')])
      eq = as.formula(parse(text=paste0(out, ' ~ ', ses, ' + AGE + NK + B + CD4 + CD8 + MO + SMOKE + PC1 + PC2 + PC3 + PC4 + (1|row) + (1|col) + (1|Sample_Plate)'))[[1]])
    }
  }
 
  lmer_result = summary(lmer(eq, data=temp))
  
  results=rep(NA,7)
  results[1] = out
  results[2] = dim(temp)[1]
  results[3:7] = lmer_result$coefficients[ses, c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")]

  return (results) 
}





path = '/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/results/HRS Enet/EWAS/'



model1 = matrix(NA, n, 7)
for (i in 1:n){
  out = cpg[i]
  model1[i,] = snpmod(out)
  print(i)
}

colnames(model1) = c('CPG', 'N', 'BETA', 'SE', 'DF', 't', 'P')




#import annotation file and combine with model results
anno = read.csv('/home/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/QC_processing_Smith/infinium-methylationepic-v-1-0-b5-manifest-file NOEXTRAS.csv', header=T)[,c(1,12,13)]
colnames(anno) = c('CPG', 'CHR', 'POS')

final = merge(anno, model1, by='CPG')


strata2 = ifelse(strata=='1', 'EA', ifelse(strata=='2', 'AA', ifelse(strata=='3', 'HA', strata)))


output = paste0(path, 'HRS_', ses, '_', strata2, '_Model', model, '_', filenum, '.csv')
write.csv(final, output , row.names=F, quote=F)







