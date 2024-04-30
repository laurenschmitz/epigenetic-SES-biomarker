
library(glmnet)
library(lme4)
library(lmerTest)
library(zoo)
library(pROC)
library(data.table)



#set SES outcome variable
out='SESindex'


#read in phenotype data and convert factors
phenos = read.csv('/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/data/HRS_methyset_phenos.csv', header=T)
phenos$smoke = as.factor(phenos$smoke)
phenos$peduc = as.factor(phenos$peduc)
phenos$drinkscat = as.factor(phenos$drinkscat)
phenos$RAGENDER = as.factor(phenos$RAGENDER)
phenos$race = as.factor(phenos$race)
phenos$Slide = as.factor(phenos$Slide)
phenos$SMOKE2=ifelse(phenos$smoke==2, 1, ifelse(phenos$smoke %in% c(0,1), 0, NA))


#####ALTERNATIVELY DONE IN RESIDUAL PROGRAM#########
#header = read.csv("/home/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer04.13.2021/allbeta.csv", header=F, nrows=1)
#meth = read.csv("/home/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer04.13.2021/allbeta.csv", header=F, skip=1)

#methy = data.frame(sample=unlist(header)[-1], t(meth[,-c(1)]), row.names = NULL)
#colnames(methy)[-1] = as.vector(meth$V1)
#rm(meth)



#png(filename = '/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/results/HRS Enet/histogram_cg17486446.png' , width = 1500, height = 1500,bg = "white",pointsize=24)
#hist(methy$cg17486446)
#dev.off()



#methy2 = methy[,colSums(is.na(methy))<201]
#rm(methy)



#methyimp = na.aggregate(as.matrix(methy2[-1]))
#methy3 = data.frame(methy2[1], methyimp)
#rm(methy2)



#methy3 = read.csv('/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/data/HRS_methyation_resids4_1.csv', header=T)

#for (i in 2:9){
#  file = paste0('/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/data/HRS_methyation_resids4_', i, '.csv')
#  resid = read.csv(file, header=T)[-1]
#  methy3 = cbind(methy3, resid)
#} 

#save(methy3, file='/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/data/HRS_methyation_resids4.RData')
##############################


#load the correct version of residuals
load('/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/data/HRS_methyation_resids3.RData')




#read in EWAS results and filter methylation data to significant CpGs
ewas = read.csv(paste0('/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/results/HRS Enet/EWAS/HRS_', out, '_ALL_Model2.csv'), header=T)
ewasx = subset(ewas, P<0.05)
ewasx$CPG = as.character(ewasx$CPG)

#remove cross-reactive probes
xrprobes = read.csv('/home/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/QC_processing_Smith/EPICcrossreactiveprobelist.csv', header=T)[,c(1,7)]
colnames(xrprobes) = c('CPG', 'XRprobe')
ewasx = subset(ewasx, !(CPG %in% xrprobes$CPG))

methy3x = data.frame(methy3[1], methy3[,ewasx$CPG])



#merge filtered methyaltion data with phenotypes
merge = merge(methy3x, phenos, by='sample')
mergex = subset(merge, !is.na(merge[out]))




cpg = names(methy3x)[-1]
n=dim(methy3x)[2]-1

#rm(methy3)




#####ENET analysis##########  

#be sure to set the appropriate family for the given outcome
set.seed(765432)
CV = cv.glmnet(data.matrix(mergex[,2:(n+1)]),mergex[,out], nfolds=5,alpha=0.5, family="gaussian", keep=T)
coefs = as.data.frame(as.matrix(coef(CV, s=CV$lambda.min)))

coefsx = subset(coefs, s1!=0)



write.csv(coefsx, paste0('/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/results/HRS Enet/HRS_Enet_XValid_', out, '_betas_lambdamin_final.csv'), quote=F)


pred = predict(CV, data.matrix(mergex[,2:(n+1)]), s=CV$lambda.min)



#coefs450 = coefs
#coefs450$s1 = ifelse(row.names(coefs450) %in% c('cg21327712','cg11047325','cg23248055','cg23491443','cg03362418','cg03077610','cg10404730','cg24445316'), 0, coefs450$s1)
#coefs450x = subset(coefs450, s1!=0)
#pred450 = as.vector(as.matrix(mergex[,2:(n+1)])%*%as.matrix(coefs450)[-1] + coefs450$s1[1])
#clockname450=paste0(out, '_mclock450')
#mergex[, clockname450] = pred450
#write.csv(mergex[,c('hhid', 'PN', clockname450)], paste0('/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/results/HRS Enet/HRS_', clockname450, '.csv'), quote=F, row.names=F)




rr = cor(pred, mergex[,out])

png(filename = paste0('/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/results/HRS Enet/HRS ',out, ' XValid Elastic Net predicted_v_actual Plot Final.png') , width = 1800, height = 1500,bg = "white",pointsize=24)
par(mar=c(5,6,4,2)+.1)
plot(pred, mergex[,out], main=paste0("XValidated Elastic Net Predicted ", out, " vs. Actual \n r=", round(rr,3)), cex.main=2, cex.lab=2, cex.axis=2, xlab='Predicted SES Index', ylab='Actual SES Index')
abline(lm(mergex[,out] ~ pred), col="red")
dev.off()


#roc = roc(response=mergex[,out], predictor=pred, ci=T)

#png(filename = paste0('/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/results/HRS Enet/HRS ',out, ' XValid Elastic Net predicted ROC.png') , width = 1500, height = 1500,bg = "white",pointsize=24)
#plot.roc(roc, main=paste0("XValidated Elastic Net Predicted ", out, " ROC"), print.auc=T, print.thres=T)
#dev.off()


clockname=paste0(out, '_mclock')
mergex[, clockname] = pred
write.csv(mergex[,c('hhid', 'PN', clockname)], paste0('/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/results/HRS Enet/HRS_', clockname, '.csv'), quote=F, row.names=F)






cvm = data.frame(lambda=CV$lambda, cvm=CV$cvm)
write.csv(cvm, paste0('/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/results/HRS Enet/HRS_Enet_XValid_', out, '_lambdas_cvm.csv'), quote=F, row.names=F)









#create 450K subset from valid MESA probes
mesacpgs = read.csv('/home/skardia_lab/clubhouse/research/projects/MESA/Methylation_081413/ValidProbe_DMRcate_MESA.csv', header=T, stringsAsFactors=F)

ewasx450 = subset(ewasx, CPG %in% mesacpgs$x)
methy3x450 = data.frame(methy3[1], methy3[,ewasx450$CPG])


#merge 450k filtered methylation data with phenotypes
merge = merge(methy3x450, phenos, by='sample')
mergex450 = subset(merge, !is.na(merge[out]))



cpg = names(methy3x450)[-1]
n=dim(methy3x450)[2]-1



######450K ENET ANALYSIS##############

set.seed(765432)   #SESindex
set.seed(9999999)  #cSES
CV = cv.glmnet(data.matrix(mergex450[,2:(n+1)]),mergex450[,out], nfolds=5,alpha=0.5, family="gaussian", keep=T)
coefs = as.data.frame(as.matrix(coef(CV, s=CV$lambda.min)))

coefsx = subset(coefs, s1!=0)


write.csv(coefsx, paste0('/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/results/HRS Enet/HRS_Enet_XValid_', out, '_betas_lambdamin_450K.csv'), quote=F)


pred = predict(CV, data.matrix(mergex450[,2:(n+1)]), s=CV$lambda.min)
rr = cor(pred, mergex[,out])


png(filename = paste0('/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/results/HRS Enet/HRS ',out, ' XValid Elastic Net predicted_v_actual Plot Final 450K.png') , width = 1800, height = 1500,bg = "white",pointsize=24)
par(mar=c(5,6,4,2)+.1)
plot(pred, mergex450[,out], main=paste0("XValidated Elastic Net Predicted ", out, " vs. Actual \n r=", round(rr,3)), cex.main=2, cex.lab=2, cex.axis=2, xlab='Predicted SES Index', ylab='Actual SES Index')
abline(lm(mergex[,out] ~ pred), col="red")
dev.off()



clockname=paste0(out, '_mclock')
mergex450[, clockname] = pred
write.csv(mergex450[,c('hhid', 'PN', clockname)], paste0('/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/results/HRS Enet/HRS_', clockname, '_450K.csv'), quote=F, row.names=F)



cvm = data.frame(lambda=CV$lambda, cvm=CV$cvm)
write.csv(cvm, paste0('/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/results/HRS Enet/HRS_Enet_XValid_', out, '_lambdas_cvm_450K.csv'), quote=F, row.names=F)







########Does biomarker predict in MESA?################

mesameth = as.data.frame(fread("gunzip -c /net/orion/skardia_lab/clubhouse/research/projects/MESA/Methylation_081413/MESA_Epi_METH_data.txt.gz", header=T))
weights = read.csv(paste0('/net/orion/skardia_lab/treehouse/science/projects/MESA/Methylation/Lauren_DNAmAge/results/HRS_Enet/HRS_Enet_XValid_', out, '_betas_lambdamin_450K.csv'), header=T)
weights = weights[order(weights$X),]

mesamethx = subset(mesameth, TargetID %in% weights$X)
#remove(mesameth)
mesamethx = mesamethx[order(mesamethx$TargetID),]
mesamethx2 = mesamethx[,-c(1)]
mesamethy = data.frame(idno=colnames(mesamethx2), t(mesamethx2), row.names=NULL)
colnames(mesamethy)[-1] = mesamethx$TargetID
#methy$idno = as.numeric(rownames(methy))
#remove(mesamethx)
#remove(mesamethx2)



mesamethy = cbind(mesamethy[1], 2**(mesamethy[-1])/(1+2**mesamethy[-1]))


cpg = names(mesamethy)[-1]
n=dim(mesamethy)[2]-1


#import MESA phenotypes
ptype = read.csv('/net/orion/skardia_lab/treehouse/science/projects/MESA/Methylation/Lauren_DNAmAge/data/MESA_clock_phenos_v3.csv', header=T)
ptype$chip_meth = as.factor(ptype$chip_meth)

#PCs
PC.comb = read.table('/net/orion/skardia_lab/clubhouse/research/projects/MESA/MESA_EPI_PCA/Combined.txt')[,1:5]
colnames(PC.comb) = c('idno', 'PC1', 'PC2', 'PC3', 'PC4')

data = merge(ptype, PC.comb, by='idno')

mesa = merge(mesamethy, data, by='idno')
mesax = subset(mesa, !is.na(mesa$smoke))



#Residualize MESA methylation beta values
mesares = data.frame(mesax[1])

for (i in 1:n){
  M = cpg[i]
  temp=mesax[,c(M,'race', 'age5c', 'gender', 'smoke', 'bcell', 'tcell', 'nkcell', 'neutro', 'PC1', 'PC2', 'PC3', 'PC4', 'chip_meth', 'pos_meth')]
  MODEL = as.formula(parse(text=paste(M,' ~ race + age5c + gender + smoke + bcell + tcell + nkcell + neutro + PC1 + PC2 + PC3 + PC4 + (1|chip_meth) + (1|pos_meth)', sep=''))[[1]])
  mesares[i+1] = residuals(lmer(MODEL, data=temp, na.action=na.exclude))
  colnames(mesares)[i+1] = M
  print(i)
}


mesa2 = merge(mesares, data, by='idno')



pred = as.vector(as.matrix(mesa2[,2:(n+1)])%*%as.matrix(weights$s1)[-1] + mean(mesa2$SESindex, na.rm=T))
rr = cor(pred, mesa2[,out], use='complete.obs')

pred2 = as.vector(as.matrix(mesax[,2:(n+1)])%*%as.matrix(weights$s1)[-1] + mean(mesax[,out], na.rm=T))
rr2 = cor(pred2, mesax[,out], use='complete.obs')


png(filename = paste0('/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/results/HRS Enet/HRS ',out, ' Elastic Net 450K predicting MESA.png') , width = 1800, height = 1500,bg = "white",pointsize=24)
par(mar=c(5,6,4,2)+.1)
plot(pred, mesa2[,out], main=paste0("XValidated Elastic Net Predicted ", out, " vs. Actual \n r=", round(rr,3)), cex.main=2, cex.lab=2, cex.axis=2, xlab='Predicted SES Index', ylab='Actual SES Index')
abline(lm(mesa2[,out] ~ pred), col="red")
dev.off()

png(filename = paste0('/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/results/HRS Enet/HRS ',out, ' Elastic Net 450K predicting MESA RAW BETAS.png') , width = 1800, height = 1500,bg = "white",pointsize=24)
par(mar=c(5,6,4,2)+.1)
plot(pred2, mesax[,out], main=paste0("XValidated Elastic Net Predicted ", out, " vs. Actual \n r=", round(rr2,3)), cex.main=2, cex.lab=2, cex.axis=2, xlab='Predicted SES Index', ylab='Actual SES Index')
abline(lm(mesax[,out] ~ pred2), col="red")
dev.off()


clockname=paste0(out, '_mclock')
mesa2[, clockname] = pred
write.csv(mesa2[,c('idno', clockname)], paste0('/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/results/HRS_Enet/MESA_', clockname, '_450K.csv'), quote=F, row.names=F)

mesax[, clockname] = pred2
write.csv(mesax[,c('idno', clockname)], paste0('/net/orion/skardia_lab/treehouse/science/projects/MESA/Methylation/Lauren_DNAmAge/results/HRS_Enet/MESA_', clockname, '_450K_raw.csv'), quote=F, row.names=F)



mesatest = lm(SESindex ~ SESindex_mclock, data=mesa2)

sink('/treehouse/skardia_lab/science/projects/MESA/Methylation/Lauren_DNAmAge/results/HRS Enet/MESA_SESindex_Biomarker_Regression.txt')
summary(mesatest)
sink()











