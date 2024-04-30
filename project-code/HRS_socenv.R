library(lme4)
library(lmerTest)
library(sas7bdat)


#2006
lb = read.csv('S:/MESA/Methylation/Lauren_DNAmAge/data/lb06.csv', header=T)


so_co = lmer(so_co ~ KAGE + GENDER + as.factor(STFIPS10) + (1|LINKCEN2010) + (1|LINKCEN2010:hhidpn), data=lb)
aes_qual = lmer(aes_qual ~ KAGE + GENDER + as.factor(STFIPS10) + (1|LINKCEN2010) + (1|LINKCEN2010:hhidpn), data=lb)
safe = lmer(safe ~ KAGE + GENDER + as.factor(STFIPS10) + (1|LINKCEN2010), data=lb)


so_co_re = ranef(so_co)[['LINKCEN2010']]
so_co_re$LINKCEN2010 = ifelse(nchar(row.names(so_co_re))==11, row.names(so_co_re), paste0('0', row.names(so_co_re)))
so_co_re = data.frame(so_co_re[2], so_co=so_co_re[,1], row.names=NULL)                           

aes_qual_re = ranef(aes_qual)[['LINKCEN2010']]
aes_qual_re$LINKCEN2010 = ifelse(nchar(row.names(aes_qual_re))==11, row.names(aes_qual_re), paste0('0', row.names(aes_qual_re)))
aes_qual_re = data.frame(aes_qual_re[2], aes_qual=aes_qual_re[,1], row.names=NULL)                           

safe_re = ranef(safe)[['LINKCEN2010']]
safe_re$LINKCEN2010 = ifelse(nchar(row.names(safe_re))==11, row.names(safe_re), paste0('0', row.names(safe_re)))
safe_re = data.frame(safe_re[2], safe=safe_re[,1], row.names=NULL)    


socenv06 = merge(merge(so_co_re, aes_qual_re, by="LINKCEN2010", all=T), safe_re, by="LINKCEN2010", all=T)
socenv06$so_co = scale(socenv06$so_co)
socenv06$aes_qual = scale(socenv06$aes_qual)
socenv06$safe = scale(socenv06$safe)
socenv06$socenv = rowSums(socenv06[2:4], na.rm=T)
socenv06$year = 2006 






#2008
lb = read.csv('S:/MESA/Methylation/Lauren_DNAmAge/data/lb08.csv', header=T)


so_co = lmer(so_co ~ LAGE + GENDER + as.factor(STFIPS10) + (1|LINKCEN2010) + (1|LINKCEN2010:hhidpn), data=lb)
aes_qual = lmer(aes_qual ~ LAGE + GENDER + as.factor(STFIPS10) + (1|LINKCEN2010) + (1|LINKCEN2010:hhidpn), data=lb)
safe = lmer(safe ~ LAGE + GENDER + as.factor(STFIPS10) + (1|LINKCEN2010), data=lb)


so_co_re = ranef(so_co)[['LINKCEN2010']]
so_co_re$LINKCEN2010 = ifelse(nchar(row.names(so_co_re))==11, row.names(so_co_re), paste0('0', row.names(so_co_re)))
so_co_re = data.frame(so_co_re[2], so_co=so_co_re[,1], row.names=NULL)                           

aes_qual_re = ranef(aes_qual)[['LINKCEN2010']]
aes_qual_re$LINKCEN2010 = ifelse(nchar(row.names(aes_qual_re))==11, row.names(aes_qual_re), paste0('0', row.names(aes_qual_re)))
aes_qual_re = data.frame(aes_qual_re[2], aes_qual=aes_qual_re[,1], row.names=NULL)                           

safe_re = ranef(safe)[['LINKCEN2010']]
safe_re$LINKCEN2010 = ifelse(nchar(row.names(safe_re))==11, row.names(safe_re), paste0('0', row.names(safe_re)))
safe_re = data.frame(safe_re[2], safe=safe_re[,1], row.names=NULL)    


socenv08=merge(merge(so_co_re, aes_qual_re, by="LINKCEN2010", all=T), safe_re, by="LINKCEN2010", all=T)
socenv08$so_co = scale(socenv08$so_co)
socenv08$aes_qual = scale(socenv08$aes_qual)
socenv08$safe = scale(socenv08$safe)
socenv08$socenv = rowSums(socenv08[2:4], na.rm=T)
socenv08$year = 2008








#2010
lb = read.csv('S:/MESA/Methylation/Lauren_DNAmAge/data/lb10.csv', header=T)


so_co = lmer(so_co ~ MAGE + GENDER + as.factor(STFIPS10) + (1|LINKCEN2010) + (1|LINKCEN2010:hhidpn), data=lb)
aes_qual = lmer(aes_qual ~ MAGE + GENDER + as.factor(STFIPS10) + (1|LINKCEN2010) + (1|LINKCEN2010:hhidpn), data=lb)
safe = lmer(safe ~ MAGE + GENDER + as.factor(STFIPS10) + (1|LINKCEN2010), data=lb)


so_co_re = ranef(so_co)[['LINKCEN2010']]
so_co_re$LINKCEN2010 = ifelse(nchar(row.names(so_co_re))==11, row.names(so_co_re), paste0('0', row.names(so_co_re)))
so_co_re = data.frame(so_co_re[2], so_co=so_co_re[,1], row.names=NULL)                           

aes_qual_re = ranef(aes_qual)[['LINKCEN2010']]
aes_qual_re$LINKCEN2010 = ifelse(nchar(row.names(aes_qual_re))==11, row.names(aes_qual_re), paste0('0', row.names(aes_qual_re)))
aes_qual_re = data.frame(aes_qual_re[2], aes_qual=aes_qual_re[,1], row.names=NULL)                           

safe_re = ranef(safe)[['LINKCEN2010']]
safe_re$LINKCEN2010 = ifelse(nchar(row.names(safe_re))==11, row.names(safe_re), paste0('0', row.names(safe_re)))
safe_re = data.frame(safe_re[2], safe=safe_re[,1], row.names=NULL)    


socenv10=merge(merge(so_co_re, aes_qual_re, by="LINKCEN2010", all=T), safe_re, by="LINKCEN2010", all=T)
socenv10$so_co = scale(socenv10$so_co)
socenv10$aes_qual = scale(socenv10$aes_qual)
socenv10$safe = scale(socenv10$safe)
socenv10$socenv = rowSums(socenv10[2:4], na.rm=T)
socenv10$year = 2010






#2012
lb = read.csv('S:/MESA/Methylation/Lauren_DNAmAge/data/lb12.csv', header=T)


so_co = lmer(so_co ~ NAGE + GENDER + as.factor(STFIPS10) + (1|LINKCEN2010) + (1|LINKCEN2010:hhidpn), data=lb)
aes_qual = lmer(aes_qual ~ NAGE + GENDER + as.factor(STFIPS10) + (1|LINKCEN2010) + (1|LINKCEN2010:hhidpn), data=lb)
safe = lmer(safe ~ NAGE + GENDER + as.factor(STFIPS10) + (1|LINKCEN2010), data=lb)


so_co_re = ranef(so_co)[['LINKCEN2010']]
so_co_re$LINKCEN2010 = ifelse(nchar(row.names(so_co_re))==11, row.names(so_co_re), paste0('0', row.names(so_co_re)))
so_co_re = data.frame(so_co_re[2], so_co=so_co_re[,1], row.names=NULL)                           

aes_qual_re = ranef(aes_qual)[['LINKCEN2010']]
aes_qual_re$LINKCEN2010 = ifelse(nchar(row.names(aes_qual_re))==11, row.names(aes_qual_re), paste0('0', row.names(aes_qual_re)))
aes_qual_re = data.frame(aes_qual_re[2], aes_qual=aes_qual_re[,1], row.names=NULL)                           

safe_re = ranef(safe)[['LINKCEN2010']]
safe_re$LINKCEN2010 = ifelse(nchar(row.names(safe_re))==11, row.names(safe_re), paste0('0', row.names(safe_re)))
safe_re = data.frame(safe_re[2], safe=safe_re[,1], row.names=NULL)    


socenv12=merge(merge(so_co_re, aes_qual_re, by="LINKCEN2010", all=T), safe_re, by="LINKCEN2010", all=T)
socenv12$so_co = scale(socenv12$so_co)
socenv12$aes_qual = scale(socenv12$aes_qual)
socenv12$safe = scale(socenv12$safe)
socenv12$socenv = rowSums(socenv12[2:4], na.rm=T)
socenv12$year = 2012







#2014
lb = read.csv('S:/MESA/Methylation/Lauren_DNAmAge/data/lb14.csv', header=T)


so_co = lmer(so_co ~ OAGE + GENDER + as.factor(STFIPS10) + (1|LINKCEN2010) + (1|LINKCEN2010:hhidpn), data=lb)
aes_qual = lmer(aes_qual ~ OAGE + GENDER + as.factor(STFIPS10) + (1|LINKCEN2010) + (1|LINKCEN2010:hhidpn), data=lb)
safe = lmer(safe ~ OAGE + GENDER + as.factor(STFIPS10) + (1|LINKCEN2010), data=lb)


so_co_re = ranef(so_co)[['LINKCEN2010']]
so_co_re$LINKCEN2010 = ifelse(nchar(row.names(so_co_re))==11, row.names(so_co_re), paste0('0', row.names(so_co_re)))
so_co_re = data.frame(so_co_re[2], so_co=so_co_re[,1], row.names=NULL)                           

aes_qual_re = ranef(aes_qual)[['LINKCEN2010']]
aes_qual_re$LINKCEN2010 = ifelse(nchar(row.names(aes_qual_re))==11, row.names(aes_qual_re), paste0('0', row.names(aes_qual_re)))
aes_qual_re = data.frame(aes_qual_re[2], aes_qual=aes_qual_re[,1], row.names=NULL)                           

safe_re = ranef(safe)[['LINKCEN2010']]
safe_re$LINKCEN2010 = ifelse(nchar(row.names(safe_re))==11, row.names(safe_re), paste0('0', row.names(safe_re)))
safe_re = data.frame(safe_re[2], safe=safe_re[,1], row.names=NULL)    


socenv14=merge(merge(so_co_re, aes_qual_re, by="LINKCEN2010", all=T), safe_re, by="LINKCEN2010", all=T)
socenv14$so_co = scale(socenv14$so_co)
socenv14$aes_qual = scale(socenv14$aes_qual)
socenv14$safe = scale(socenv14$safe)
socenv14$socenv = rowSums(socenv14[2:4], na.rm=T)
socenv14$year = 2014







#2016
lb = read.csv('S:/MESA/Methylation/Lauren_DNAmAge/data/lb16.csv', header=T)


so_co = lmer(so_co ~ PAGE + GENDER + as.factor(STFIPS10) + (1|LINKCEN2010) + (1|LINKCEN2010:hhidpn), data=lb)
aes_qual = lmer(aes_qual ~ PAGE + GENDER + as.factor(STFIPS10) + (1|LINKCEN2010) + (1|LINKCEN2010:hhidpn), data=lb)
safe = lmer(safe ~ PAGE + GENDER + as.factor(STFIPS10) + (1|LINKCEN2010), data=lb)


so_co_re = ranef(so_co)[['LINKCEN2010']]
so_co_re$LINKCEN2010 = ifelse(nchar(row.names(so_co_re))==11, row.names(so_co_re), paste0('0', row.names(so_co_re)))
so_co_re = data.frame(so_co_re[2], so_co=so_co_re[,1], row.names=NULL)                           

aes_qual_re = ranef(aes_qual)[['LINKCEN2010']]
aes_qual_re$LINKCEN2010 = ifelse(nchar(row.names(aes_qual_re))==11, row.names(aes_qual_re), paste0('0', row.names(aes_qual_re)))
aes_qual_re = data.frame(aes_qual_re[2], aes_qual=aes_qual_re[,1], row.names=NULL)                           

safe_re = ranef(safe)[['LINKCEN2010']]
safe_re$LINKCEN2010 = ifelse(nchar(row.names(safe_re))==11, row.names(safe_re), paste0('0', row.names(safe_re)))
safe_re = data.frame(safe_re[2], safe=safe_re[,1], row.names=NULL)    


socenv16=merge(merge(so_co_re, aes_qual_re, by="LINKCEN2010", all=T), safe_re, by="LINKCEN2010", all=T)
socenv16$so_co = scale(socenv16$so_co)
socenv16$aes_qual = scale(socenv16$aes_qual)
socenv16$safe = scale(socenv16$safe)
socenv16$socenv = rowSums(socenv16[2:4], na.rm=T)
socenv16$year = 2016






socenv = rbind(socenv06, socenv08, socenv10, socenv12, socenv14, socenv16)


geo = read.sas7bdat('R:/Health_Retirement_Study/GeoCode/Detail/HRSxDetail2018v2/data/SAS/hrsxgeo18v8b_r.sas7bdat')
geo2 = merge(geo[c(1:3,14)], socenv, by.x=c('LINKCEN2010','YEAR'), by.y=c('LINKCEN2010','year'))

geo3 = reshape(geo2[-1], v.names=c("so_co", "aes_qual", "safe", "socenv"), idvar=c('HHID', 'PN'), timevar='YEAR', direction='wide', sep='_')
geo3$so_co = rowMeans(geo3[c(3,7,11,15,19,23)], na.rm=T)
geo3$aes_qual = rowMeans(geo3[c(4,8,12,16,20,24)], na.rm=T)
geo3$safe = rowMeans(geo3[c(5,9,13,17,21,25)], na.rm=T)
geo3$socenv = rowMeans(geo3[c(6,10,14,18,22,26)], na.rm=T)




write.csv(geo3[c(1,2,27:30)], 'S:/MESA/Methylation/Lauren_DNAmAge/data/HRS_socenv_06_16.csv', row.names = F, na='')




