library(qqman)
library(bacon)




#race = c('M1', 'M2', 'M3', 'F1', 'F2', 'F3')
race = c('ALL', 'EA', 'AA', 'HA', 'F', 'M')
ses = c('COLLEDUC', 'EDYEARS', 'HSEDUC', 'PEDUC', 'HHINC', 'NDS', 'POV', 'socenv', 'peducyrs', 'rafeduc2', 'rameduc2', 'occ_prest', 'lnhhincx', 'lnwealth2x', 'SESindex', 'cSES')



#Put EPIC results together

for (i in 1:14){
  for (j in 1:6){

      final = data.frame(CPG=character(0), CHR=integer(0), POS=integer(0), N=integer(0), BETA=numeric(0), SE=numeric(0), P=numeric(0))
      for (n in 1:9){
        file = paste0('S:/MESA/Methylation/Lauren_DNAmAge/results/HRS Enet/EWAS/HRS_', ses[i], '_', race[j], '_Model2_', n, '.csv')
        fileout = paste0('S:/MESA/Methylation/Lauren_DNAmAge/results/HRS Enet/EWAS/HRS_', ses[i], '_', race[j], '_Model2.csv')
        qqout = paste0('S:/MESA/Methylation/Lauren_DNAmAge/results/HRS Enet/EWAS/HRS_', ses[i], '_', race[j], '_Model2_QQplot.png')
        qqout.bacon = paste0('S:/MESA/Methylation/Lauren_DNAmAge/results/HRS Enet/EWAS/HRS_', ses[i], '_', race[j], '_Model2_QQplot_bacon.png')
        
        results = read.csv(file, header=T)
        final = rbind(final, results)
      } 
      
      trimtest = data.frame(trim=numeric(0), lambda=numeric(0))
      for (k in 1:30){
        set.seed(2736459)
        trimit=(.97+.001*(k-1))
        z.bacon = bacon(final$t, trim=trimit)
        final$P.bacon = pval(z.bacon)
        chisq2=qchisq(1-final$P.bacon,1)
        lambda.bacon=round((median(chisq2)/qchisq(0.5,1)), digits=3)
        trimtest[k,1]=trimit
        trimtest[k,2]=lambda.bacon
      }
      
      set.seed(2736459)
      trimit=0.99
      z.bacon = bacon(final$t, trim=trimit)
      final$P.bacon = pval(z.bacon)
      
      chisq=qchisq(1-final$P,1)
      lambda=round((median(chisq)/qchisq(0.5,1)), digits=3)
      
      chisq2=qchisq(1-final$P.bacon,1)
      lambda.bacon=round((median(chisq2)/qchisq(0.5,1)), digits=3)
      
      final$P.FDR = p.adjust(final$P)
      final$P.bacon.FDR = p.adjust(final$P.bacon)
      
      write.csv(final, fileout, row.names=F)
      
      title = paste0('QQ Plot of HRS ', ses[i], ' EWAS Results - ', race[j], ' Model 2\nLambda=', lambda)
      title.bacon = paste0('Bacon (trim=', trimit, ') Adjusted QQ Plot of HRS ', ses[i], ' EWAS Results - ', race[j], ' Model 2\nLambda=', lambda.bacon)
      
      png(filename = qqout , width = 1500, height = 1500,bg = "white",pointsize=24)
      par(mar=c(5,6,4,2)+.1)
      qq(final$P, main=title, cex.main=2, cex.lab=2, cex.axis=2)
      dev.off()
      
      png(filename = qqout.bacon , width = 2000, height = 1500,bg = "white",pointsize=24)
      par(mar=c(5,10,4,10)+.1)
      qq(final$P.bacon, main=title.bacon, cex.main=2, cex.lab=2, cex.axis=2)
      dev.off()
      
    
  } 
}








#For Manhattan plots of SESindex and cSES
final$CHRn = ifelse(final$CHR=='X', 23, ifelse(final$CHR=='Y', 24, as.numeric(final$CHR)))
title.man = paste0('Manhattan Plot of HRS ', ses[i], ' EWAS Results')

enetsites=read.csv('S:/MESA/Methylation/Lauren_DNAmAge/results/HRS Enet/HRS_Enet_XValid_SESindex_betas_lambdamin_final.csv', header=T)[-1,1]

png(filename = paste0('S:/MESA/Methylation/Lauren_DNAmAge/results/HRS Enet/EWAS/HRS_', ses[i], '_', race[j], '_Model2_Manhattan.png') , width = 2200, height = 1500,bg = "white",pointsize=30)
par(mar=c(5,6,4,2)+.1)
manhattan(final, chr='CHRn', bp='POS', snp='CPG', main=title.man, cex=2, cex.main=2, cex.lab=2, cex.axis=2, highlight=enetsites, chrlabs=c(1:22,'X','Y'), suggestiveline=-log10(1.2e-07), genomewideline=-log10(6e-08))
dev.off()


enetsites2=read.csv('S:/MESA/Methylation/Lauren_DNAmAge/results/HRS Enet/HRS_Enet_XValid_cSES_betas_lambdamin_final.csv', header=T)[-1,1]

png(filename = paste0('S:/MESA/Methylation/Lauren_DNAmAge/results/HRS Enet/EWAS/HRS_', ses[i], '_', race[j], '_Model2_Manhattan.png') , width = 2200, height = 1500,bg = "white",pointsize=30)
par(mar=c(5,6,4,2)+.1)
manhattan(final, chr='CHRn', bp='POS', snp='CPG', main=title.man, cex=2, cex.main=2, cex.lab=2, cex.axis=2, highlight=enetsites2, chrlabs=c(1:22,'X','Y'), suggestiveline=-log10(1.2e-07), genomewideline=-log10(6e-08))
dev.off()

















