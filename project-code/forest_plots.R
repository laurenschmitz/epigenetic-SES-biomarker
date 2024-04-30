
library(forestploter)
library(grid)
library(gridExtra)
library(stringr)
library(readxl)


############HRS##################
hrs.ases = read.csv('S:/MESA/Methylation/Lauren_DNAmAge/results/HRS_Enet/HRS_SESindex_mclock_models_table_raw.csv')

hrs.asesw = reshape(hrs.ases, direction='wide', idvar='SES', timevar='Model')

hrs.asesw$a <- paste(rep(" ", 26), collapse = " ")
hrs.asesw$` ` <- paste(rep(" ", 2), collapse = " ")
hrs.asesw$b <- paste(rep(" ", 26), collapse = " ")
hrs.asesw$`  ` <- paste(rep(" ", 2), collapse = " ")
hrs.asesw$c <- paste(rep(" ", 26), collapse = " ")
hrs.asesw$`   ` <- paste(rep(" ", 2), collapse = " ")
hrs.asesw$d <- paste(rep(" ", 26), collapse = " ")
hrs.asesw$`    ` <- paste(rep(" ", 2), collapse = " ")
hrs.asesw$e <- paste(rep(" ", 26), collapse = " ")


colnames(hrs.asesw)[c(1,32,34,36,38,40)] = c('' ,'Number of\nchronic conditions', 'Cardiometabolic\nconditions index', 'Self-reported\nhealth status', 'Mortality', 'Langa-Weir\ndementia')
hrs.asesw[1]=c('aSES-BIO', 'cSES-BIO')


tm = forest_theme(base_size = 16, colhead=list(fg_params=list(hjust=c(0,0.5,0.5,0.5,0.5,0.5,0.5), x=c(0,0.5,0.5,0.5,0.5,0.5,0.5))), core=list(bg_params=list(fill = c('white', 'white'))), xaxis_cex = .8, xlab_adjust=c('center','center','center'), legend_position = 'none')

p = forest(hrs.asesw[,c(1,32:36)], 
           est=list(hrs.asesw$estimate.r14conde.1,hrs.asesw$estimate.zMSI18.1,hrs.asesw$estimate.R14SHLT.1,hrs.asesw$estimate.r14conde.2,hrs.asesw$estimate.zMSI18.2,hrs.asesw$estimate.R14SHLT.2), 
           lower=list(hrs.asesw$lowercl.r14conde.1,hrs.asesw$lowercl.zMSI18.1,hrs.asesw$lowercl.R14SHLT.1,hrs.asesw$lowercl.r14conde.2,hrs.asesw$lowercl.zMSI18.2,hrs.asesw$lowercl.R14SHLT.2), 
           upper=list(hrs.asesw$uppercl.r14conde.1,hrs.asesw$uppercl.zMSI18.1,hrs.asesw$uppercl.R14SHLT.1,hrs.asesw$uppercl.r14conde.2,hrs.asesw$uppercl.zMSI18.2,hrs.asesw$uppercl.R14SHLT.2), 
           ci_column = c(2,4,6), 
           ref_line = c(1,0,0),
           xlim=list(c(0.98,1.18),c(-.02,.22),c(-.02,.22)),
           ticks_at = list(c(1.0,1.04,1.08,1.12,1.16),c(0.0,0.05,0.1,0.15,0.2),c(0.0,0.05,0.1,0.15,0.2)),
           xlab = c('IRR','Beta','Beta'),
           nudge_y = .4,
           theme = tm)

plot(p)



tm2 = forest_theme(base_size = 16, colhead=list(fg_params=list(hjust=c(0,0.5,0.5,0.5,0.5), x=c(0,0.5,0.5,0.5,0.5))), core=list(bg_params=list(fill = c('white', 'white'))), xaxis_cex = .8, xlab_adjust='center', legend_name = '', legend_value = c('Model 1  ','Model 2'), legend_position = 'bottom')

p2 = forest(hrs.asesw[,c(1,38:40)], 
           est=list(hrs.asesw$estimate.dead18.1,hrs.asesw$estimate.langawei.1,hrs.asesw$estimate.dead18.2,hrs.asesw$estimate.langawei.2), 
           lower=list(hrs.asesw$lowercl.dead18.1,hrs.asesw$lowercl.langawei.1,hrs.asesw$lowercl.dead18.2,hrs.asesw$lowercl.langawei.2), 
           upper=list(hrs.asesw$uppercl.dead18.1,hrs.asesw$uppercl.langawei.1,hrs.asesw$uppercl.dead18.2,hrs.asesw$uppercl.langawei.2), 
           ci_column = c(2,4), 
           ref_line = c(1,1),
           xlim=list(c(0.9,2.1),c(0.7,1.7)),
           ticks_at = list(c(1.0,1.25,1.5,1.75,2.0),c(0.8,1.0,1.2,1.4,1.6)),
           xlab = c('OR','OR'),
           nudge_y = .4,
           theme = tm2)

plot(p2)




png(filename = 'S:/MESA/Methylation/Lauren_DNAmAge/results/HRS_Enet/HRS_SESindex_mclock_models_forest_plot.png' , width = 850, height = 350,bg = "white",pointsize=24)
grid.arrange(p, p2)
dev.off()






#########MESA##############
mesa.ses = read.csv('S:/MESA/Methylation/Lauren_DNAmAge/results/HRS_Enet/MESA_SESindex_mclock_models_table_raw.csv')

mesa.sesw = reshape(mesa.ses, direction='wide', idvar='SES', timevar='Model')

mesa.sesw$a <- paste(rep(" ", 26), collapse = " ")
mesa.sesw$` ` <- paste(rep(" ", 2), collapse = " ")
mesa.sesw$b <- paste(rep(" ", 26), collapse = " ")
mesa.sesw$`  ` <- paste(rep(" ", 2), collapse = " ")
mesa.sesw$c <- paste(rep(" ", 26), collapse = " ")
mesa.sesw$`   ` <- paste(rep(" ", 2), collapse = " ")
mesa.sesw$d <- paste(rep(" ", 26), collapse = " ")


colnames(mesa.sesw)[c(1,26,28,30,32)] = c('' , 'Cardiometabolic\nconditions index', 'Self-reported\nhealth status', 'Mortality', 'ICD-based all cause\ndementia')
mesa.sesw[1]=c('aSES-BIO450', 'cSES-BIO450')


tm = forest_theme(base_size = 16, colhead=list(fg_params=list(hjust=c(0,0.5,0.5,0.5,0.5), x=c(0,0.5,0.5,0.5,0.5))), core=list(bg_params=list(fill = c('white', 'white'))), xaxis_cex = .8, xlab_adjust=c('center','center'), legend_position = 'none')

p = forest(mesa.sesw[,c(1,26:28)], 
           est=list(mesa.sesw$estimate.zMSI18.1,mesa.sesw$estimate.genhlth6.1,mesa.sesw$estimate.zMSI18.2,mesa.sesw$estimate.genhlth6.2), 
           lower=list(mesa.sesw$lowercl.zMSI18.1,mesa.sesw$lowercl.genhlth6.1,mesa.sesw$lowercl.zMSI18.2,mesa.sesw$lowercl.genhlth6.2), 
           upper=list(mesa.sesw$uppercl.zMSI18.1,mesa.sesw$uppercl.genhlth6.1,mesa.sesw$uppercl.zMSI18.2,mesa.sesw$uppercl.genhlth6.2), 
           ci_column = c(2,4), 
           ref_line = c(0,0),
           xlim=list(c(-.07,.22),c(-.07,.22)),
           ticks_at = list(c(-0.05,0.0,0.05,0.1,0.15,0.2),c(-0.05,0.0,0.05,0.1,0.15,0.2)),
           xlab = c('Beta','Beta'),
           nudge_y = .4,
           theme = tm)

plot(p)



tm2 = forest_theme(base_size = 16, colhead=list(fg_params=list(hjust=c(0,0.5,0.5,0.5,0.5), x=c(0,0.5,0.5,0.5,0.5))), core=list(bg_params=list(fill = c('white', 'white'))), xaxis_cex = .8, xlab_adjust='center', legend_name = '', legend_value = c('Model 1  ','Model 2'), legend_position = 'bottom')

p2 = forest(mesa.sesw[,c(1,30:32)], 
            est=list(mesa.sesw$estimate.dth.1,mesa.sesw$estimate.icddemen.1,mesa.sesw$estimate.dth.2,mesa.sesw$estimate.icddemen.2), 
            lower=list(mesa.sesw$lowercl.dth.1,mesa.sesw$lowercl.icddemen.1,mesa.sesw$lowercl.dth.2,mesa.sesw$lowercl.icddemen.2), 
            upper=list(mesa.sesw$uppercl.dth.1,mesa.sesw$uppercl.icddemen.1,mesa.sesw$uppercl.dth.2,mesa.sesw$uppercl.icddemen.2), 
            ci_column = c(2,4), 
            ref_line = c(1,1),
            xlim=list(c(0.69,1.71),c(0.64,1.86)),
            ticks_at = list(c(0.8,1.0,1.2,1.4,1.6),c(0.75,1.0,1.25,1.5,1.75)),
            xlab = c('OR','OR'),
            nudge_y = .4,
            theme = tm2)

plot(p2)




png(filename = 'S:/MESA/Methylation/Lauren_DNAmAge/results/HRS_Enet/MESA_SESindex_mclock_models_forest_plot.png' , width = 850, height = 350,bg = "white",pointsize=24)
grid.arrange(p, p2)
dev.off()



