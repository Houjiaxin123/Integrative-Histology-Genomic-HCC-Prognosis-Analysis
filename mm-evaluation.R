rm(list = ls())
options(stringsAsFactors = F)
setwd('C:/Users/12sigma/Desktop/paper-revised-exp')


library(glmnet)
library(rms)
library(VIM)
library(survival)
library(Hmisc)
library(timeROC)
library(ggplot2)
library(survminer)
library(survcomp)
library(ggDCA)
library(pec)

load('trainset_all_modality.RData')
load('testset_all_modality.RData')

##c-index
train_model_mm<-coxph(Surv(OS_MONTHS,OS_STATUS)~ DL_RISK_SCORE+DCAF13+ELAC2+ZNF320+KIF18B+FERMT3+SEX+AGE, data = train_all_modal)
train_model_patho<-coxph(Surv(OS_MONTHS,OS_STATUS)~ DL_RISK_SCORE+SEX+AGE, data = train_all_modal)
train_model_gene<-coxph(Surv(OS_MONTHS,OS_STATUS)~ DCAF13+ELAC2+ZNF320+KIF18B+FERMT3+SEX+AGE, data = train_all_modal)
c_index_mm = concordance.index(predict(train_model_mm),surv.time = train_all_modal$OS_MONTHS, surv.event = train_all_modal$OS_STATUS,method = "noether")
c_index_patho = concordance.index(predict(train_model_patho),surv.time = train_all_modal$OS_MONTHS, surv.event = train_all_modal$OS_STATUS,method = "noether") 
c_index_gene = concordance.index(predict(train_model_gene),surv.time = train_all_modal$OS_MONTHS, surv.event = train_all_modal$OS_STATUS,method = "noether") 

# sum.surv.mm<-summary(train_model_mm)
# c_index_mm<-sum.surv.mm$concordance
# sum.surv.patho<-summary(train_model_patho)
# c_index_patho<-sum.surv.patho$concordance
# sum.surv.gene<-summary(train_model_gene)
# c_index_gene<-sum.surv.gene$concordance

##modal test
pre_mm<-predict(train_model_mm, newdata=test_all_modal)
f_test_mm <- coxph(Surv(OS_MONTHS, OS_STATUS)~pre_mm, data=test_all_modal)
c_index_test_mm = concordance.index(predict(f_test_mm),surv.time = test_all_modal$OS_MONTHS, surv.event = test_all_modal$OS_STATUS,method = "noether")

# sum.surv.test.mm = summary(f_test_mm)
# c_index_test_mm = sum.surv.test.mm$concordance

pre_patho<-predict(train_model_patho, newdata=test_all_modal)
f_test_patho <- coxph(Surv(OS_MONTHS, OS_STATUS)~pre_patho, data=test_all_modal)
c_index_test_patho = concordance.index(predict(f_test_patho),surv.time = test_all_modal$OS_MONTHS, surv.event = test_all_modal$OS_STATUS,method = "noether")

# sum.surv.test.patho = summary(f_test_patho)
# c_index_test_patho = sum.surv.test.patho$concordance

pre_gene<-predict(train_model_gene, newdata=test_all_modal)
f_test_gene <- coxph(Surv(OS_MONTHS, OS_STATUS)~pre_gene, data=test_all_modal)
c_index_test_gene = concordance.index(predict(f_test_gene),surv.time = test_all_modal$OS_MONTHS, surv.event = test_all_modal$OS_STATUS,method = "noether")

# sum.surv.test.gene = summary(f_test_gene)
# c_index_test_gene = sum.surv.test.gene$concordance

##cindex compare
p_train_1 = cindex.comp(c_index_mm, c_index_gene)
p_train_2 = cindex.comp(c_index_mm, c_index_patho)

p_test_1 = cindex.comp(c_index_test_mm, c_index_test_gene)
p_test_2 = cindex.comp(c_index_test_mm, c_index_test_patho)

##Forest plot
png('hazard-ratio-forest-new.png', width = 800, height = 800)
ggforest(train_model_mm, 
         main = 'Hazard Ratio',
         data = test_all_modal,
         cpositions = c(0.05, 0.20, 0.35),
         fontsize = 1.2, 
         refLabel = 'reference',
         noDigits = 3)
dev.off()


##Nomogram (internal)
dd<-datadist(train_all_modal) 
options(datadist='dd') 
coxm_train = cph(Surv(OS_MONTHS,OS_STATUS)~ DL_RISK_SCORE+DCAF13+ELAC2+ZNF320+KIF18B+FERMT3+SEX+AGE,data = train_all_modal,x=T, y=T, surv=T)
cox.zph(coxm_train)
surv <- Survival(coxm_train)
surv1 <- function(x)surv(1*12,lp=x)
surv2 <- function(x)surv(3*12,lp=x)

png('nomogram-trainset.png', width = 1000, height = 600)
plot(nomogram(coxm_train,
              fun=list(surv1,surv2),
              lp= F,
              funlabel=c('1-year survival','3-year survival'),
              maxscale=100,
              fun.at=c('0.95','0.9','0.8','0.7','0.6','0.4','0.2','0.1')),
     cex.axis=0.9, label.every=1)
dev.off()

##calibration curve
cal1 = calibrate(coxm_train, cmethod = 'KM',method = 'boot' , u = 12*1 , m = 77, B = 100)
cal2 = calibrate(coxm_train, cmethod = 'KM',method = 'boot' , u = 12*3 , m = 77, B = 100)

par(mar = c(8,5,3,2), cex = 1.0)
png('C-index-1year-3year-trainset.png', width = 800, height = 600)
plot(cal1, lwd = 2, lty = 2, 
     add = F,
     subtitles = F, 
     errbar.col = 'salmon',
     xlim = c(0,1), 
     ylim = c(0,1),
     xlab = "Nomogram-Predicted Probability", 
     ylab = "Actual OS (proportion)",
     col = 'salmon')
plot(cal2, lwd = 2, lty = 2,
     add = T,
     subtitles = F,
     errbar.col = 'deepskyblue',
     xlim = c(0,1),
     ylim = c(0,1),
     xlab = "Nomogram-Predicted Probability",
     ylab = "Actual OS (proportion)",
     col = 'deepskyblue')
legend('bottomright',c(paste('1-year curve'),
                       paste('3-year curve')),
       col = c("salmon","deepskyblue"),lwd = 2.5, cex = 1.5)
abline(0, 1, lty = 3, lwd = 2, col = 'grey')
dev.off()

##Nnomogram (External)
pre_mm_score<-predict(coxm_train, newdata=test_all_modal)
coxm_test = cph(Surv(OS_MONTHS,OS_STATUS)~ pre_mm_score,data = test_all_modal,x=T, y=T, surv=T)
surv <- Survival(coxm_test)
surv1 <- function(x)surv(1*12,lp=x)
surv2 <- function(x)surv(3*12,lp=x)


##calibration curve
cal4 = calibrate(coxm_test, cmethod = 'KM',method = 'boot', u = 12*1, m = 38, B = 10)
cal5 = calibrate(coxm_test, cmethod = 'KM',method = 'boot', u = 12*3, m = 38, B = 10)

par(mar = c(8,5,3,2), cex = 1.0)
png('C-index-1year-3year-testset.png', width = 800, height = 600)
plot(cal4, lwd = 2, lty = 2, 
     add = F,
     subtitles = F, 
     errbar.col = 'salmon',
     xlim = c(0,1), 
     ylim = c(0,1),
     xlab = "Nomogram-Predicted Probability", 
     ylab = "Actual OS (proportion)",
     col = 'salmon')
plot(cal5, lwd = 2, lty = 2,
     add = T,
     subtitles = F,
     errbar.col = 'deepskyblue',
     xlim = c(0,1),
     ylim = c(0,1),
     xlab = "Nomogram-Predicted Probability",
     ylab = "Actual OS (proportion)",
     col = 'deepskyblue')
legend('bottomright',c(paste('1-year curve'),
                       paste('3-year curve')),
       col = c("salmon","deepskyblue"),lwd = 2.5, cex = 1.5)
abline(0, 1, lty = 3, lwd = 2, col = 'grey')
dev.off()

##Risk score calculation
cox_train_model_patho = step(train_model_patho, direction = 'both')
cox_train_model_gene = step(train_model_gene, direction = 'both')
cox_train_model_mm = step(train_model_mm, direction = 'both')

risk_score_patho = predict(cox_train_model_patho, type = 'risk', newdata = test_all_modal)
risk_score_gene = predict(cox_train_model_gene, type = 'risk', newdata = test_all_modal)
risk_score_mm = predict(cox_train_model_mm, type = 'risk', newdata = test_all_modal)

risk_level_patho = as.vector(ifelse(risk_score_patho > median(risk_score_patho),'High','Low'))
write.table(cbind(id=rownames(test_all_modal), cbind(test_all_modal[,4:5],risk_score_patho,risk_level_patho)), "risk_score_patho.txt",sep='\t',quote=F,row.names=F)

risk_level_gene = as.vector(ifelse(risk_score_gene > median(risk_score_gene),'High','Low'))
write.table(cbind(id=rownames(test_all_modal), cbind(test_all_modal[,4:5],risk_score_gene,risk_level_gene)), "risk_score_gene.txt",sep='\t',quote=F,row.names=F)

risk_level_mm = as.vector(ifelse(risk_score_mm > median(risk_score_mm),'High','Low'))
write.table(cbind(id=rownames(test_all_modal), cbind(test_all_modal[,4:5],risk_score_mm,risk_level_mm)), "risk_score_mm.txt",sep='\t',quote=F,row.names=F)

##ROC plot (test set)
risk_patho=read.table("risk_score_patho.txt",header = T)
risk_gene=read.table("risk_score_gene.txt",header = T)
risk_mm=read.table("risk_score_mm.txt",header = T)

predict_1_year = 1*12
predict_3_year = 3*12

ROC_patho =timeROC(T=risk_patho$OS_MONTHS,delta = risk_patho$OS_STATUS,marker = risk_patho$risk_score, cause=1,
                   weighting = 'marginal',times = c(predict_1_year,predict_3_year),ROC = T)
ROC_gene =timeROC(T=risk_gene$OS_MONTHS,delta = risk_gene$OS_STATUS,marker = risk_gene$risk_score, cause=1,
                  weighting = 'marginal',times = c(predict_1_year,predict_3_year),ROC = T)
ROC_mm =timeROC(T=risk_mm$OS_MONTHS,delta = risk_mm$OS_STATUS,marker = risk_mm$risk_score, cause=1,
                weighting = 'marginal',times = c(predict_1_year,predict_3_year),ROC = T)
png('ROC-1year-testset.png', width = 600, height = 600)
plot(ROC_patho,time = predict_1_year,title = F,lwd=2.5,col = 'salmon')
plot(ROC_gene,time = predict_1_year,title = F,lwd=2.5,col = 'deepskyblue',add = T)
plot(ROC_mm,time = predict_1_year,title = F,lwd=2.5,col = 'limegreen',add = T)
legend('bottomright',c(paste('1-year pathology AUC =',round(ROC_patho$AUC[1],3)),
                       paste('1-year gene AUC =',round(ROC_gene$AUC[1],3)),
                       paste('1-year multi-modality AUC =',round(ROC_mm$AUC[1],3))),
       col = c("salmon","deepskyblue","limegreen"),lwd = 2.5, cex = 1.5)
dev.off()

png('ROC-3year-testset.png', width = 600, height = 600)
plot(ROC_patho,time = predict_3_year,title = F,lwd=2.5, col = 'salmon')
plot(ROC_gene,time = predict_3_year,title = F,lwd=2.5,col = 'deepskyblue',add = T)
plot(ROC_mm,time = predict_3_year,title = F,lwd=2.5,col = 'limegreen',add = T)
legend('bottomright',c(paste('3-year pathology AUC =',round(ROC_patho$AUC[2],3)),
                       paste('3-year gene AUC =',round(ROC_gene$AUC[2],3)),
                       paste('3-year multi-modality AUC =',round(ROC_mm$AUC[2],3))),
       col = c("salmon","deepskyblue","limegreen"),lwd = 2.5, cex = 1.5)
dev.off()

##Risk score calculation (training set)
risk_score_patho_trainset = predict(cox_train_model_patho, type = 'risk', newdata = train_all_modal)
risk_score_gene_trainset = predict(cox_train_model_gene, type = 'risk', newdata = train_all_modal)
risk_score_mm_trainset = predict(cox_train_model_mm, type = 'risk', newdata = train_all_modal)

risk_level_patho_trainset = as.vector(ifelse(risk_score_patho_trainset > median(risk_score_patho_trainset),'High','Low'))
write.table(cbind(id=rownames(train_all_modal), cbind(train_all_modal[,4:5],risk_score_patho_trainset,risk_level_patho_trainset)), "risk_score_patho_trainset.txt",sep='\t',quote=F,row.names=F)

risk_level_gene_trainset = as.vector(ifelse(risk_score_gene_trainset > median(risk_score_gene_trainset),'High','Low'))
write.table(cbind(id=rownames(train_all_modal), cbind(train_all_modal[,4:5],risk_score_gene_trainset,risk_level_gene_trainset)), "risk_score_gene_trainset.txt",sep='\t',quote=F,row.names=F)

risk_level_mm_trainset = as.vector(ifelse(risk_score_mm_trainset > median(risk_score_mm_trainset),'High','Low'))
write.table(cbind(id=rownames(train_all_modal), cbind(train_all_modal[,4:5],risk_score_mm_trainset,risk_level_mm_trainset)), "risk_score_mm_trainset.txt",sep='\t',quote=F,row.names=F)

##ROC plot (training set)
risk_patho_trainset=read.table("risk_score_patho_trainset.txt",header = T)
risk_gene_trainset=read.table("risk_score_gene_trainset.txt",header = T)
risk_mm_trainset=read.table("risk_score_mm_trainset.txt",header = T)

ROC_patho_train =timeROC(T=risk_patho_trainset$OS_MONTHS,delta = risk_patho_trainset$OS_STATUS,marker = risk_patho_trainset$risk_score, cause=1,
                         weighting = 'marginal',times = c(predict_1_year,predict_3_year),ROC = T)
ROC_gene_train =timeROC(T=risk_gene_trainset$OS_MONTHS,delta = risk_gene_trainset$OS_STATUS,marker = risk_gene_trainset$risk_score, cause=1,
                        weighting = 'marginal',times = c(predict_1_year,predict_3_year),ROC = T)
ROC_mm_train =timeROC(T=risk_mm_trainset$OS_MONTHS,delta = risk_mm_trainset$OS_STATUS,marker = risk_mm_trainset$risk_score, cause=1,
                      weighting = 'marginal',times = c(predict_1_year,predict_3_year),ROC = T)
png('ROC-1year-trainset.png', width = 600, height = 600)
plot(ROC_patho_train,time = predict_1_year,title = F,lwd=2.5,col = 'salmon')
plot(ROC_gene_train,time = predict_1_year,title = F,lwd=2.5,col = 'deepskyblue',add = T)
plot(ROC_mm_train,time = predict_1_year,title = F,lwd=2.5,col = 'limegreen',add = T)
legend('bottomright',c(paste('1-year pathology AUC =',round(ROC_patho_train$AUC[1],3)),
                       paste('1-year gene AUC =',round(ROC_gene_train$AUC[1],3)),
                       paste('1-year multi-modality AUC =',round(ROC_mm_train$AUC[1],3))),
       col = c("salmon","deepskyblue","limegreen"),lwd = 2.5, cex = 1.5)
dev.off()

png('ROC-3year-trainset.png', width = 600, height = 600)
plot(ROC_patho_train,time = predict_3_year,title = F,lwd=2.5, col = 'salmon')
plot(ROC_gene_train,time = predict_3_year,title = F,lwd=2.5,col = 'deepskyblue',add = T)
plot(ROC_mm_train,time = predict_3_year,title = F,lwd=2.5,col = 'limegreen',add = T)
legend('bottomright',c(paste('3-year pathology AUC =',round(ROC_patho_train$AUC[2],3)),
                       paste('3-year gene AUC =',round(ROC_gene_train$AUC[2],3)),
                       paste('3-year multi-modality AUC =',round(ROC_mm_train$AUC[2],3))),
       col = c("salmon","deepskyblue","limegreen"),lwd = 2.5, cex = 1.5)
dev.off()


##KM curve
kms_patho = survfit(Surv(risk_patho$OS_MONTHS,risk_patho$OS_STATUS)~risk_patho$risk_level_patho,data = risk_patho)
kms_gene = survfit(Surv(risk_gene$OS_MONTHS,risk_gene$OS_STATUS)~risk_gene$risk_level_gene,data = risk_gene)
kms_mm = survfit(Surv(risk_mm$OS_MONTHS,risk_mm$OS_STATUS)~risk_mm$risk_level_mm,data = risk_mm)

png('survival-curve-KM-test-patho.png', width = 600, height = 500)
ggsurvplot(kms_patho,
           pval = TRUE,
           pval.method = TRUE,
           conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = 'strata',
           xlab = "survival time in months",
           ylab = "survival probabilities",
           break.x.by = 12,
           surv.median.line = 'hv',
           ggtheme = theme_bw(),
           palette = c('salmon', 'deepskyblue'),
           legend.labs = c('High risk', 'Low risk'),
           legand.title = 'Risk',
           size = 1)
dev.off()

png('survival-curve-KM-test-gene.png', width = 600, height = 500)
ggsurvplot(kms_gene,
           pval = TRUE,
           pval.method = TRUE,
           conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = 'strata',
           xlab = "survival time in months",
           ylab = "survival probabilities",
           break.x.by = 12,
           surv.median.line = 'hv',
           ggtheme = theme_bw(),
           palette = c('salmon', 'deepskyblue'),
           legend.labs = c('High risk', 'Low risk'),
           legand.title = 'Risk',
           size = 1)
dev.off()

png('survival-curve-KM-test-mm.png', width = 600, height = 500)
ggsurvplot(kms_mm,
           pval = TRUE,
           pval.method = TRUE,
           conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = 'strata',
           xlab = "survival time in months",
           ylab = "survival probabilities",
           break.x.by = 12,
           surv.median.line = 'hv',
           ggtheme = theme_bw(),
           palette = c('salmon', 'deepskyblue'),
           legend.labs = c('High risk', 'Low risk'),
           legand.title = 'Risk',
           size = 1)
dev.off()

##DCA (training set)
png('decision-curve-train.png', width = 800, height = 600)
fig = dca(train_model_patho, train_model_gene, train_model_mm,
          new.data = train_all_modal)
ggplot(fig,
       linetype = F, 
       lwd = 1,
       col = c("salmon","deepskyblue","limegreen", "grey", "black"),
       labs = c('Patho model', 'Gene model', 'Multi-modality model', 'All', 'None'))+xlim(0,0.5)+theme(legend.position = c(0.85,0.9), legend.text = element_text(size = 20))
dev.off()


##DCA (test set)
test_model_patho = train_model_patho
test_model_gene = train_model_gene
test_model_mm = train_model_mm
png('decision-curve-test.png', width = 800, height = 600)
fig = dca(test_model_patho, test_model_gene, test_model_mm,
          new.data = test_all_modal)
ggplot(fig,
       linetype = F, 
       lwd = 1,
       col = c("salmon","deepskyblue","limegreen", "grey", "black"),
       labs = c('Patho model', 'Gene model', 'Multi-modality model', 'All', 'None'))+xlim(0,0.5)+theme(legend.position = c(0.85,0.9), legend.text = element_text(size = 20))
dev.off()



