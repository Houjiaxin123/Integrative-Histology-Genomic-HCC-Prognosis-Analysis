library(glmnet)
library(rms)
library(VIM)
library(survival)
library(Hmisc)
library(timeROC)


load('trainset_eigen_survival_231.RData')

eigengene = eigengene[, -9]

# survival$OS_STATUS[survival$OS_STATUS == 'died']<- 1
# survival$OS_STATUS[survival$OS_STATUS == 'alive']<- 0
# survival$OS_STATUS = as.numeric(survival$OS_STATUS)

x <- data.matrix(eigengene)
y <- data.matrix(Surv(survival$OS_MONTHS, survival$OS_STATUS))

#LASSO
set.seed(123)
fit <- glmnet(x, y, family = "cox", alpha = 1, nlambda = 100)
png("LASSO-1.png", width = 1000, height = 800)
plot(fit,label=T)
dev.off()
png('LASSO-2.png', width = 1000, height = 800)
plot(fit,xvar="lambda",label=T)
dev.off()

fitcv <- cv.glmnet(x, y, family="cox", alpha=1, nlambda = 100, nfolds=5)
png("LASSO-lambda.min.png",width = 1000,height = 800)
plot(fitcv)
coef(fitcv, s="lambda.min")
dev.off()

dd<-datadist(eigengene)
options(datadist='dd')
coxm = cph(Surv(survival$OS_MONTHS,survival$OS_STATUS) ~ MEturquoise+MEyellow+MEred+MEpink+MEgreen,data = eigengene,x=T, y=T, surv=T)
cox.zph(coxm)

#C-index
f<-coxph(Surv(survival$OS_MONTHS,survival$OS_STATUS)~ MEturquoise+MEyellow+MEred+MEpink+MEgreen,data = eigengene)
sum.surv<-summary(f)
c_index<-sum.surv$concordance
