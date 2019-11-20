## Chapter 9

## ------------------------------------
## 例9-1
## ------------------------------------
rm(list=ls())
## 读入数据
X1=read.csv("eg9-1.csv")
X1

## ----------
## (1)预分析
## 生成列联表
X<-matrix(c(53,434,111,38,108,48),3,2,
          dimnames=list(A=c('高','中','低'),B=c('满意','不满意')))
t<-as.table(X)
t

## 频率列联表
prop.table(t)

## 计算期望频数
E=matrix(c(1:6),3,dimnames=list(A=c('高','中','低'),B=c('满意','不满意')))
for(i in 1:3){
  for(j in 1:2){
    E[i,j]=rowSums(t)[i]*colSums(t)[j]/sum(t)
  }
}
E

## 似然比统计量
lhr=2*sum(t*log(t/E))
p.value.lhr=pchisq(lhr,2,lower.tail = F)
## perason统计量
psn=sum((t-E)^2/E)
p.value.psn=pchisq(psn,2,lower.tail = F)
STAT<-list("lhr"=lhr,"p.value.lhr"=p.value.lhr,"psn"=psn,"p.value.psn"=p.value.psn)
STAT


## ----------
## (2)非饱和模型的分析
library(MASS)
X.m1<-loglm(~A+B,data=t)
X.m1$param
X.m1$lrt
X.m1$pearson

y=X1$频数
x1=X1$收入情况
x2=X1$满意情况
log.glm<-glm(y~x1*x2,family=poisson(link=log),data=X1)
log.glm$null.deviance

## ---仿照SPSS输出K-way and Higher-Order Effects---##
LHR<-c(log.glm$null.deviance+lhr,lhr,log.glm$null.deviance,lhr)
df<-c(5,2,3,2)
P.VALUE=c()
for(i in 1:4)
{
  p=pchisq(LHR[i],df[i],lower.tail = F)
  P.VALUE<-c(P.VALUE,p)
}
K.way<-data.frame(df,LHR,P.VALUE)
rownames(K.way)<-c("K-way & Higher Order Effects 1","K-way & Higher Order Effects 2",
                   "K-Way Effects 1","K-way Effects 2")
K.way


## ----------
## (3)饱和模型的分析
X.m2<-loglm(~A*B,data=t)
X.m2
X.m2$param



## ------------------------------------
## 例9-2
## ------------------------------------
## 读入数据
Y=read.csv("eg9-2.csv")
Y

## ----------
## (1)预分析
## 画图观察
attach(Y)
library(ggplot2)
ggplot(Y,aes(x=x, y=p))+
  geom_point(size=2)+
  labs(x = "x", y = "p") 

## ----------
## (2)普通最小二乘回归
y=Y$p
x=Y$x
Y.m1<-lm(y~x)
summary(Y.m1)
## 进行预测
pre<-predict(Y.m1,data.frame(x=8))
p_hat=exp(pre)/(1+exp(pre))
p_hat

## ----------
## (3)加权最小二乘回归
## 诊断异方差
res=abs(residuals(Y.m1))
library(psych)
corr.test(data.frame(x,res),use = "complete")
## WLS回归
w=Y$WLS
Y.m2<-lm(y~x,weights = w)
summary(Y.m2)
## 进行预测
pre<-predict(Y.m2,data.frame(x=8))
p_hat<-exp(pre)/(1+exp(pre))
p_hat



## ------------------------------------
## 例9-3
## ------------------------------------
## 读入数据
Z=read.csv("eg9-3.csv")
head(Z)
tail(Z)

## ----------
## (1)调用glm()函数进行估计
## 估计Logistics模型
y=Z$Y
sex=Z$SEX
age=Z$AGE
x2=Z$X2
Z.m1<-glm(y~sex+age+x2,family = binomial(link = logit))
summary(Z.m1)

## 逐步回归
Z.m2<-step(Z.m1)
summary(Z.m2)
exp(Z.m2$coefficients)


## ----------
## (2)多重共线性的检验
library(DAAG)
vif(Z.m2)


## ----------
## (3)计算判别精度
logit=predict(Z.m2,data.frame(sex=Z$SEX,age=Z$AGE))
probibility=exp(logit)/(1+exp(logit))

newY=c()
for(i in 1:length(probibility)){
  if(probibility[i]<0.5)
    newY=c(newY,0)
  if(probibility[i]>0.5)
    newY=c(newY,1)
  if(probibility[i]==0.5)
    newY=c(newY,NA)
}

new.Z=cbind(Z,probibility,newY)
new.Z
table(y,newY)
sum(diag(prop.table(table(y,newY))))