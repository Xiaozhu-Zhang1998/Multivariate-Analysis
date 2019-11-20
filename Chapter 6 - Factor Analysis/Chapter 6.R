## Chapter 6

## ------------------------------------
## 例6-1
## ------------------------------------
rm(list=ls())

## 相关系数及其检验
## 读入数据
X=read.csv("eg6_1.csv",header=T)
rownames(X)=X[,1]
X=X[,-1]
X1=X[,c(3,5:8)]
p=ncol(X1)
head(X1)

## ----------
## (1)计算变量之间的相关性
library(psych)
corr.test(X1,use = "complete")


## ----------
## (2)主成分法：求解因子载荷
## 先进行主成分分析
source("Pca.R")
X.pca=Pca(X1,cor=T,1)
k=X.pca$`Number of Chosen Components`
X.pca[-7]

## 特殊度
1-X.pca$"Communalities"

## 计算的是“原始变量标准化的系数矩阵”
A=c()
for (i in 1:k){
  a=X.pca$`Coefficient of Components`[,i]/sqrt(X.pca$`Variance of Component`[i])
  A=cbind(A,a)
}
colnames(A)=paste("factor",c(1:k),sep="")
A


## ----------
## (3)因子旋转
## Step1：正交旋转
X.var=principal((X1),nfactors = 3,rotate="varimax",scores=T)
X.var
X.var$weights  #Component Score Coefficient Matrix
head(X.var$scores)

## Step2：斜交旋转
X.pro=principal((X1),nfactors = 3,rotate="promax",scores=T)
X.pro
X.pro$loadings[1:p,1:3]%*%X.pro$r.scores #计算Structure
X.pro$weights
head(X.pro$scores)

## Step3：验证
## 对因子的"标准化"进行验证
## 求均值
apply(X.pro$scores,2,mean)
## 求标准差
apply(X.pro$scores,2,sd)


## ----------
## (4)可视化：散点图
library(ggplot2)
ggplot(X1,aes(x=X.pro$scores[,1], y = X.pro$scores[,2]))+
  geom_point()+
  labs(x = "RC1", y = "RC2")



## ------------------------------------
## 例6-2
## ------------------------------------
## (1)函数的封装
Fa.pca<-function(X,cor,z0){
  source("Pca.R")
  #相关性检验
  cor_test=corr.test(X,use = "complete")
  #主成分分析
  pca<-Pca(X,cor,z0)
  #因子个数
  k=pca$"Number of Chosen Components"
  #总方差的解释
  va=pca$"Total Variance Explained"
  #方差
  lambda=pca$"Variance of Component"
  #载荷矩阵
  A=pca$"Component Matrix"
  #共同度&特殊度
  C=data.frame(pca$Communalities,1-pca$Communalities)
  colnames(C)=c("Communalities","Specificities")
  #成分得分系数矩阵
  Coeff=c()
  for (i in 1:k){
    a=pca$`Coefficient of Components`[,i]/sqrt(pca$`Variance of Component`[i])
    Coeff=cbind(Coeff,a)
  }
  colnames(Coeff)=paste("factor",c(1:k),sep="")
  
  #正交旋转
  X.var=principal((X),nfactors = k,rotate="varimax",scores=T)
  
  list("Correlation"=cor_test,"Number of Chosen Components"=k,"Total Variance Explained"=va,
       "Variance of Component"=lambda,"Component Matrix"=A,"Communities & Specialities"=C,
       "Component Score Coefficient Matrix"=Coeff,"Varimax"=X.var,"weights of Varimax"=X.var$weights,
       "scores"=X.var$scores
  )
  
}


## ----------
## (2)调用封装函数Fa.pca()分析例6.2
Y=read.csv("eg6_2.csv",header=T)
rownames(Y)=Y[,1]
Y=Y[,-1]
Y[,6]=1/Y[,6]
Y[,7]=1/Y[,7]
head(Y)

Y.fa.pca=Fa.pca(Y,cor=T,0.8)
Y.fa.pca


## ----------
## (3)可视化
## 先分组
type=c()
for(i in 1:nrow(Y)){
  if(Y.fa.pca$scores[i,1]>0 && Y.fa.pca$scores[i,2]>0){
    type=c(type,1)
  }
  if(Y.fa.pca$scores[i,1]>0 && Y.fa.pca$scores[i,2]<0){
    type=c(type,2)
  }
  if(Y.fa.pca$scores[i,1]<0 && Y.fa.pca$scores[i,2]<0){
    type=c(type,3)
  }
  if(Y.fa.pca$scores[i,1]<0 && Y.fa.pca$scores[i,2]>0){
    type=c(type,4)
  }
}

## 再作图
Y=cbind(Y,type)
colnames(Y)=c("X1","X2","X3","X4","X5","X6","X7","X8","type")
Y.var<-prcomp(Y[,1:8], scale. = TRUE)
library(ggbiplot)
ggbiplot(Y.var, obs.scale = 1, var.scale = 1,
         groups = as.factor(type), ellipse = TRUE, circle = TRUE,size=2) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')



## ------------------------------------
## 例6-3
## ------------------------------------
## (1)加权最小二乘法
## Step1：求解初始因子
Z=read.csv("eg6_3.csv",header=T)
rownames(Z)=Z[,1]
Z=Z[,-1]
p=ncol(Z)
k=Pca(Z,cor=T,1)$"Number of Chosen Components"
k
## 加权最小二乘法
library(GPArotation)
Z.wls<-fa(cor(Z),nfactors = k,n.obs = nrow(Z),rotate = "none",scores = T,fm="wls")
Z.wls

## Step2：正交旋转
Z.wls.var<-fa(cor(Z),nfactors = k,n.obs = nrow(Z),rotate = "varimax",scores = T,fm="wls")
Z.wls.var
Z.wls.var$weights

## Step 3：斜交旋转
Z.wls.pro<-fa(cor(Z),nfactors = k,n.obs = nrow(Z),rotate = "promax",scores = T,fm="wls")
Z.wls.pro
Struc=Z.wls.pro$loadings[1:p,1:k]%*%Z.wls.pro$score.cor #计算Structure
colnames(Struc)=c("WLS1","WLS2","WLS3")
Struc

## Step4：因子分析图
fa.diagram(Z.wls.var,simple=FALSE)


## ----------
## (2)最大似然法
## Step1：求解初始因子
Z.ml<-fa(cor(Z),nfactors = 3,n.obs = nrow(Z),rotate = "none",scores = T,fm="ml")
Z.ml

## Step2: 正交旋转
Z.ml.var<-fa(cor(Z),nfactors = 3,n.obs = nrow(Z),rotate = "varimax",scores = T,fm="ml")
Z.ml.var
Z.ml.var$weights

## Step3：斜交旋转
Z.ml.pro<-fa(cor(Z),nfactors = 3,n.obs = nrow(Z),rotate = "promax",scores = T,fm="ml")
Z.ml.pro
Struc2=Z.ml.pro$loadings[1:p,1:3]%*%Z.ml.pro$score.cor #计算Structure
colnames(Struc2)=c("ML1","ML2","ML3")
Struc2

## Step4：因子分析图
Z.ml.2 <- fa(cor(Z),nfactors = 2,n.obs = nrow(Z),rotate = "varimax",scores = T,fm="ml")
factor.plot(Z.ml.2,labels = rownames(Z.ml.2$loadings))


## ----------
## (3)主成分法
Z.fa.pca<-Fa.pca(Z,cor=T,1)
Z.fa.pca[c(8,9)]


## ----------
## (4)因子得分
F=(54.381*Z.fa.pca$scores[,1]+22.077*Z.fa.pca$scores[,2]+10.647*Z.fa.pca$scores[,3])/87.105
Z.Score=data.frame(Z.fa.pca$scores,F)
Z.Score

## 因子得分图
ggplot(Z,aes(x=Z.Score[,1], y = Z.Score[,2]))+
  geom_point()+
  labs(x = "RC1", y = "RC2")
