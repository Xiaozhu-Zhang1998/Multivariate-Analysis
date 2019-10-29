## Chapter 5

## ------------------------------------
## 例5-1
## ------------------------------------
rm(list=ls())
## 读入数据

## ----------
## (1)选择主成分个数
X=read.csv("eg5_1.csv",header=T)
rownames(X)=X[,1]
X=X[,-1]
X1=X[,c(3,5:8)]
p=ncol(X1)
head(X1)

## Step1：主成分方差(特征根)
PCA.X=princomp(X1,cor=TRUE)
summary(PCA.X)

lambda=PCA.X$sdev^2
lambda

## Step2：碎石图
library(ggplot2)
ggplot(as.data.frame(lambda),aes(x=c("Comp1","Comp2","Comp3","Comp4","Comp5"), y=as.numeric(lambda), group=1)) +
  geom_line(colour="skyblue3", size=0.7)+geom_point(size=2,colour="skyblue3")+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black"))
## 选择前3个主成分


## ----------
## (2)系数矩阵，因子负荷量，和主成分对原始变量的方差贡献率
## Step1: 系数矩阵
k=3
eigV=PCA.X$loadings[1:p,1:k]
eigV

## Step2: 因子负荷量
ComM=c()
for (i in 1:k){
  ComM=cbind(ComM,eigV[,i]*sqrt(lambda[i]))
}
ComM

## Step3: 主成分对原始变量的方差贡献率
Communal=c()
for (i in 1:p){
  ext=sum(ComM[i,]^2)
  Communal=rbind(Communal,ext)
}
rownames(Communal)=colnames(X1)
Communal


## ----------
## (3)计算得分
## 计算各个样品的主成分值
PCA.X.Score=as.matrix(scale(X1))%*%as.matrix(eigV)
head(PCA.X.Score)

## 计算得分
scores=c()
for (i in 1:nrow(X1)){
  scores1=0
  for (j in 1:k){
    w=lambda[j]/sum(lambda[1:k])
    scores1=scores1+PCA.X.Score[i,j]*w
  }
  scores=c(scores,scores1)
}
X.score=cbind(X1,scores)
head(X.score)



## ------------------------------------
## 例5-2
## ------------------------------------
## 封装函数
Pca <- function(X,cor,z0){
  p=ncol(X)
  library(psych)
  PCA.X=princomp(X,cor) #cor=TRUE or cor=FALSE
  
  #主成分方差(特征根)
  lambda=PCA.X$sdev^2
  #碎石图
  library(ggplot2)
  pl<- ggplot(as.data.frame(lambda),aes(x=as.character(c(1:p)), y=as.numeric(lambda), group=1)) +
    geom_line(colour="skyblue3", size=0.7)+geom_point(size=2,colour="skyblue3")+
    geom_hline(yintercept = z0, linetype = 2)+
    theme(panel.background = element_blank(),axis.line = element_line(colour = "black"))
  
  #选择主成分个数
  k=0
  for(i in 1:p){
    if(lambda[i]>z0){
      k=k+1
    }
  }
  
  #系数矩阵(特征向量)
  eigV=PCA.X$loadings[1:p,1:k]
  
  #Component Matrix(因子负荷量)
  ComM=c()
  for (i in 1:k){
    ComM=cbind(ComM,eigV[,i]*sqrt(lambda[i]))
  }
  
  #Communalities(提取信息量)
  Communal=c()
  for (i in 1:p){
    ext=sum(ComM[i,]^2)
    Communal=rbind(Communal,ext)
  }
  rownames(Communal)=colnames(X)
  
  #计算得分
  PCA.X.Score=as.matrix(scale(X))%*%as.matrix(eigV)
  scores=c()
  for (i in 1:nrow(X)){
    scores1=0
    for (j in 1:k){
      w=lambda[j]/sum(lambda[1:k])
      scores1=scores1+PCA.X.Score[i,j]*w
    }
    scores=c(scores,scores1)
  }
  
  X.score=cbind(PCA.X.Score,scores)
  
  
  
  return(list("Total Variance Explained"=summary(PCA.X),"Variance of Component"=lambda,"Number of Chosen Components"=k,
              "Component Matrix"=ComM,"Communalities"=Communal,"Coefficient of Components"=eigV,
              "scores"=X.score,"Scree Plot"=pl))
}


## 调用Pca()分析例5.2
## 读入数据
Y=read.csv("eg5-2.csv",header=T)
rownames(Y)=Y[,1]
Y=Y[,-1]

## 调用函数
Y.Pca<-Pca(Y,cor=T,1)
Y.Pca

## 可视化:散点图
## 先分类
type=c()
for(i in 1:nrow(Y)){
  if(Y.Pca$scores[i,1]>0 && Y.Pca$scores[i,2]>0){
    type=c(type,1)
  }
  if(Y.Pca$scores[i,1]>0 && Y.Pca$scores[i,2]<0){
    type=c(type,2)
  }
  if(Y.Pca$scores[i,1]<0 && Y.Pca$scores[i,2]<0){
    type=c(type,3)
  }
  if(Y.Pca$scores[i,1]<0 && Y.Pca$scores[i,2]>0){
    type=c(type,4)
  }
}
type

## 作图
Y=cbind(Y,type)
colnames(Y)=c("X1","X2","X3","X4","X5","X6","X7","X8","X9","type")
y.pca <- prcomp(Y[,1:9], scale. = TRUE)
library(ggbiplot)
ggbiplot(y.pca, obs.scale = 1, var.scale = 1,
         groups = as.factor(type), ellipse = TRUE, circle = TRUE,size=2) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')



## ------------------------------------
## 例5-3
## ------------------------------------
## 读入数据
Z=read.csv("eg5-3.csv",header=T)
attach(Z)
rownames(Z)=Z[,1]
Z=Z[,-1]
Z[,6]=1/Z[,6]
Z[,7]=1/Z[,7]
head(Z)

## 调用函数Pca()
Z.Pca<-Pca(Z,cor=T,0.8)
Z.Pca

## 计算得分排名
y1=as.numeric((-1)*Z.Pca$scores[,1])
rk=nrow(Z)-rank(y1)+1 
demo.rk=data.frame(y1,rk)
colnames(demo.rk)=c("y1","名次")
demo.rk

