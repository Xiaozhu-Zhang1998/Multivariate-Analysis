## Chapter 8

## ------------------------------------
## 例8-1
## ------------------------------------
rm(list=ls())
## 输入数据
X<-data.frame(
  X1=c(191, 193, 189, 211, 176, 169, 154, 193, 176, 156,
       189, 162, 182, 167, 154, 166, 247, 202, 157, 138),
  X2=c(36, 38, 35, 38, 31, 34, 34, 36, 37, 33,
       37, 35, 36, 34, 33, 33, 46, 37, 32, 33),
  X3=c(50, 58, 46, 56, 74, 50, 64, 46, 54, 54,
       52, 62, 56, 60, 56, 52, 50, 62, 52, 68),
  Y1=c( 5, 12, 13, 8, 15, 17, 14, 6, 4, 15,
        2, 12, 4, 6, 17, 13, 1, 12, 11, 2),
  Y2=c(162, 101, 155, 101, 200, 120, 215, 70, 60, 225,
       110, 105, 101, 125, 251, 210, 50, 210, 230, 110),
  Y3=c(60, 101, 58, 38, 40, 38, 105, 31, 25, 73,
       60, 37, 42, 40, 250, 115, 50, 120, 80, 43)
)

## ----------
## (1)变量间的相关关系
X=scale(X,center=T,scale = T)
## 第一组变量间的相关系数
cor(X[,c(1,2,3)])
## 第二组变量间的相关系数
cor(X[,c(4,5,6)])
## 第一、二组变量间的相关系数
cor(X[,c(1,2,3)],X[,c(4,5,6)])


## ----------
## (2)典型相关系数的检验
can.x<-cancor(X[,1:3],X[,4:6])
can.x$cor

## 编写自定义函数进行检验
cor.test<-function(X,Y,n,p,q,alpha){
  can<-cancor(X,Y)
  j=min(p,q)
  count=0
  QSTA=c()
  PVAL=c()
  for(i in 1:j){
    Lambda=prod(1-can$cor[i:j]^2)
    Q=-(n-i-1/2*(p+q+1))*log(Lambda)
    pvalue=pchisq(Q,(p-i+1)*(q-i+1),lower.tail = F)
    if(pvalue<alpha){
      count=count+1
    }
    QSTA=c(QSTA,Q)
    PVAL=c(PVAL,pvalue)
  }
  return(list("count"=count,"Q-statistics"=QSTA,"pvalue"=PVAL))
}

## 调用
cor.test(X[,1:3],X[,4:6],20,3,3,0.05)
cor.test(X[,1:3],X[,4:6],20,3,3,0.10)


## ----------
## (3)典型权重、典型载荷和交叉载荷
## (3.1)典型权重
## 第1组变量
can.x$xcoef
## 第2组变量
can.x$ycoef

U<-X[,1:3]%*%can.x$xcoef
V<-X[,4:6]%*%can.x$ycoef
colnames(U)<-c("U1","U2","U3")
colnames(V)<-c("V1","V2","V3")
Score<-data.frame(U,V)
Score

## (3.2)典型载荷
## 第1组变量及其典型变量的相关系数
cor(X[,1:3],U)
## 第2组变量及其典型变量的相关系数
cor(X[,4:6],V)

## (3.3)典型交叉载荷
## 第1组变量与第2组典型变量的相关系数
cor(X[,1:3],V)
## 第2组变量与第1组典型变量的相关系数
cor(X[,4:6],U)


## ----------
## (4)冗余分析
## 共同方差的比例
L1=cor(X[,4:6],V)
L1^2  #典型载荷的平方
apply(L1^2,1,sum)  #验证，按行求和为1
R2=colMeans(L1^2)  #计算典型变量所解释的共同方差的比例
R2
cumsum(R2)

## 解释的方差比例
P=can.x$cor^2

## 冗余指数
R2*P
cumsum(R2*P)
data.frame("R2"=R2,"Cumulative R2"=cumsum(R2),
           "Cor2"=P,"Redundancy"=R2*P,"Cumulative Redundancy"=cumsum(R2*P))
R2=colMeans(cor(X[,1:3],U)^2)
R2*P
data.frame("R2"=R2,"Cumulative R2"=cumsum(R2),
           "Cor2"=P,"Redundancy"=R2*P,"Cumulative Redundancy"=cumsum(R2*P))



## ------------------------------------
## 例8-2
## ------------------------------------
## 封装函数
cca<-function(X,p,q,alpha){
  X=scale(X,center = T,scale = T)
  #变量间的相关性
  X1=X[,1:p]
  X2=X[,(p+1):(p+q)]
  cor1=cor(X1)
  cor2=cor(X2)
  cor12=cor(X1,X2)
  
  #典型相关系数及其检验
  can.x<-cancor(X1,X2)   #进行典型相关分析
  test<-cor.test(X1,X2,nrow(X),p,q,alpha)  #进行检验
  
  U<-as.matrix(X1)%*%as.matrix(can.x$xcoef)
  V<-as.matrix(X2)%*%as.matrix(can.x$ycoef)
  colnames(U)<-paste("U",c(1:p),sep = "",collapse = NULL)
  colnames(V)<-paste("V",c(1:q),sep = "",collapse = NULL)
  U=U[,1:min(p,q)]
  V=V[,1:min(p,q)]
  #典型载荷
  L11=cor(X1,U)
  L22=cor(X2,V)
  #交叉载荷
  L12=cor(X1,V)
  L21=cor(X2,U)
  
  #冗余分析
  #第一组变量
  R2_1=colMeans(L11^2)
  P=can.x$cor^2
  redun1<-data.frame("Proportion"=R2_1,"Cumulative Proportion"=cumsum(R2_1),
                     "R-Squared"=P,"Canonical Proportion"=R2_1*P,"Cumulative Proportion"=cumsum(R2_1*P))
  #第二组变量
  R2_2=colMeans(L22^2)
  redun2<-data.frame("R2"=R2_2,"Cumulative R2"=cumsum(R2_2),
                     "Cor2"=P,"Redundancy"=R2_2*P,"Redundancy"=cumsum(R2_2*P))
  
  return(list("cor1"=cor1,"cor2"=cor2,"cor12"=cor12,"CanonicalCor"=can.x$cor,"cor.test"=test,
              "xcoef"=can.x$xcoef,"ycoef"=can.x$ycoef,
              "Loadings_xU"=L11,"Loadings_yV"=L22,"LoadingsxV"=L12,"LoadingsyU"=L21,
              "RedunAna_x"=redun1,"RedunAna_y"=redun2,"scores_U"=U,"scores_V"=V))
}


## 调用cca()分析例8.2
## 读入数据
Y=read.csv("eg8-2.csv")
rownames(Y)=Y[,1]
Y=Y[,-1]
head(Y)

## 调用函数
cca.Y=cca(Y,4,6,0.05)
cca.Y[c(1:13)]

## 可视化:散点图
library(ggplot2)
Yplot=data.frame(cca.Y$scores_U[,1],cca.Y$scores_U[,2])
ggplot(Yplot,aes(x=Yplot[,1], y = Yplot[,2]))+
  geom_point()+
  labs(x = "U1", y = "V1")



## ------------------------------------
## 例8-3
## ------------------------------------
## 读入数据
Z=read.csv("eg8-3.csv")
rownames(Z)=Z[,1]
Z=Z[,-1]
head(Z)

## 调用函数Pca()
cca.Z=cca(Z,6,4,0.05)
cca.Z[c(1:13)]

## 可视化：散点图
Zplot=data.frame(cca.Z$scores_U[,1],cca.Z$scores_U[,2])
ggplot(Zplot,aes(x=Zplot[,1], y = Zplot[,2]))+
  geom_point()+
  labs(x = "U1", y = "V1") 

