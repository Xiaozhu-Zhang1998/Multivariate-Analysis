## Chapter 7

## ------------------------------------
## 例7-1
## ------------------------------------
rm(list=ls())
## 读入数据
X=read.csv("eg7-1.csv")
X=X[,-1]
X=subset(X,degree<5)
X$degree<-factor(X$degree,labels=c("Less than high school","High school","Junior college","Bachelor","Graduate"))
X$race<-factor(X$race,labels=c("white","black","others"))
head(X)


## ----------
## (1)列联表及其检验
addmargins(table(X),c(1,2))
addmargins(prop.table(table(X)))
chisq.test(table(X))


## ----------
## (2)利用函数ca()进行因子分析
library(ca)
ca.x=ca(table(X))
ca.x


## ----------
## (3)获取更多信息
## 二维图
ca.x$rowcoord
ca.x$colcoord
plot(ca.x)

## 因子载荷矩阵
## (第一个变量)
dim11=c(-1.217441,-0.206512,-0.801486,1.902870,1.129272)
dim12=c(-1.688145,0.785718,0.786794,-0.828163,-0.167158)
a11=dim11*ca.x$sv[1]  #第1个因子与第1个变量
a12=dim12*ca.x$sv[2]  #第2个因子与第1个变量
loadings1<-data.frame(a11,a12)
rownames(loadings1)=c("Less than high school","High school","Junior college","Bachelor","Graduate")
colnames(loadings1)=c("dim1","dim2")
loadings1

## (第二个变量)
dim21=c(0.297415,-2.767419,1.189338)
dim22=c(0.323306,-0.547246,-4.187140)
a21=dim21*ca.x$sv[1]  #第1个因子与第2个变量
a22=dim22*ca.x$sv[2]  #第2个因子与第2个变量
loadings2<-data.frame(a21,a22)
rownames(loadings2)=c("white","black","others")
colnames(loadings2)=c("dim1","dim2")
loadings2

## CTR(i)
## (第一个变量)
CTR11=ca.x$rowmass*dim11^2  #第1个因子与第1个变量
CTR12=ca.x$rowmass*dim12^2  #第2个因子与第1个变量
CTR1=data.frame(CTR11,CTR12)
rownames(CTR1)=c("Less than high school","High school","Junior college","Bachelor","Graduate")
colnames(CTR1)=c("dim1","dim2")
CTR1

## (第二个变量)
CTR21=ca.x$colmass*dim21^2  #第1个因子与第1个变量
CTR22=ca.x$colmass*dim22^2  #第2个因子与第1个变量
CTR2=data.frame(CTR21,CTR22)
rownames(CTR2)=c("white","black","others")
colnames(CTR2)=c("dim1","dim2")
CTR2


## ----------
## (4)行剖面和列剖面
## 行剖面
RowPro=c()
for(i in 1:5){
  z1=table(X)[i,]/rowSums(table(X))[i]
  RowPro=rbind(RowPro,z1)
}
RowPro=rbind(RowPro,ca.x$colmass)
rownames(RowPro)=c("Less than high school","High school","Junior college","Bachelor","Graduate","MASS")
RowPro

Tab2=table(X)
colnames(Tab2)<-NULL
plot(ca(Tab2))

## 列剖面
ColPro=c()
for(i in 1:3){
  z2=table(X)[,i]/colSums(table(X))[i]
  ColPro=cbind(ColPro,z2)
}
ColPro=cbind(ColPro,ca.x$rowmass)
colnames(ColPro)=c("white","black","others","MASS")
ColPro

Tab1=table(X)
rownames(Tab1)<-NULL
plot(ca(Tab1))



## ------------------------------------
## 例7-2
## ------------------------------------
## (1)利用ca()函数进行分析
## 读入数据
Y=read.csv("eg7-2.csv")
rownames(Y)=Y[,1]
Y=Y[,-1]
head(Y)

## 调用ca()
ca.y=ca(Y,nd=2)
ca.y


## ----------
## (2)计算因子载荷
dim1=c(-0.826507,1.196315,-1.281792,-0.648884)
dim2=c(0.704335,0.018102,-1.725246,-2.557861)
a1=dim1*ca.y$sv[1]  #第1个因子与第2个变量
a2=dim2*ca.y$sv[2]  #第2个因子与第2个变量
loadings<-data.frame(a1,a2)
rownames(loadings)=c("工资性收入","家庭经营纯收入","财产性收入","转移性收入")
colnames(loadings)=c("dim1","dim2")
loadings


## ----------
## (3)CTR(i)
CTR1=ca.y$colmass*dim1^2  #第1个因子与第2个变量
CTR2=ca.y$colmass*dim2^2  #第2个因子与第2个变量
CTR=data.frame(CTR1,CTR2)
rownames(CTR)=c("工资性收入","家庭经营纯收入","财产性收入","转移性收入")
colnames(CTR)=c("dim1","dim2")
CTR


## ----------
## (4)二维图
ca.y$rowcoord
ca.y$colcoord
plot(ca.y)
