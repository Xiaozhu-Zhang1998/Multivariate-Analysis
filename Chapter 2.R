## 加载数据
X=read.table("clipboard",header=T)
## ----------------------------------------
library(mvnormtest)

## 检验每一个变量是否服从正态分布
par(mfrow=c(4,2))
for(i in 4:11){
  k1=mshapiro.test(t(X[,i]))
  qqnorm(X[,i])
    if(k1$p.value>0.05){
    print(list(i-3,k1$p.value))
  }
}

## 对于服从正态分布的单个变量，检验其向量是否服从正态分布
mshapiro.test(t(X[,c(4,5,6,10)]))
## 4个分量构成的向量不服从正态分布

## 检验3个分量构成的向量是否服从正态分布
for(i in c(4,5,6,10)){
  for(j in c(5,6,10)){
    for(k in c(6,10)){
      if(i<j && j<k){
        te=mshapiro.test(t(X[,c(i,j,k)]))
        if(te$p.value>0.05){
          print(list(c(i-3,j-3,k-3),te$p.value))
        }
      }
    }
  }
}

## 得到了两组向量，服从正态分布

## ----------------------------------------
## 进行方差齐性检验(自编程序)

## attach(X)
r=3
G1 <- subset(X,"行业"=="电力、煤气及水的生产和供应业")
G2 <- subset(X,"行业"=="房地行业")
G3 <- subset(X,"行业"=="信息技术业")
n1=nrow(G1)
n2=nrow(G2)
n3=nrow(G3)
n=n1+n2+n3

cov.test <- function(p,nr){ #输入c(检验变量个数,检验变量编号(向量))
  nr=nr+3
  L1=(n1-1)*cov(G1[,nr])
  L2=(n2-1)*cov(G2[,nr])
  L3=(n3-1)*cov(G3[,nr])
  L=L1+L2+L3
  
  ## 构造M统计量
  A=(n-r)*log(exp(1),det(L/(n-r)))
  B1=(n1-1)*log(exp(1),det(L1/(n1-1)))
  B2=(n2-1)*log(exp(1),det(L2/(n2-1)))
  B3=(n3-1)*log(exp(1),det(L3/(n3-1)))
  M=A-B1-B2-B3
  ## 寻找分布
  if(n1==n2 && n2==n3){
    d1=(2*p^2+3*p-1)*(r-1)/(6*(p+1)*r*(n-1)) 
    d2=(p-1)*(p+2)*(r^2+r+1)/(6*r^2*(n-1)^2)
  }
  
  if(n1!=n2 || n1!=n3 || n2!=n3){
    d1=(2*p^2+3*p-1)/(6*(p+1)*(r-1))*(1/(n1-1)+1/(n2-1)+1/(n3-1)-1/(n-r)) 
    d2=(p-1)*(p+2)/(6*(r-1))*(1/(n1-1)^2+1/(n2-1)^2+1/(n3-1)^2-1/(n-r)^2)
  }
  
  f1=p*(p+1)*(r-1)/2
  f2=(f1+2)/(d2-d1^2)
  b=f1/(1-d1-f1/f2)
  FS=M/b
  pvalue=pf(FS,f1,f2,ncp=0,lower.tail=FALSE,log.p=FALSE)
  list("M"=M,"b"=b,"p-value"=pvalue)
}

cov.test(3,c(1,2,3))
cov.test(3,c(1,3,7))

## ----------------------------------------
## 方差分析
## 多元方差分析
attach(X)
aggregate(cbind(净资产收益率,总资产报酬率,资产负债率,销售增长率),by=list(行业),FUN=mean)

y.1 <- as.data.frame(cbind(行业,净资产收益率,总资产报酬率,资产负债率))
type1 <- as.factor(y.1$行业)
fit1 <- manova(cbind(净资产收益率,总资产报酬率,资产负债率)~type1,data=y.1)
summary(fit1,test=c("Wilks"))

y.2 <- as.data.frame(cbind(行业,净资产收益率,资产负债率,销售增长率))
type2 <- as.factor(y.2$行业)
fit2 <- manova(cbind(净资产收益率,资产负债率,销售增长率)~type2,data=y.2)
summary(fit2,test=c("Wilks"))

## 单因素方差分析
y.3 <- as.data.frame(cbind(行业,净资产收益率,总资产报酬率,资产负债率,销售增长率))
type3 <- as.factor(y.3$行业)
fit3 <- manova(cbind(净资产收益率,总资产报酬率,资产负债率,销售增长率)~type3)
summary.aov(fit3)

## 事后比较
library(agricolae)
back.test <- function(i){
  attach(X)
  fit <- aov(i~行业,data=X)
  #bonferroni检验
  out.LSD <- LSD.test(fit,"行业",p.adj="bonferroni")
  #SNK检验
  out.SNK <- SNK.test(fit,"行业")
  #TukeyHSD检验
  out.TUK=TukeyHSD(fit)
  #Scheffe检验
  out.SHF <- scheffe.test(fit,"行业")
  
  par(mfrow=c(2,2))
  plot(out.LSD)
  plot(out.SNK)
  plot(out.TUK)
  plot(out.SHF)
  list(LSD=out.LSD$group,SNK=out.SNK$group,TUK=out.TUK$行业,SHF=out.SHF$group)
}

back.test(净资产收益率)
back.test(总资产报酬率)
