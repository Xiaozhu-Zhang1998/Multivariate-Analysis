#读取数据
setwd("C:\\Users\\小竹子\\Documents\\Multivariate Statistical Analysis\\期末论文")
data=read.csv("forestfires.csv")
head(data,5)
attach(data)
#---坐标变量,说明观测点分布较为均匀,数据具有说服力---#
#--------------------------------------------------
#--------------------------------------------------
#--------------------------------------------------
library(ggplot2)
position=data.frame(X,Y)
ggplot(position,aes(x = X, y = Y))+
  geom_point(size=2)+
  labs(x = "x-axis spatial coordinate", y = "y-axis spatial coordinate")+
  scale_x_continuous(breaks=seq(1,9,1))+
  scale_y_continuous(breaks=seq(2,9,1))

#---1. 对area进行分类---#
#--------------------------------------------------
#--------------------------------------------------
#--------------------------------------------------
#想法,boxplot
ggplot(data,aes(x=area))+geom_histogram()
#code
A=sort(area)
k=1.5 #设定临界值
cluster=c()
type=1
start=1
for(i in 1:length(A)){
  q=quantile(A[start:i])
  if(q[5]-q[4]<=k*(q[4]-q[2])){
    cluster=c(cluster,type)
  }else{
    type=type+1
    cluster=c(cluster,type)
    start=i+1
  }
}

degree=c()
for(i in 1:length(area)){
  for(j in 1:length(area)){
    if(area[i]==A[j]){
      degree[i]=cluster[j]
    }
  }
}

data<-cbind(data,degree)
#boxplot after cluster
ggplot(data,aes(x=as.factor(degree),y=area))+
  geom_boxplot()+
  labs(x = "degree", y = "area")

#---2. 典型相关分析---#
#--------------------------------------------------
#--------------------------------------------------
#--------------------------------------------------
data.X=data[,5:8]
data.Y=data[,9:12]
can=cancor(data.X,data.Y)
can$cor
#序贯检验
cor.test<-function(X,Y,alpha){
  n=nrow(X)
  p=ncol(X)
  q=ncol(Y)
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
cor.test(data.X,data.Y,alpha=0.05)
#计算U,V
can$xcoef
can$ycoef
U=as.matrix(data.X)%*%as.matrix(can$xcoef[,1:3])
colnames(U)=paste("U",1:3,sep="")
V=as.matrix(data.Y)%*%as.matrix(can$ycoef[,1:3])
colnames(V)=paste("V",1:3,sep="")
#典型载荷
cor(data.X,U)
cor(data.Y,V)

#---3. 对应分析---#
#--------------------------------------------------
#--------------------------------------------------
#--------------------------------------------------
t=table(month,degree)
t
#行剖面
addmargins(prop.table(t,1),2)
tmonth=c() #行剖面作图
tD=c()
freq=c()
pt=prop.table(t,1)
for(i in 1:nrow(t)){
  for(j in 1:ncol(t)){
    tmonth=c(tmonth,rownames(t)[i])
    tD=c(tD,colnames(t)[j])
    freq=c(freq,t[i,j])
  }
}
ggplot(data.frame(month=tmonth,degree=tD,percentage=freq),
       aes(month,percentage,fill=degree))+
  geom_bar(stat="identity",position="fill")

#列剖面
addmargins(prop.table(t,2),1)
ggplot(data.frame(month=tmonth,degree=tD,percentage=freq), #列剖面作图
       aes(degree,percentage,fill=month))+
  geom_bar(stat="identity",position="fill")

#对应分析
library(ca)
ca.month=ca(table(month,degree))
plot(ca.month)

#---3.1 对数线性模型---#
#--------------------------------------------------
#--------------------------------------------------
#--------------------------------------------------
table(degree,day)
serious<-c()
weekend<-c()
for(i in 1:nrow(data)){
  if(day[i]=="fri" || day[i]=="sat" || day[i]=="sun"){
    weekend=c(weekend,1)
  }else{
    weekend=c(weekend,0)
  }
  if(degree[i]==4 || degree[i]==5 || degree[i]==6){
    serious=c(serious,1)
  }else{
    serious=c(serious,0)
  }
}
day.test<-data.frame(weekend,serious)
table(day.test)
summary(table(day.test))
library(MASS)
logmol=loglm(~weekend*serious,table(day.test))
logmol$param


#---4. 判别分析(距离判别法)---#
#--------------------------------------------------
#--------------------------------------------------
#--------------------------------------------------
r=4
#子总体
#将4、5、6进行合并(因为样本数量过少)
G1 <- subset(data,degree==1)
G2 <- subset(data,degree==2)
G3 <- subset(data,degree==3)
G4 <- rbind(subset(data,degree==4),
            subset(data,degree==5),
            subset(data,degree==6))
n1=nrow(G1)
n2=nrow(G2)
n3=nrow(G3)
n4=nrow(G4)
n=n1+n2+n3+n4
#判断协方差阵是否同质
cov.test <- function(p,nr){ #输入c(检验变量个数,检验变量编号(向量))
  #构造M统计量
  L1=(n1-1)*cov(G1[,nr])
  L2=(n2-1)*cov(G2[,nr])
  L3=(n3-1)*cov(G3[,nr])
  L4=(n4-1)*cov(G4[,nr])
  L=L1+L2+L3+L4
  A=(n-r)*log(det(L/(n-r)))
  B1=(n1-1)*log(det(L1/(n1-1)))
  B2=(n2-1)*log(det(L2/(n2-1)))
  B3=(n3-1)*log(det(L3/(n3-1)))
  B4=(n4-1)*log(det(L4/(n4-1)))
  M=A-B1-B2-B3-B4
  #寻找分布
  if(n1==n2 && n2==n3 && n3==n4){
    d1=(2*p^2+3*p-1)*(r-1)/(6*(p+1)*r*(n-1)) 
    d2=(p-1)*(p+2)*(r^2+r+1)/(6*r^2*(n-1)^2)
  }
  if(n1!=n2 || n1!=n3 || n2!=n3 || n3!=n4){
    d1=(2*p^2+3*p-1)/(6*(p+1)*(r-1))*(1/(n1-1)+1/(n2-1)+1/(n3-1)-1/(n-r)) 
    d2=(p-1)*(p+2)/(6*(r-1))*(1/(n1-1)^2+1/(n2-1)^2+1/(n3-1)^2-1/(n-r)^2)
  }
  f1=p*(p+1)*(r-1)/2
  f2=(f1+2)/(d2-d1^2)
  b=f1/(1-d1-f1/f2)
  FS=M/b
  pvalue=pf(FS,f1,f2,lower.tail=FALSE,log.p=FALSE)
  list("M"=M,"b"=b,"p-value"=pvalue)
}
#test
cov.test(4,5:8)
#因此采用方差不等计算
#计算协方差阵和均值向量
S1=cov(G1[,5:8])
S2=cov(G2[,5:8])
S3=cov(G3[,5:8])
S4=cov(G4[,5:8])

m1=colMeans(G1[,5:8])
m2=colMeans(G2[,5:8])
m3=colMeans(G3[,5:8])
m4=colMeans(G4[,5:8])
#预测及回判
d1=c()
for (i in 1:nrow(data)){
  dis1=as.matrix(data[i,5:8]-m1)%*%as.matrix(solve(S1))%*%as.matrix(t(data[i,5:8]-m1))
  d1=c(d1,dis1)
}
d2=c()
for (i in 1:nrow(data)){
  dis2=as.matrix(data[i,5:8]-m2)%*%as.matrix(solve(S2))%*%as.matrix(t(data[i,5:8]-m2))
  d2=c(d2,dis2)
}
d3=c()
for (i in 1:nrow(data)){
  dis3=as.matrix(data[i,5:8]-m3)%*%as.matrix(solve(S3))%*%as.matrix(t(data[i,5:8]-m3))
  d3=c(d3,dis3)
}
d4=c()
for (i in 1:nrow(data)){
  dis4=as.matrix(data[i,5:8]-m4)%*%as.matrix(solve(S4))%*%as.matrix(t(data[i,5:8]-m4))
  d4=c(d4,dis4)
}
#结果整理
D1=c()
dist=c()
for(i in 1:nrow(data)){
  if(d1[i]==min(d1[i],d2[i],d3[i],d4[i])){
    D1=c(D1,1)
  }
  if(d2[i]==min(d1[i],d2[i],d3[i],d4[i])){
    D1=c(D1,2)
  }
  if(d3[i]==min(d1[i],d2[i],d3[i],d4[i])){
    D1=c(D1,3)
  }
  if(d4[i]==min(d1[i],d2[i],d3[i],d4[i])){
    D1=c(D1,4)
  }
  dist=c(dist,min(d1[i],d2[i],d3[i],d4[i]))
}
#计算精度
ndegree=c()
for(i in 1:nrow(data)){
  if(degree[i]==1 || degree[i]==2 || degree[i]==3){
    ndegree[i]=degree[i]
  }else{
    ndegree[i]=4
  }
}
table(ndegree,D1)
sum(diag(prop.table(table(ndegree,D1))))

#---5. 寻找精度较低的原因---#
#--------------------------------------------------
#--------------------------------------------------
#--------------------------------------------------
#1. 正态分布
library(mvnormtest)
Glist=list(G1,G2,G3,G4)
W=c()
P.value=c()
for(i in 1:4){
  test=mshapiro.test(t(Glist[[i]][,5:8]))
  W=c(W,test$statistic)
  P.value=c(P.value,test$p.value)
}
MNtest<-data.frame(W,P.value)
rownames(MNtest)=paste("G",1:4,sep="")
MNtest

#如何改进？
ld=lda(ndegree~FFMC+DMC+DC+ISI,data)
new=predict(ld)
t2=table(ndegree,new$class)
t2
sum(diag(prop.table(t2)))


qd=qda(ndegree~FFMC+DMC+DC+ISI,data)
new=predict(qd)
t3=table(ndegree,new$class)
t3
sum(diag(prop.table(t3)))

#---6. Logistic---#
#--------------------------------------------------
#--------------------------------------------------
#--------------------------------------------------
#-------------------这个比较靠谱--------------------------------------------
data=cbind(data,ndegree,D1)
dd=subset(data,D1==1)
Glist=list(subset(data,ndegree==1),
           subset(data,ndegree==2),
           subset(data,ndegree==3),
           subset(data,ndegree==4))
A=Glist[[1]]
B=rbind(Glist[[2]],Glist[[3]],Glist[[4]])
A$ndegree=1
B$ndegree=0
mol.1=glm(ndegree~FFMC+DMC+DC+ISI,
        family = binomial(link = logit),data=rbind(A,B))
#---
A=Glist[[2]]
B=rbind(Glist[[3]],Glist[[4]])
A$ndegree=1
B$ndegree=0
mol.2=glm(ndegree~FFMC+DMC+DC+ISI,
        family = binomial(link = logit),data=rbind(A,B))
#---
A=Glist[[3]]
B=rbind(Glist[[4]])
A$ndegree=1
B$ndegree=0
mol.3=glm(ndegree~FFMC+DMC+DC+ISI,
        family = binomial(link = logit),data=rbind(A,B))

y=c()
for(k in 1:nrow(dd)){
  logit=predict(mol.1,data.frame(FFMC=dd$FFMC[k],DMC=dd$DMC[k],
                               DC=dd$DC[k],ISI=dd$ISI[k]))
  pro=exp(logit)/(1+exp(logit))
  if(pro>0.5037){
    y=c(y,1)
  }else{
    logit=predict(mol.2,data.frame(FFMC=dd$FFMC[k],DMC=dd$DMC[k],
                                 DC=dd$DC[k],ISI=dd$ISI[k]))
    pro=exp(logit)/(1+exp(logit))
    if(pro>0.7982){
      y=c(y,2)
    }else{
      logit=predict(mol.3,data.frame(FFMC=dd$FFMC[k],DMC=dd$DMC[k],
                                   DC=dd$DC[k],ISI=dd$ISI[k]))
      pro=exp(logit)/(1+exp(logit))
      if(pro>0.7249){
        y=c(y,3)
      }else{
      y=c(y,4)
      }
    }
  }
}
dd$D1<-y
ddd=subset(data,D1!=1)
N=rbind(dd,ddd)
t1=table(N$ndegree,N$D1)
t1
sum(diag(prop.table(t1)))


K1=subset(N,D1==1)
table(K1$ndegree,K1$month)
K2=subset(N,D1==2)
table(K2$ndegree,K2$month)
K3=subset(N,D1==3)
table(K3$ndegree,K3$month)
K4=subset(N,D1==4)
table(K4$ndegree,K4$month)
