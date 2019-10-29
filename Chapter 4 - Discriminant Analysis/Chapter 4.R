## Chapter 4

## --------------------------------------------------
## ��4.1, Fisher�б�
## --------------------------------------------------
## ��������
X = read.csv("eg4-1.csv",header = T)

## ������ͳ��
attach(X)
aggregate(cbind(Sepal.Length,Sepal.Width,Petal.Length,Petal.Width),by=list(y),FUN=mean)

X.1 <- as.data.frame(cbind(Sepal.Length,Sepal.Width,Petal.Length,Petal.Width,y))
type <- as.factor(X.1$y)
fit <- manova(cbind(Sepal.Length,Sepal.Width,Petal.Length,Petal.Width)~type,data=X.1)
summary(fit,test=c("Wilks"))

## ��ʽ�������
## ����ģ��
library(MASS)
ld=lda(y~.,X)
ld
ld$svd

## Ԥ��ͻ���
X.pre=predict(ld)
newy=X.pre$class
head(cbind(X,newy))
tab=table(y,newy)
tab
sum(diag(prop.table(tab)))

## ������ͼ
library(ggplot2)
y <- as.factor(X$y)
X.2=cbind(y,as.data.frame(X.pre$x))
p=ggplot((X.2),aes(x=LD1,y=LD2))
p+geom_point(aes(colour=y),alpha=0.8,size=4)

## Ԥ��
pre=predict(ld,data.frame(Sepal.Length=5.8,Sepal.Width=3.1,Petal.Length=3.8,Petal.Width=1.2))
pre
## Ԥ����ͼ
pre.1=cbind("4",as.data.frame(pre$x))
colnames(pre.1)=colnames(X.2)
y<-factor(y,levels = c(1,2,3,4))
pre.2=rbind(X.2[1:nrow(X.2),],pre.1)
p=ggplot(cbind(pre.2),aes(x=LD1,y=LD2))
p+geom_point(aes(colour=y),alpha=0.8,size=4)



## --------------------------------------------------
## ��4.2, �����б�
## --------------------------------------------------
## ��������
Z = read.csv("eg4-2.csv",header = T)
rownames(Z)=Z[,2]
Z=Z[,-1]
Z=Z[,-1]

## ���鷽������
attach(Z)
r=3
G1 <- subset(Z,Group=="1")
G2 <- subset(Z,Group=="2")
G3 <- subset(Z,Group=="3")
n1=nrow(G1)
n2=nrow(G2)
n3=nrow(G3)
n=n1+n2+n3

L1=(n1-1)*cov(G1[,1:8])
L2=(n2-1)*cov(G2[,1:8])
L3=(n3-1)*cov(G3[,1:8])
L=L1+L2+L3
Sp=L/(n-r)
Sp_1=solve(Sp)

m1=colMeans(G1[,1:8])
m2=colMeans(G2[,1:8])
m3=colMeans(G3[,1:8])

## �������Ͼ���  
d1=c()
for (i in 1:nrow(Z)){
  dis1=as.matrix(Z[i,1:8]-m1)%*%as.matrix(Sp_1)%*%as.matrix(t(Z[i,1:8]-m1))
  d1=c(d1,dis1)
}

d2=c()
for (i in 1:nrow(Z)){
  dis2=as.matrix(Z[i,1:8]-m2)%*%as.matrix(Sp_1)%*%as.matrix(t(Z[i,1:8]-m2))
  d2=c(d2,dis2)
}

d3=c()
for (i in 1:nrow(Z)){
  dis3=as.matrix(Z[i,1:8]-m3)%*%as.matrix(Sp_1)%*%as.matrix(t(Z[i,1:8]-m3))
  d3=c(d3,dis3)
}

## ����Ԥ�⼰����
newG=c()
dist=c()
for(z in 1:nrow(Z)){
  if(d1[z]==min(d1[z],d2[z],d3[z])){
    newG=c(newG,1)
  }
  if(d2[z]==min(d1[z],d2[z],d3[z])){
    newG=c(newG,2)
  }
  if(d3[z]==min(d1[z],d2[z],d3[z])){
    newG=c(newG,3)
  }
  dist=c(dist,min(d1[z],d2[z],d3[z]))
}

## ������
Z.1=data.frame(Z$Group,dist,newG)
rownames(Z.1)=rownames(Z)
colnames(Z.1)=c("ԭ���","����","�����")
Z.1

## ������ȷ��
Z.pre=cbind(Z,newG)
tab=table(Z.pre$Group[1:n],Z.pre$newG[1:n])
tab

sum(diag(prop.table(tab)))



## --------------------------------------------------
## ��4.3, Fisher�б�
## --------------------------------------------------
## ��������
Y = read.csv("eg4-3.csv",header = T)
rownames(Y)=Y[,2]
Y=Y[,-1]
Y=Y[,-1]

## ɸѡ������G1��G2
attach(Y)
G1 <- subset(Y,Group=="1")
G2 <- subset(Y,Group=="2")
n1=nrow(G1)
n2=nrow(G2)
n=n1+n2
p=8
r=2

## (1)��������G1��G2���б�����ľ�ֵ
## �����ֵ
x1 <- colMeans(G1[,1:8])
x2 <- colMeans(G2[,1:8])
x1-x2
x1+x2

## (2)����Э������Sigma�Ĺ���ֵ����
S1=cov(G1[,1:8])
S2=cov(G2[,1:8])
Sp=(S1*(n1-1)+S2*(n2-1))/(n1+n2-2)
Sp_1=solve(Sp)


## (3)����Fisher�����б���
alpha=t(x1-x2)%*%Sp_1

## (4)����������Ԫ�����ֵ���е�m�Ĺ���ֵ
m=0.5*alpha%*%(x1+x2)

## (5)����ͳ�Ƽ�����Fֵ
## ���Ͼ���Ϊ
D2=t(x1-x2)%*%Sp_1%*%(x1-x2)
## Fͳ����Ϊ
F=(n1+n2-p-1)/((n1+n2-2)*p)*n1*n2/(n1+n2)*D2
pvalue=pf(F,p,n1+n2-p-1,ncp=0,lower.tail=FALSE,log.p=FALSE)

## (6)���м�������Ʒ�Ĺ���
T=alpha%*%t(Y[,1:8])
newG=c()
for(i in 1:31){
  if(T[i]>=m){
    newG=c(newG,1)
  }
  if(T[i]<m){
    newG=c(newG,2)
  }
}

## �õ����
Y.1=data.frame(Y$Group,t(T),newG)
rownames(Y.1)=rownames(Y)
colnames(Y.1)=c("ԭ���","�б���ֵ","�����")
Y.1

## (7)������ȷ��
tab=table(Y$Group[n1+n2],newG[n1+n2])
tab
sum(diag(prop.table(tab)))

## (8)��ͼ
#������ͼ
library(ggplot2)
type <- as.factor(newG)
Y.2=cbind(type,as.data.frame(t(T)))
p<-ggplot(Y.2, aes(x = factor(rownames(Y.2)), fill = factor(Y.2$type), y = Y.2$V1))
q<-geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge")
u<-geom_hline(aes(yintercept = m),color="gray40"  ,linetype="dashed",cex=1)
p+q+u+theme(axis.text.x = element_text(angle = 90))



## --------------------------------------------------
## ��4.1, Bayes�б�
## --------------------------------------------------
## lda()
X = read.csv("eg4-1.csv",header = T)
attach(X)
ld1=lda(y~.,prior=c(1,2,1)/4,X)
pre.by=predict(ld1)
by.1<-data.frame(X$y,pre.by$posterior,pre.by$class)
colnames(by.1)=c("ԭ���","�������G1","�������G2","�������G3","�����")
by.1[c(1:5,51:55,101:105),]

## ����Э��������Ƿ����
r=3
G1 <- subset(X,y=="1")
G2 <- subset(X,y=="2")
G3 <- subset(X,y=="3")
n1=nrow(G1)
n2=nrow(G2)
n3=nrow(G3)
n=n1+n2+n3

source("cov_test.R")
cov.test(4,c(1:4))

## ��д�������Bayes�б�
discriminiant.bayes<-function(TrnX1, TrnX2, rate=1, TstX = NULL, var.equal = FALSE){
  if (is.null(TstX) == TRUE) TstX<-rbind(TrnX1,TrnX2)
  if (is.vector(TstX) == TRUE)  TstX<-t(as.matrix(TstX))
  else if (is.matrix(TstX) != TRUE)
    TstX<-as.matrix(TstX)
  if (is.matrix(TrnX1) != TRUE) TrnX1<-as.matrix(TrnX1)
  if (is.matrix(TrnX2) != TRUE) TrnX2<-as.matrix(TrnX2)
  
  nx<-nrow(TstX)
  blong<-matrix(rep(0, nx), nrow=1, byrow=TRUE, 
                dimnames=list("blong", 1:nx))
  mu1<-colMeans(TrnX1); mu2<-colMeans(TrnX2) 
  if (var.equal == TRUE  || var.equal == T){
    S<-var(rbind(TrnX1,TrnX2)); beta<-2*log(rate)
    w<-mahalanobis(TstX, mu2, S)-mahalanobis(TstX, mu1, S)
  }
  else{
    S1<-var(TrnX1); S2<-var(TrnX2)
    beta<-2*log(rate)+log(det(S1)/det(S2))
    w<-mahalanobis(TstX, mu2, S2)-mahalanobis(TstX, mu1, S1)
  }
  
  for (i in 1:nx){
    if (w[i]>beta)
      blong[i]<-1
    else
      blong[i]<-2
  }
  blong
}


d.b1<-discriminiant.bayes(G1[,1:4],G2[,1:4],rate=1,var.equal = FALSE)
d.b1
table(X$y[1:100],d.b1)
sum(diag(prop.table(table(X$y[1:100],d.b1))))

d.b2<-discriminiant.bayes(G1[,1:4],G3[,1:4],rate=1,var.equal = FALSE)
d.b2
table(X$y[1:100],d.b2)
sum(diag(prop.table(table(X$y[1:100],d.b2))))

d.b3<-discriminiant.bayes(G2[,1:4],G3[,1:4],rate=1,var.equal = FALSE)
d.b3
table(X$y[1:100],d.b3)
sum(diag(prop.table(table(X$y[1:100],d.b3))))