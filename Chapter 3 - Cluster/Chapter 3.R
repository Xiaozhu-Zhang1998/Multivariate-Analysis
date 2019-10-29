## ��������
X=read.table("clipboard",header=T)
rownames(X)=X[,1]
X=X[,-1]

## ---------------------------------------------------------
## ϵͳ���෨
## ---------------------------------------------------------
## �����ķ�װ
H.clust <- function(X,d,m,gmin,gmax,ind){ #��������c(���ݼ�X,"��������","���෽��",������С����,�����������,"����ָ��")
  X.scaled <- scale(X)
  dis <- dist(X.scaled,method=d)
  hc <- hclust(dis,method=m)
  #Ѱ����Ѿ�����
  library(NbClust)
  devAskNewPage(ask=FALSE)
  nc <- NbClust(X.scaled,distance=d,min.nc=gmin,max.nc=gmax,method=m,index=ind)
  t <- table(nc$Best.n[1,])
  tframe <- as.data.frame(t)
  tframe1 <- tframe[which(tframe$Freq==max(t)),]
  k0=as.numeric(as.character(tframe1$Var1))[1]
  #������
  clusters <- cutree(hc,k=k0)
  par(mfrow=c(1,1))
  barplot(t,xlab="Number of Clusters",ylab="Number of Criteria",main="Clusters Chosen")
  plot(hc,hang=-1)
  rect.hclust(hc,k=k0)
  return(list(clusters=table(clusters),median=aggregate(X,by=list(clusters=clusters),median)))
}


## R�;���
X1=t(X)
H.clust(X1,"euclidean","average",2,7,c("silhouette","gap"))


## Q�;���
## ��̾��뷨
H.clust(X,"euclidean","single",2,7,"all")

## ����뷨
H.clust(X,"euclidean","complete",2,7,"all")

## ���ķ�
H.clust(X,"euclidean","centroid",2,7,"all")

## ��ƽ����
H.clust(X,"euclidean","average",2,7,"all")

## Ward��
H.clust(X,"euclidean","ward.D",2,7,"all")



## ---------------------------------------------------------
## K-Means���෨
## ---------------------------------------------------------
## ȷ��������Ŀ
library(ggplot2)
library(factoextra)
df <- scale(X)
fviz_nbclust(df, kmeans, method = "wss") + geom_vline(xintercept = 3, linetype = 2)


## ���о���
set.seed(123)
km <- kmeans(df,3)
km$centers

km$iter

dd <- cbind(X, cluster = km$cluster)
table(dd$cluster)

fviz_cluster(km, data = df,
             palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
             ellipse.type = "euclid",
             star.plot = TRUE, 
             repel = TRUE,
             ggtheme = theme_minimal()
)             



## ---------------------------------------------------------
## ģ�����෨
## ---------------------------------------------------------
## ��������
Y=read.csv("eg3-7.csv")
Y1=Y
rownames(Y1)=Y1[,1]
Y1=Y1[,-1]
Y1=Y1[,-1]

## ���о���
library(cluster)
fannyz=fanny(Y1,3,metric="SqEuclidean")
list("����ֲ�"=fannyz$clustering,"����������"=fannyz$membership)

## �����Ӿ�Ч��
clusplot(fannyz)
library(e1071)
result1<-cmeans(Y1,3,50)
library(scatterplot3d)
scatterplot3d(result1$membership, color=result1$cluster, type="h", 
              angle=55, scale.y=0.7, pch=16, main="Pertinence")

## ������ȷ��
t=table(Y[,2],fannyz$clustering)
sum(diag(prop.table(t)))