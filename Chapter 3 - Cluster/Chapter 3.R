## Chapter 3

## 加载数据
X=read.table("clipboard",header=T)
rownames(X)=X[,1]
X=X[,-1]

## ---------------------------------------------------------
## 系统聚类法
## ---------------------------------------------------------
## 函数的封装
H.clust <- function(X,d,m,gmin,gmax,ind){ #输入向量c(数据集X,"距离种类","聚类方法",允许最小类数,允许最大类数,"分类指标")
  X.scaled <- scale(X)
  dis <- dist(X.scaled,method=d)
  hc <- hclust(dis,method=m)
  #寻找最佳聚类数
  library(NbClust)
  devAskNewPage(ask=FALSE)
  nc <- NbClust(X.scaled,distance=d,min.nc=gmin,max.nc=gmax,method=m,index=ind)
  t <- table(nc$Best.n[1,])
  tframe <- as.data.frame(t)
  tframe1 <- tframe[which(tframe$Freq==max(t)),]
  k0=as.numeric(as.character(tframe1$Var1))[1]
  #结果输出
  clusters <- cutree(hc,k=k0)
  par(mfrow=c(1,1))
  barplot(t,xlab="Number of Clusters",ylab="Number of Criteria",main="Clusters Chosen")
  plot(hc,hang=-1)
  rect.hclust(hc,k=k0)
  return(list(clusters=table(clusters),median=aggregate(X,by=list(clusters=clusters),median)))
}


## R型聚类
X1=t(X)
H.clust(X1,"euclidean","average",2,7,c("silhouette","gap"))


## Q型聚类
## 最短距离法
H.clust(X,"euclidean","single",2,7,"all")

## 最长距离法
H.clust(X,"euclidean","complete",2,7,"all")

## 重心法
H.clust(X,"euclidean","centroid",2,7,"all")

## 类平均法
H.clust(X,"euclidean","average",2,7,"all")

## Ward法
H.clust(X,"euclidean","ward.D",2,7,"all")



## ---------------------------------------------------------
## K-Means聚类法
## ---------------------------------------------------------
## 确定聚类数目
library(ggplot2)
library(factoextra)
df <- scale(X)
fviz_nbclust(df, kmeans, method = "wss") + geom_vline(xintercept = 3, linetype = 2)


## 进行聚类
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
## 模糊聚类法
## ---------------------------------------------------------
## 读入数据
Y=read.csv("eg3-7.csv")
Y1=Y
rownames(Y1)=Y1[,1]
Y1=Y1[,-1]
Y1=Y1[,-1]

## 进行聚类
library(cluster)
fannyz=fanny(Y1,3,metric="SqEuclidean")
list("分类分布"=fannyz$clustering,"样本隶属度"=fannyz$membership)

## 提升视觉效果
clusplot(fannyz)
library(e1071)
result1<-cmeans(Y1,3,50)
library(scatterplot3d)
scatterplot3d(result1$membership, color=result1$cluster, type="h", 
              angle=55, scale.y=0.7, pch=16, main="Pertinence")

## 计算正确率
t=table(Y[,2],fannyz$clustering)
sum(diag(prop.table(t)))
