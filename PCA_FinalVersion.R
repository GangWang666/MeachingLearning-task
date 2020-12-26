# By Gang Wang, College of Informatics, Huazhong Agricultural University
# 2020-12-25

# 需要用到下列R包，如果未安装，则运行以下代码
# install.packages("tidyverse")
# install.packages("data.table")
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("BiocInstaller")
# BiocManager::install("sva")
# BiocManager::install("bladderbatch")
# 
# #PCA包
# install.packages("FactoMineR")
# install.packages("factoextra")

library(tidyverse)
library(data.table)
library(tibble)

rm(list=ls())

# 获取数据集
library(sva)
library(bladderbatch)
data(bladderdata)
# bladderdata的属性是EsetExpressionSet，所以可以分别用pData和exprs方法获得注释信息和表达矩阵
edata <- exprs(bladderEset) # 表达矩阵
# sum(is.na(edata)) # 数据检查

pheno <- pData(bladderEset) # 注释信息


pheno$hasCancer <- as.numeric(pheno$cancer == "Cancer")
pheno[,3]<-as.character(pheno[,3])
for(x in 1:57) if(pheno$hasCancer[x]==1){
  pheno$hasCancer[x]<-'cancar'
} else {
  pheno$hasCancer[x]<-'normal'
}

#转置edata表达矩阵
pca_edata<- as.data.frame(t(edata))

# 校正前PCA
library(FactoMineR)
library(factoextra)
bd.pca <- PCA(pca_edata, graph = FALSE)
met_tag <- get_pca_var(bd.pca)
eig.val <- get_eigenvalue(bd.pca)
fviz_screeplot(bd.pca,addlabels = T, ylim=c(0,50),ncp = 20)
# PCA作图:校正前数据
ind.group <- fviz_pca_ind(bd.pca, geom = "point", addEllipses = TRUE,col.ind = pheno$`hasCancer`)
ggpubr::ggpar(ind.group,
              title = "Principal Component Analysis",
              xlab = "PC1", ylab = "PC2",
              # 标题名字位置
              legend.title = "Species", legend.position = "top",
              # 主题和配色
              ggtheme = theme_gray(), palette = "jco")

ind.batch <- fviz_pca_ind(bd.pca, geom = "point", addEllipses = TRUE, col.ind = pheno$`batch`)
ggpubr::ggpar(ind.batch,
              title = "Principal Component Analysis",
              xlab = "PC1", ylab = "PC2",
              # 标题名字位置
              legend.title = "batches", legend.position = "top",
              # 主题和配色
              ggtheme = theme_gray(), palette = "jco")


# 较正批次效应
model <- model.matrix(~hasCancer, data = pheno)
combat_edata <- ComBat(dat = edata, batch = pheno$batch, mod = model)
combat_edata<- as.data.frame(t(combat_edata))

# 较正后PCA
combat.pca <- PCA(combat_edata, graph = FALSE)
combat_tag <- get_pca_var(combat.pca)
eig.combat <- get_eigenvalue(combat.pca)
fviz_screeplot(combat.pca,addlabels = T, ylim=c(0,50),ncp = 20)
# PCA作图:校正后数据
combat.group <- fviz_pca_ind(combat.pca, geom = "point", addEllipses = TRUE,col.ind = pheno$`hasCancer`)
ggpubr::ggpar(combat.group,
              title = "Principal Component Analysis",
              xlab = "PC1", ylab = "PC2",
              # 标题名字位置
              legend.title = "Species", legend.position = "top",
              # 主题和配色
              ggtheme = theme_gray(), palette = "jco")

combat.batch <- fviz_pca_ind(combat.pca, geom = "point", addEllipses = TRUE, col.ind = pheno$`batch`)
ggpubr::ggpar(combat.batch,
              title = "Principal Component Analysis",
              xlab = "PC1", ylab = "PC2",
              # 标题名字位置
              legend.title = "batches", legend.position = "top",
              # 主题和配色
              ggtheme = theme_gray(), palette = "jco")
