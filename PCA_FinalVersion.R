# By Gang Wang, College of Informatics, Huazhong Agricultural University
# 2020-12-25

# ��Ҫ�õ�����R�������δ��װ�����������´���
# install.packages("tidyverse")
# install.packages("data.table")
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("BiocInstaller")
# BiocManager::install("sva")
# BiocManager::install("bladderbatch")
# 
# #PCA��
# install.packages("FactoMineR")
# install.packages("factoextra")

library(tidyverse)
library(data.table)
library(tibble)

rm(list=ls())

# ��ȡ���ݼ�
library(sva)
library(bladderbatch)
data(bladderdata)
# bladderdata��������EsetExpressionSet�����Կ��Էֱ���pData��exprs�������ע����Ϣ�ͱ������
edata <- exprs(bladderEset) # �������
# sum(is.na(edata)) # ���ݼ��

pheno <- pData(bladderEset) # ע����Ϣ


pheno$hasCancer <- as.numeric(pheno$cancer == "Cancer")
pheno[,3]<-as.character(pheno[,3])
for(x in 1:57) if(pheno$hasCancer[x]==1){
  pheno$hasCancer[x]<-'cancar'
} else {
  pheno$hasCancer[x]<-'normal'
}

#ת��edata�������
pca_edata<- as.data.frame(t(edata))

# У��ǰPCA
library(FactoMineR)
library(factoextra)
bd.pca <- PCA(pca_edata, graph = FALSE)
met_tag <- get_pca_var(bd.pca)
eig.val <- get_eigenvalue(bd.pca)
fviz_screeplot(bd.pca,addlabels = T, ylim=c(0,50),ncp = 20)
# PCA��ͼ:У��ǰ����
ind.group <- fviz_pca_ind(bd.pca, geom = "point", addEllipses = TRUE,col.ind = pheno$`hasCancer`)
ggpubr::ggpar(ind.group,
              title = "Principal Component Analysis",
              xlab = "PC1", ylab = "PC2",
              # ��������λ��
              legend.title = "Species", legend.position = "top",
              # �������ɫ
              ggtheme = theme_gray(), palette = "jco")

ind.batch <- fviz_pca_ind(bd.pca, geom = "point", addEllipses = TRUE, col.ind = pheno$`batch`)
ggpubr::ggpar(ind.batch,
              title = "Principal Component Analysis",
              xlab = "PC1", ylab = "PC2",
              # ��������λ��
              legend.title = "batches", legend.position = "top",
              # �������ɫ
              ggtheme = theme_gray(), palette = "jco")


# ��������ЧӦ
model <- model.matrix(~hasCancer, data = pheno)
combat_edata <- ComBat(dat = edata, batch = pheno$batch, mod = model)
combat_edata<- as.data.frame(t(combat_edata))

# ������PCA
combat.pca <- PCA(combat_edata, graph = FALSE)
combat_tag <- get_pca_var(combat.pca)
eig.combat <- get_eigenvalue(combat.pca)
fviz_screeplot(combat.pca,addlabels = T, ylim=c(0,50),ncp = 20)
# PCA��ͼ:У��������
combat.group <- fviz_pca_ind(combat.pca, geom = "point", addEllipses = TRUE,col.ind = pheno$`hasCancer`)
ggpubr::ggpar(combat.group,
              title = "Principal Component Analysis",
              xlab = "PC1", ylab = "PC2",
              # ��������λ��
              legend.title = "Species", legend.position = "top",
              # �������ɫ
              ggtheme = theme_gray(), palette = "jco")

combat.batch <- fviz_pca_ind(combat.pca, geom = "point", addEllipses = TRUE, col.ind = pheno$`batch`)
ggpubr::ggpar(combat.batch,
              title = "Principal Component Analysis",
              xlab = "PC1", ylab = "PC2",
              # ��������λ��
              legend.title = "batches", legend.position = "top",
              # �������ɫ
              ggtheme = theme_gray(), palette = "jco")