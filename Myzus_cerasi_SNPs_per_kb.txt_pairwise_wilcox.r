
data<-read.table("Myzus_cerasi_SNPs_per_kb.txt",header=TRUE)
attach(data)
X <- as.numeric(VARIANTS.KB)
C <- as.factor(plant)
dat <- c(X,C)
nam <- c(rep(1,length(X)),rep(2,length(C)))
kruskal.test(dat,nam)
pairwise.wilcox.test(VARIANTS.KB, plant, p.adj="bonferroni", exact=F)