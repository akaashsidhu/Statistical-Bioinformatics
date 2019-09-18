This is a collaboration with:
#Akaash Sidhu
#Bram Ratz
#Muskan Aggarwal

#*********************************************************************
#Analysis of Genotype Data Available from the 1000 Genomes Project
#*********************************************************************

#Install packages, load libraries, and read in the OCA file. 

source("https://bioconductor.org/biocLite.R")
#biocLite("trio")
library(tidyverse)
library(dplyr)
library(ggplot2)
library(vcfR)
library(adegenet)
library(adegraphics)
library(trio)
library(lattice)
library(ape)
library(reshape2)
library(snpStats)
library(ggpubr)

setwd()

vcf_OAC1 <- read.vcfR("OAC1gene.vcf")
#Read VCF file into R.

table_OAC1 <- read.table("OAC1gene.vcf", stringsAsFactors = FALSE)
#Create a tabular version

vcf_OAC1
head(vcf_OAC1)
#Check the vcfr object

head(getFIX(vcf_OAC1))
#look at the fix region - contains info about each variantin the sample

vcf_OAC1@gt[1:6, 1:4]
#look at the gt region - contains info genotype info about each variant in the sample

chrom <- create.chromR(name='OAC1_data', vcf=vcf_OAC1)
plot(chrom) 
#This will show important statistics summed over the entire vcf file. 

dp_vcf <- extract.gt(vcf_OAC1, element = 'DP', as.numeric = TRUE)
#Extract allele depths per each sample. This is a quick check for read depth distribution per individual. 

OAC1.genlight <- vcfR2genlight(vcf_OAC1, n.cores = 1)
#Convert into vcfr into genlight object

excel.data <- read.csv("sample.csv")
#Read in the sample excel file that contains the sample ids.

excel.sub <- as.character(excel.data[,1])
names.sub <- indNames(OAC1.genlight)
#Subset the data

SameSamples <- which(excel.sub %in% names.sub)
#Compare the two dataframes to see what is present in both

excel.final <- excel.data[SameSamples,]
#Include only the similar population ids

head(excel.final)
pop(OAC1.genlight) <- excel.final$Population.code
#Add population names

OAC1.genlight
#Check the basic information of the object to make sure everything is working correctly. 

indNames(OAC1.genlight)
#Check individual names

as.matrix(OAC1.genlight)[1:16, 1:10]
#Observe a small portion of the data

pop(OAC1.genlight)
#Look at the populations associated with the sample

#****************************************************************
#                   ~~~~PCA~~~~
#****************************************************************

pca.OAC1 <- glPca(OAC1.genlight, nf = 300, n.cores=1)
#Perform PCA

#Proportion of variance explained by the first three axes

#First axis
pca.OAC1$eig[1]/sum(pca.OAC1$eig) 

#Second axis
pca.OAC1$eig[2]/sum(pca.OAC1$eig) 

#Third axis
pca.OAC1$eig[3]/sum(pca.OAC1$eig) 
#First three explain ~70% of the variance found

pdf("PCA_all_SNPs_ax12.pdf", width = 14, height = 7)
col <- funky(5)

g1 <- s.class(pca.OAC1$scores, 
              pop(OAC1.genlight), 
              xax = 1,
              yax = 2,
              col = transp(col,.6),
              ellipseSize = 0,
              starSize = 0,
              ppoints.cex = 4,
              pgrid.draw = F,
              plot = FALSE)

g2 <- s.label(pca.OAC1$scores,
              xax = 1,
              yax = 2,
              ppoints.col = "red",
              plabels = list(box = list(draw = FALSE), optim = TRUE), paxes.draw = T, pgrid.draw = F, plabels.cex = 1, plot = FALSE)
ADEgS(c(g1, g2), layout = c(1, 2))
dev.off()

pca.OAC1.scores <- as.data.frame(pca.OAC1$scores)
pca.OAC1.scores$pop <- pop(OAC1.genlight)
#An alternative way to obtain PCA

cols <- colours()[500:526]
#Colour variable

p1 <- ggplot(pca.OAC1.scores, aes(x=PC1, 
                                  y=PC2, 
                                  colour=pop(OAC1.genlight)))
p1 <- p1 + geom_point(size=2)
p1 <- p1 + scale_color_manual(values = cols)
p1 <- p1 + theme_bw()
p1
#Here PC1 is plotted against PC2

#****************************************************************
#                   ~~~~K-means Clustering~~~~
#****************************************************************

head(vcf_OAC1)
head(excel.final)

OAC.K.genlight <- vcfR2genlight(vcf_OAC1)
#Create a genlight object

maxK <- 5
#Set maximum number to 5

myMat <- matrix(nrow=5, ncol=maxK)
#Create a 5 by 5 matrix

colnames(myMat) <- 1:ncol(myMat)
#Set column names

for(i in 1:nrow(myMat)){
  grp <- find.clusters(OAC.K.genlight, n.pca = 5, choose.n.clust = FALSE,  max.n.clust = maxK)
  myMat[i,] <- grp$Kstat
}
#Created a for loop to...

my_df <- melt(myMat)
#In order to visualize the clusters, a new dataframe is created.

colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)
#Set column names

head(my_df)
#View the new dataframe

p1 <- ggplot(my_df, aes(x = K, y = BIC))
p1 <- p1 + geom_boxplot()
p1 <- p1 + theme_bw()
p1 <- p1 + xlab("Number of groups (K)")
p1
#Plot the K against BIC

#FOR: DAPC
my_k <- 2:5

grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(9)
  grp_l[[i]] <- find.clusters(OAC.K.genlight, n.pca = 5, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(OAC.K.genlight, pop = grp_l[[i]]$grp, n.pca = 5, n.da = my_k[i])
}

my_df <- as.data.frame(dapc_l[[ length(dapc_l) ]]$ind.coord)
my_df$Group <- dapc_l[[ length(dapc_l) ]]$grp
head(my_df)
#Viewing differentiation of our data into different clusters

my_pal <- RColorBrewer::brewer.pal(n=8, name = "Dark2")

p2 <- ggplot(my_df, aes(x = LD1, y = LD2, color = Group, fill = Group))
p2 <- p2 + geom_point(size = 4, shape = 21)
p2 <- p2 + theme_bw()
p2 <- p2 + scale_color_manual(values=c(my_pal))
p2 <- p2 + scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))
p2
#Plot the data

tmp <- as.data.frame(dapc_l[[1]]$posterior)
tmp$K <- my_k[1]
tmp$Isolate <- rownames(tmp)
tmp <- melt(tmp, id = c("Isolate", "K"))
names(tmp)[3:4] <- c("Group", "Posterior")
tmp$Region <- excel.final$Population.name
my_df <- tmp
#View as a barplot, use facets to separate the different values of K. Combine data into single long form data.frame and add population data

for(i in 2:length(dapc_l)){
  tmp <- as.data.frame(dapc_l[[i]]$posterior)
  tmp$K <- my_k[i]
  tmp$Isolate <- rownames(tmp)
  tmp <- melt(tmp, id = c("Isolate", "K"))
  names(tmp)[3:4] <- c("Group", "Posterior")
  tmp$Region <- excel.final$Population.name
  
  my_df <- rbind(my_df, tmp)
}

grp.labs <- paste("K =", my_k)
names(grp.labs) <- my_k

p3 <- ggplot(my_df, aes(x = Isolate, y = Posterior, fill = Group))
p3 <- p3 + geom_bar(stat = "identity")
p3 <- p3 + facet_grid(K ~ Region, scales = "free_x", space = "free", 
                      labeller = labeller(K = grp.labs))
p3 <- p3 + theme_bw()
p3 <- p3 + ylab("Posterior membership probability")
p3 <- p3 + theme(legend.position='none')
p3 <- p3 + scale_fill_manual(values=c(my_pal))
p3 <- p3 + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p3
#Build the ggplot

ggarrange(ggarrange(p1,
                    p2,
                    ncol = 2, labels = c("A", "B")),
          p3,
          nrow = 2,
          labels = c("", "C"),
          heights = c(1, 2)
)
#Put this all together into one plot using ggpubr



#####

myFreq <- glMean(aa.genlight)
myFreq <- c(myFreq, 1-myFreq)
hist(myFreq, proba=TRUE, col="darkseagreen3", xlab="Allele frequencies",
     main="Distribution of allele frequencies", nclass=20)
temp <- density(myFreq, bw=.05)
lines(temp$x, temp$y*2,lwd=3)
#Generally the allele frequencies are fixed close to 0

temp <- density(glNA(aa.genlight), bw=10)
plot(temp, type="n", xlab="Position in the alignment", main="Location of the missing values (NAs)",
     xlim=c(0,1701))
polygon(c(temp$x,rev(temp$x)), c(temp$y, rep(0,length(temp$x))), col=transp("blue",.3))
points(glNA(flu), rep(0, nLoc(flu)), pch="|", col="blue")
#Here, the few missing values are all located at the beginning at the alignment, probably reflecting heterogeneity in DNA amplification during the sequencing process. In larger datasets, such simple investigation can give crucial insights about the quality of the data and the existence of possible sequencing biases

pca1test <- glPca(aa.genlight)
scatter(pca1test, posi="bottomright")
title("PCA of the OCA1")

myCol <- colorplot(pca1test$scores,pca1test$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")
add.scatter.eig(pca1test$eig[1:40],2,1,2, posi="topright", inset=.05, ratio=.3)


dapc1test <- dapc(aa.genlight, n.pca=5, n.da=1)
scatter(dapc1test,scree.da=FALSE, bg="white", posi.pca="topright", legend=TRUE,
        txt.leg=paste("group", 1:5), col=c("red","blue"))

compoplot(dapc1test, col=c("red","blue", "green", "yellow", "orange"),lab="", txt.leg=paste("group", 1:5), ncol=2)
