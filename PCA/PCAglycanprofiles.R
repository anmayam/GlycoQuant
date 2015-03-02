# Anoop 
# program to perform PCA on glycomic profiles

#rm(list =ls()); 
library(ggplot2)
library(MASS)



data.in <- read.csv("C:/development/glycan_graphs/data/NLF_SY5Y_Expression_Labels.csv", header = TRUE) ;
#data.in <- read.csv("C:/development/glycan_graphs/data/NLF_SY5Y_Composition_Labels.csv", header = TRUE) ; 
num.classes <- c(length(which(data.in$Labels ==-1)), length(which(data.in$Labels==0)), length(which(data.in$Labels==1))) 
colors = c(rep("red", 1, num.classes[1]), rep("gray", 1, num.classes[2]), rep("blue", 1, num.classes[3]));
#g.profiles <- (data.in[,3:6]);
#g.profiles <- (data.in[,7:12]);
g.profiles <- (data.in[,2:4]) ; 

N <- dim(g.profiles)[1] ; 
k <- dim(g.profiles)[2];

g.dist <- mds.edm1(g.profiles); 

g.cmds <- cmdscale(g.dist, k =2);
plot(g.cmds[,1], g.cmds[,2],  col = colors, cex = 1.5, pch = 16, main = "CMDS plot", xlab = "g.cmds[,1]", ylab = "g.cmds[,2]")



browser();
g.pca <- prcomp(g.profiles);
scores = as.data.frame(g.pca$x) ; 
plot.new()

ggplot(data = scores, aes(x = PC1, y = PC2,  label = rownames(scores)))+ 
 geom_hline(yintercept = 0, colour = "gray65") +
 geom_vline(xintercept = 0, colour = "gray65") +
 geom_text(colour = colors, alpha = 0.8, size = 6) +
 ggtitle("PCA plot")



