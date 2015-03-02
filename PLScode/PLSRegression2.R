# Anoop Mayampurath
# program to perform PLS regression on glycan abundance
library("plsdepot") ; 
data.expression <-  read.csv("C:/development/glycan_graphs/data/NLF_SY5Y_Expression_Labels.csv", header = TRUE) ;
data.profiles <- read.csv("C:/development/glycan_graphs/data/NLF_SY5Y_Subtrees_Labels.csv", header = TRUE) ; 

g.fc <- as.numeric(data.expression$FoldChange) ;
g.profiles.composition <- as.data.frame(data.profiles[, 3:6]) ; 
g.profiles.seq <- as.data.frame(data.profiles[, 7:12]);
g.profiles.sub <- as.data.frame(data.profiles[, 15:26]); 

y <- array(log2(1/g.fc)) ; 


model.comp <- plsreg1(g.profiles.composition, y, comps = 4) ; 
model.seq <- plsreg1(g.profiles.seq, y, comps = 3) ; 
model.sub <- plsreg1(g.profiles.sub, y, comps = 6) ; 


