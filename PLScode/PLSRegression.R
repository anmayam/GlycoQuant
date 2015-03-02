# Anoop Mayampurath
#program to perform PLS regression on glycan abundance
library("pls") ; 
data.expression <-  read.csv("C:/development/glycan_graphs/data/NLF_SY5Y_Expression_Labels.csv", header = TRUE) ;
data.profiles <- read.csv("C:/development/glycan_graphs/data/NLF_SY5Y_Subtrees_Labels.csv", header = TRUE) ; 

g.fc <- as.numeric(data.expression$FoldChange) ;
g.profiles.composition <- as.data.frame(data.profiles[, 3:6]) ; 
g.profiles.seq <- as.data.frame(data.profiles[, 7:12]);
g.profiles.sub <- as.data.frame(data.profiles[, 15:26]); 

y <- log2(g.fc) ; 
data.fit.comp <- cbind.data.frame(y, g.profiles.composition);
data.fit.seq <- cbind.data.frame(y, g.profiles.seq) ; 
data.fit.sub <- cbind.data.frame(y, g.profiles.sub) ;

model.comp <- plsr(y ~.,  data = data.fit.comp, valdation = "LOO") ; 
model.seq <- plsr(y ~., data = data.fit.seq, validation = "LOO") ; 
model.sub <- plsr(y ~., data= data.fit.sub, validation = "LOO") ; 

yvar.model.comp <- 100*drop(R2(model.comp, estimate = "train", intercept=FALSE)$val)
yvar.model.seq <- 100*drop(R2(model.seq, estimate = "train", intercept=FALSE)$val)
yvar.model.sub <- 100*drop(R2(model.sub, estimate = "train", intercept=FALSE)$val)

ylimits = c(20,80); 
x_marks <- seq(1, 4, by =1) ; 
plot(x_marks, yvar.model.comp, type = "b", col = "red", lwd = 3,  ylim = ylimits, xlab = "Number of components", ylab = "% explained variance", cex.lab = 1.5, cex.axis = 1.2); 
lines(x_marks, yvar.model.seq[1:4], type = "b", col = "blue", lwd = 3) ;
lines(x_marks, yvar.model.sub[1:4], type = "b", col = "black", lwd = 3) ;
legend("topright", legend =c("Composition", "Monosaccharide", "Sequence"), bty = c(1,1,1), lwd=c(3,3,3), cex = 1.5,  col = c("red", "blue", "black"))
