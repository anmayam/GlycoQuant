# Anoop
# Program to plot correlation and heat map



#Correlation
panel.cor <- function(x, y, digits=3, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = (cor(x,y, method = "pearson"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex * abs(r))
}


data.in <- read.csv("C:\\development\\NBL_glycomics\\data\\NFL_SY5Y.csv", header = TRUE, row.names = 1) ; 
g.profiles <- t(data.in[,-1]);
pairs(g.profiles,upper.panel=panel.cor)


#HeatMap
heatmap.2(t(g.profiles), trace="none", density="none", 
          scale="row",
          hclust=function(x) hclust(x,method="complete"),margins=c(5,9),
          distfun=function(x) as.dist((1-cor(t(x)))/2))



