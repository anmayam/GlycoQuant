# Program to parse an expression file and change it into an expression 

setup = function(file){
  data <- read.csv(file, header=TRUE)
  #file<- "/Users/stan/Desktop/SummerLink/SummerLink-2014-type-4/#NFL_SY5Y.csv"
  row<- dim(data)[2]-1
  col <- dim(data)[1]
  test <- matrix(nrow=row, ncol=col)
  j = 1;
  k = 0;
  for(j in 2:col){
    for(k in 1:row){
      test[k,j-1]<- as.numeric(data[j-1,k+1])
    }}
  l= 0; 
  for(l in 1:row){
    test[l,col] <- as.numeric(data[col,l+1])
  }
  return(test)

