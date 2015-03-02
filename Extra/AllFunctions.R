mds.edm1 <- function(X) {
#
#  Computes an EDM-1 from a data matrix.
#
return(as.matrix(dist(X)))
}


graph.knn <- function(D,k) {
#
#  Returns an nxn 0-1 matrix KNN.
#    KNN[i,j]=1 iff j is a neighbor of i.
#    Note that KNN may be asymmetric.
#  We assume that D[i,j]>0 except when i=j.
#
n <- nrow(D)
KNN <- matrix(0,nrow=n,ncol=n)
near <- 2:(k+1)
for (i in 1:n) {
  v <- D[i,]
  j <- order(v)
  j <- j[near]
  KNN[i,j] <- 1
  }
return(KNN)
}


graph.adj <- function(KNN) {
#
#  Uses the output of graph.knn to construct an adjacency matrix.
#  Vertices i & j are connected iff either j is a neighbor of i
#    or i is a neighbor of j.
#
return(pmax(KNN,t(KNN)))
}


graph.unit <- function(KNN) {
#
#  Assigns unit edge weights based on adjacency.
#  W[i,j]=1 if i & j are connected, W[i,i]=0, otherwise W[i,j] = Inf
#
n <- nrow(KNN)
A <- graph.adj(KNN)
i <- which(A==0)
A[i] <- Inf
i <- seq(from=1,to=n^2,by=n+1)
A[i] <- 0
return(A)
}


graph.dis <- function(KNN,D) {
#
#  Assigns dissimilarity edge weights based on adjacency.
#  W[i,j] = D[i,j] if i & j are connected, W[i,i]=0, 
#    otherwise W[i,j] = Inf.
#
n <- nrow(KNN)
A <- graph.adj(KNN)
i <- which(A==0)
D[i] <- Inf
i <- seq(from=1,to=n^2,by=n+1)
D[i] <- 0
return(D)
}


graph.heat <- function(D,sigma) {
#
#  Assigns edge weights equal to similarities based on the heat kernel.
#
return(exp(-sigma*D^2))
}


graph.short <- function(W) {
#
#  Computes all shortest path distances for a weighted graph.
#
n <- nrow(W)
E <- matrix(0,nrow=n,ncol=n)
m <- 1
while (m < n-1) {
  for (i in 1:n) {
    for (j in 1:n) {
      E[i,j] <- min(W[i,]+W[,j])
      }
    }
  W <- E
  m <- 2*m
  }
return(W)
}


graph.laplacian <- function(Gamma) {
#
#  Computes the Laplacian matrix of a weighted graph.
#
tot <- apply(Gamma,1,sum)
return(diag(tot)-Gamma)
}


mds.tau <- function(H)
{
#
#  This function returns the double centering of the inputted matrix.
#  See Critchley for details.
#
        n <- nrow(H)
        P <- diag(n) - 1/n
        return(-0.5 * P %*% H %*% P)
}


