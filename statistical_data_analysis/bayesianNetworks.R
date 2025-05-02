################################################################# 
#####            Project1: Bayesian Networks                #####
#################################################################

#
# 0) Setup environment
# 
### Install the pakages "deal" and "yeastCC"
install.packages("deal",repos="http://lib.stat.cmu.edu/R/CRAN",dependencies=TRUE)
BiocManager::install("yeastCC")

### Load the packages
library("deal")
library("yeastCC")
library(network)
#
# 1) Data Preprocessing
#
### Load the data
expr <- as.data.frame(t(exprs(yeastCC)[orf800,]))
cat("Observations:", nrow(expr), "\n")
cat("Genes:", ncol(expr), "\n")

### Replace missing data with gene expression median
# for-loop
for (j in 1:ncol(expr)) {
  index_na <- which(is.na(expr[,j]))  # indices of NAs
  expr[index_na, j] <- median(expr[,j], na.rm=T)  # replace by median
}

# apply() alternative 
medianfill <- function(exprCol) {
  exprCol[which(is.na(exprCol))] <- median(exprCol, na.rm=T)
  return(exprCol)
}
expr <- as.data.frame(apply(expr, 2, medianfill))

### Filter genes based on inter quartile range iqr
# for-loop
iqr <- vector(length=ncol(expr))
for (j in 1:ncol(expr)) {
  iqr[j] <- quantile(expr[,j], c(0.75)) - quantile(expr[,j], c(0.25))
}

# apply() alternative
gene.iqr <- function(exprCol) {
  quantile(exprCol, c(0.75)) - quantile(exprCol, c(0.25))
}
iqr <- apply(expr, 2, gene.iqr)

# Keep only genes with iqr > 1.6
expr = expr[, iqr > 1.6]   # keep only genes with large variation
names(expr)  # print selected genes
#YBR054W	YRO2
#YBR088C	POL30
#YER124C	DSE1
#YGL028C	SCW11
#YLR286C	CTS1
#YHR143W	DSE2
#YNL327W	EGT2
#YGR108W	CLB1
#YNR067C	DSE4
#YOL007C	CSI2

# deal package seems not to work if all nodes are continuous.
# Thus, we need to add a discrete "dummy" node to work around this issue. 
# Adding this node will not have any influence on the final results.
expr$dummy <- factor(rep("1", nrow(expr)))
head(expr, 3)

#
# 2-6) Build Bayesian Network
#
### 2) Create prior structure 
# EITHER: Using default 
#G0 <- network(expr) # prior structure with no edges

# OR: Specify your own prior network manually:
# To insert an arrow from node 'A' to node 'B',
# first click node 'A' and then click node 'B'.
# Add an arrow from YOL007C to YBR088C, 
# and from YNL327W to YER124C, YHR143W and YNR067C.
# When the graph is finished, click 'stop',
# Then, inspect the local probability distribution of node i,
# by clicking on node i and then on the background.,
# Finish by clicking the center button 'Stop'.
G0  <- network(expr, specifygraph=TRUE, inspectprob=TRUE)

# We don't want any arrows starting from the "dummy" node, thus we construct a list of banned dependencies:
banlist(G0) <- matrix(c(11,11,11,11,11,11,11,11,11,11,1,2,3,4,5,6,7,8,9,10),ncol=2)
plot(G0)


### 3) Show local probability distribution
localprob(G0)
localprob(G0)$YBR088C

### 4) Compute joint prior distribution
prior <- jointprior(G0, 5)  # equivalent to imaginary sample size = 5

### 5) Learn the initial network

G0 <- getnetwork(learn(G0, expr, prior))
print(G0$score)

### 6) Search for optimal network (takes some time)
nwSearch <- autosearch(G0, expr, prior, removecycles=FALSE, trace=FALSE)
G <- getnetwork(nwSearch)
plot(G)
# TU SKOŃCZYŁEM


sigmas = apply(expr, MARGIN = 2, var)

dane = replicate(30, expr, simplify= FALSE)

for (i in 1:length(dane)){
  for (j in 1:(length(dane[[1]][1,])-1)){
    dane[[i]][,j] = dane[[i]][,j] + rnorm(77, mean = 0, sd = sqrt(sigmas[j]/10))
  }
}


gen1=sapply(dane, function(x) x[,6])
gen1 = t(gen1)
#dim(gen1) = c(30,77)
boxplot(gen1)


pbns = lapply(dane, build.optimal.network)

plot.bn(pbns[[5]], file="/home/jasiekst/fuw/sadII/pbn5.pdf")

genes=head(colnames(expr),10)

build.adjacency.matrix <- function(exprData) { 
  M = matrix(numeric(100), nrow=10)
  p = exprData$nodes$YBR054W$parents
  for (i in 1:length(p)){
    M[p[i],1] = 1
  }
  p = exprData$nodes$YBR088C$parents
  for (i in 1:length(p)){
    M[p[i],2] = 1
  }
  p = exprData$nodes$YER124C$parents
  for (i in 1:length(p)){
    M[p[i],3] = 1
  }
  p = exprData$nodes$YGL028C$parents
  for (i in 1:length(p)){
    M[p[i],4] = 1
  }
  p = exprData$nodes$YGR108W$parents
  for (i in 1:length(p)){
    M[p[i],5] = 1
  }
  p = exprData$nodes$YHR143W$parents
  for (i in 1:length(p)){
    M[p[i],6] = 1
  }
  p = exprData$nodes$YLR286C$parents
  for (i in 1:length(p)){
    M[p[i],7] = 1
  }
  p = exprData$nodes$YNL327W$parents
  for (i in 1:length(p)){
    M[p[i],8] = 1
  }
  p = exprData$nodes$YNR067C$parents
  for (i in 1:length(p)){
    M[p[i],9] = 1
  }
  p = exprData$nodes$YOL007C$parents
  for (i in 1:length(p)){
    M[p[i],10] = 1
  }
  return(M)
}

pbns.adj = lapply(pbns, build.adjacency.matrix)
summaric.edges = Reduce('+', pbns.adj)/30
idx = which(summaric.edges>0)
summaric.edges[idx]

g1 = build.adjacency.matrix(G)
idx.present = which(g1==1)
idx.absent = which(g1==0)

present.final = intersect(idx.present, idx)
absent.final = intersect(idx.absent, idx)


build.ticks <- function(a){
  chart.ticks = c()
  for (i in 1:length(a)){
    parent.name = genes[(a[i]%/%10)+1]
    child.name = genes[a[i]%%10+1]
    s = paste(parent.name, '->', child.name)
    chart.ticks = c(chart.ticks, s)
  }
  return(chart.ticks)
}

t1 = build.ticks(present.final)
t2 = build.ticks(absent.final)

par(mar=c(8,4,2,2))
plot(summaric.edges[present.final], xaxt='n', xlab='', ylab='edge frequency')
axis(1, at=1:length(t1), labels=t1, las=2, padj=1, cex.lab=0.5, cex.axis=0.775)

par(mar=c(8,4,2,2))
plot(summaric.edges[absent.final], xaxt='n', xlab='', ylab='edge frequency')
axis(1, at=1:length(t2), labels=t2, las=2, padj=1, cex.lab=0.5, cex.axis=0.775)

q1 = quantile(summaric.edges[present.final], c(0.25))
edges1 = build.ticks(intersect(which(summaric.edges<q1),present.final))

q2 = quantile(summaric.edges[absent.final], c(0.75))
edges2 = build.ticks(intersect(which(summaric.edges>q2),absent.final))


summaric.edges[]


g = c(1,2)
expr[1,c]



build.adjacency.matrix(pbns[[1]])

for (el in pbns[[1]]$nodes){
  print(pbns[[1]]$nodes$el)
  }






g1 = get.edges(pbns[[5]], 1, neighborhood = ("combined"))


get.adjacency(graph = pbns[[1]])


pbn1=build.optimal.network(dane[[1]])


rnorm 


### Function for building an optimal network from expression data (from above)
build.optimal.network <- function(exprData, N=5) { 
  #N0 <- network(exprData)
  #print(head(exprData))
  #banlist(N0) <- matrix(c(11,11,11,11,11,11,11,11,11,11,1,2,3,4,5,6,7,8,9,10),ncol=2)
  #plot(N0)
  prior <- jointprior(G0, N)
  N0 <- getnetwork(learn(G0, exprData, prior))
  nwasarch <- autosearch(N0, exprData, prior, removecycles=FALSE, trace=FALSE)
  getnetwork(nwasarch)
}

### Custom function for plotting a BN 
plot.bn <- function(BN, file=NULL) {
  par(mar=c(0,0,0,0))
  plot(BN, cexscale=13, unitscale=27, arrowlength=0.1, xr=c(0, 350), yr=c(20,370))
  if (!is.null(file)) {
    plt <- recordPlot()	
    pdf(file)
    replayPlot(plt)
    dev.off()
  }
}
BN <- build.optimal.network(expr)
plot(BN)
plot.bn(BN, file="~/Desktop/BNstar.pdf")
