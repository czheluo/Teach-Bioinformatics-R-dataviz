##########################################
#  http://www.majorbio.com/
#  Copyright (C) 2019 in Majorbio workshop
#  contact: meng.luo@majorbio.com
##########################################

################## Correlation based network analysis ###############################

# #install.packages("vegan")
# install.packages("igraph")
# install.packages("Hmisc")

library(vegan)
library(igraph)
library(Hmisc)


coRnetwork <- function(matrix, cor.cutoff, p.cutoff) {
  # load packages
  library(vegan)
  library(igraph)
  library(Hmisc)
  
  matrix1 <- matrix
  matrix1[matrix1 > 0] <- 1

  # correlation analysis based on spearman's or pearson co-efficient
  matrix.dist <- rcorr(t(matrix), type = "spearman")
  # matrix.dist<-rcorr(t(matrix),type="pearson")
  matrix.cor <- matrix.dist$r
  matrix.cor.p <- matrix.dist$P


  # Multiple testing correction using Benjamini-Hochberg standard false discovery rate correction
  matrix.cor.p <- p.adjust(matrix.cor.p, method = "BH")

  # Consider positive cooccurence at given coefficient (cor.cutoff) and p-value cutoffs
  matrix.cor1 <- matrix.cor
  matrix.cor1.p <- matrix.cor.p
  matrix.cor1[which(matrix.cor1 <= cor.cutoff)] <- 0
  matrix.cor1[which(matrix.cor1.p > p.cutoff)] <- 0
  # delete those rows and columns with sum = 0
  matrix.cor1 <- matrix.cor1[which(rowSums(matrix.cor1) != 1), ]
  matrix.cor1 <- matrix.cor1[, which(colSums(matrix.cor1) != 0)]

  # Consider netagive cooccurence at given coefficient (-cor.cutoff) and p-value cutoffs
  matrix.cor2 <- matrix.cor
  matrix.cor2.p <- matrix.cor.p
  matrix.cor2[which(matrix.cor2 > (-cor.cutoff))] <- 0
  matrix.cor2[which(matrix.cor2.p > p.cutoff)] <- 0
  # delete those rows and columns with sum = 0
  matrix.cor2 <- matrix.cor2[which(rowSums(matrix.cor2) != 0), ]
  matrix.cor2 <- matrix.cor2[, which(colSums(matrix.cor2) != 0)]

  # Consider both positive and netagive cooccurence at given coefficient (cor.cutoff) and p-value cutoffs
  matrix.cor3 <- matrix.cor
  matrix.cor3.p <- matrix.cor.p
  matrix.cor3[which(matrix.cor3 >= (-cor.cutoff) & matrix.cor3 <= cor.cutoff)] <- 0
  matrix.cor3[which(matrix.cor3.p > p.cutoff)] <- 0

  # delete those rows and columns with sum = 0
  matrix.cor3 <- matrix.cor3[which(rowSums(matrix.cor3) != 1), ]
  matrix.cor3 <- matrix.cor3[, which(colSums(matrix.cor3) != 0)]

  # get pairs r
  # This is to remove redundancy as upper correlation matrix == lower
  ma1 <- matrix.cor1
  ma2 <- matrix.cor2
  ma3 <- matrix.cor3
  ma1[upper.tri(matrix.cor1, diag = TRUE)] <- NA
  pair.r1 <- reshape2::melt(ma1, na.rm = TRUE, value.name = "cor")
  ma2[upper.tri(ma2, diag = TRUE)] <- NA
  pair.r2 <- reshape2::melt(ma2, na.rm = TRUE, value.name = "cor")
  ma3[upper.tri(ma3, diag = TRUE)] <- NA
  pair.r3 <- reshape2::melt(ma3, na.rm = TRUE, value.name = "cor")
  pair.r1<-pair.r1[which(pair.r1[,3]!=0),]
  pair.r2<-pair.r2[which(pair.r2[,3]!=0),]
  pair.r3<-pair.r3[which(pair.r3[,3]!=0),]
  write.csv(pair.r1, file = "Pos_otu.csv",row.names = F)
  write.csv(pair.r2, file = "Neg_otu.csv",row.names = F)
  write.csv(pair.r3, file = "PosNeg_otu.csv",row.names = F)

  # generating graph using igraph
  g1 <- graph.adjacency(matrix.cor1, weight = T, mode = "undirected")
  g1 <- simplify(g1)
  V(g1)$label <- V(g1)$name
  V(g1)$degree <- degree(g1)

  g2 <- graph.adjacency(matrix.cor2, weight = T, mode = "undirected")
  g2 <- simplify(g2)
  V(g2)$label <- V(g2)$name
  V(g2)$degree <- degree(g2)

  g3 <- graph.adjacency(matrix.cor3, weight = T, mode = "undirected")
  g3 <- simplify(g3)
  V(g3)$label <- V(g3)$name
  V(g3)$degree <- degree(g3)

  # append the output into results
  result <- list()
  result$matrix.cor <- matrix.cor
  result$matrix.cor.p <- matrix.cor.p

  result$matrix.cor1 <- matrix.cor1
  result$graph1 <- g1

  result$matrix.cor2 <- matrix.cor2
  result$graph2 <- g2

  result$matrix.cor3 <- matrix.cor3
  result$graph3 <- g3
  return(result)
}



# Co-occurrence-network-analysis
## OTU filtering, network generation, topological analysis and export OTU table
library(igraph)
library(Hmisc)

setwd("I:\\MAJORBIO\\out-work\\Material\\data")


Abu <- read.table("oturelative.txt", header = T)
Abu <- read.table("otuabu.txt", header = T)
Abu <- as.matrix(Abu)

### Filtering OTUs
table <- Abu
table[table > 0] <- 1
table.generalist <- Abu[which(rowSums(table) >= 12), ]
Abu <- table.generalist

## Creating gml files of network (to be visulized in Gephi or Cytoscape)

## cutoffs for correlation coefficient and P-value
pattern <- coRnetwork(Abu, 0.6, 0.01)

write.graph(pattern$graph1, "Pos0.6-rela.gml", format = "gml") # network file for positive association
write.graph(pattern$graph2, "Neg0.6-rela.gml", format = "gml") # network file for negative association
write.graph(pattern$graph3, "PosNeg0.6-rela.gml", format = "gml") # network file for all association

write.graph(pattern$graph1, "Pos0.6-abu.gml", format = "gml") # network file for positive association
write.graph(pattern$graph2, "Neg0.6-abu.gml", format = "gml") # network file for negative association
write.graph(pattern$graph3, "PosNeg0.6-abu.gml", format = "gml") # network file for all association


## Calculating network topological properties and viz plot

g <- pattern$graph1 ### positive network
# g<-pattern$graph1   ###negative network
g <- pattern$graph3 ### positive&&negative network
plot(g)

c <- cluster_walktrap(g)

## Global toplogical features
modularity(c)
md <- modularity(g, membership(c), weights = NULL)
cc <- transitivity(g,
  vids = NULL,
  weights = NULL
)
spl <- average.path.length(g, directed = FALSE, unconnected = TRUE)
gd <- graph.density(g, loops = FALSE)
nd <- diameter(g, directed = FALSE, unconnected = TRUE, weights = NA)

node.degree <- degree(g, v = V(g), mode = "all")
ad <- mean(node.degree)

e <- ecount(g)
v <- vcount(g)
global.topology <- data.frame(e, v, cc, spl, md, gd, nd, ad)
write.csv(global.topology, file = "Pos0.6-rela-global.topology.csv")
write.csv(global.topology, file = "Pos0.6-abu-global.topology.csv")
# Node toplogical features
betweenness.centrality <- betweenness(g,
  v = V(g),
  directed = FALSE, weights = NA,
  nobigint = TRUE, normalized = FALSE
)
closeness.centrality <- closeness(g,
  vids = V(g),
  weights = NA, normalized = FALSE
)
node.transitivity <- transitivity(g,
  type = c("local"), vids = NULL,
  weights = NA
)

node.topology <- data.frame(node.degree, betweenness.centrality, closeness.centrality, node.transitivity)
write.csv(node.topology, file = "Pos0.6-rela-node.topology.csv")
write.csv(node.topology, file = "Pos0.6-abu-node.topology.csv")
# Ploting node degreee distribution in a log-log plot
degree.df <- data.frame(table(degree = factor(node.degree, levels = seq_len(max(node.degree)))))
degree.df$degree <- as.numeric(as.character(degree.df$degree))

### Creating an abundance table for OTUs present in the positive and negative network
my.list1 <- row.names(pattern$matrix.cor1)
### my.list2 <- row.names(pattern$matrix.cor2)

logical1 <- row.names(Abu) %in% my.list1
### logical2 <- row.names(Abu)  %in% my.list2

tab.subset1 <- subset(Abu, logical1)
### tab.subset2 <- subset(Abu,logical2)

write.table(tab.subset1, "Pos0.6-Abu.txt", sep = "\t")
### write.table(tab.subset2,'Neg0.6-NW.txt',sep="\t")
