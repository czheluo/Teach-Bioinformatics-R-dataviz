##########################################
#  http://www.majorbio.com/
#  Copyright (C) 2019 in Majorbio workshop
#  contact: meng.luo@majorbio.com
##########################################

# requrie packages
library(ComplexHeatmap)
library(nVennR)
library(RColorBrewer)
color <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = T)]
rcolor <- color[sample(1:length(color), length(color))]

df <- data.frame(
  A = sample(c(0, 1), 100, replace = T),
  B = sample(c(0, 1), 100, replace = T),
  C = sample(c(0, 1), 100, replace = T),
  D = sample(c(0, 1), 100, replace = T),
  E = sample(c(0, 1), 100, replace = T)
)

fromBin <- function(binList) {
  result <- 0
  for (b in binList) {
    result <- bitwShiftL(result, 1)
    result <- result + b
  }
  return(result)
}
interpretBin <- function(dff) {
  nels <- bitwShiftL(1, ncol(dff))
  result <- vector(mode = "numeric", length = nels)
  for (r in rownames(dff)) {
    n <- fromBin(dff[r, ]) + 1
    result[n] <- result[n] + 1
  }
  return(result)
}

regs <- interpretBin(df)
myV <- createVennObj(nSets = ncol(df), sNames = colnames(df), sSizes = regs)
myV <- plotVenn(nVennObj = myV, outFile = "mnV.svg")
myV <- plotVenn(nVennObj = myV, outFile = "mnV.svg")

myV4 <- plotVenn(list(a = c(1, 2, 3), b = c(3, 4, 5), c = c(3, 6, 1)), nCycles = 2000, setColors = c("red", "green", "blue"), labelRegions = F, fontScale = 2, opacity = 0.2, borderWidth = 2)

showSVG(myV4, outFile = "mnV4.svg")


library(nVennR)

setwd("")
v1 <- read.csv("1.csv", header = T)
v2 <- read.csv("2.csv", header = T)
v3 <- read.csv("3.csv", header = T)
v4 <- read.csv("4.csv", header = T)

myV4 <- plotVenn(list(PTLD4 = v1$circbase_ID, PTLD5 = v2$circbase_ID, PTLD6 = v3$circbase_ID, PTLD9 = v4$circbase_ID),
  nCycles = 2000, setColors = c("red", "green", "blue", "yellow"),
  labelRegions = F, fontScale = 2, opacity = 0.2, borderWidth = 2, outFile = "mnVR.svg"
)

showSVG(myV4, opacity = 0.8, systemShow = T)




library(VennDiagram)
venn.diagram(
  x = list(
    PTLD4 = v1$circbase_ID, PTLD5 = v2$circbase_ID,
    PTLD6 = v3$circbase_ID, PTLD9 = v4$circbase_ID
  ),
  category.names = c("PTLD4", "PTLD5", "PTLD6", "PTLD9"),
  filename = "4_venn_diagramm.png",
  output = TRUE,
  imagetype = "png",
  height = 1000,
  width = 800,
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  lty = "blank",
  fill = c("yellow", "purple", "green", "blue"),
  cex = 0.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer"
  # cat.pos = c(-27, 27, 135),
  # cat.dist = c(0.055, 0.055, 0.085),
  # cat.fontfamily = "sans",
  # rotation = 1
)

#
library(ComplexHeatmap)
png(paste("upset",".png",sep=""),width=1000, height=800)
lt <- list(
  PTLD4 = v1$circbase_ID, PTLD5 = v2$circbase_ID,
  PTLD6 = v3$circbase_ID, PTLD9 = v4$circbase_ID
)
m <- make_comb_mat(lt)
UpSet(m,
      queries = list(list(query=intersects, 
                          params=list("PTLD6", "PTLD9"), color="yellow", active=T)))


dev.off()

setwd("G:\\MAJORBIO\\R Course\\Intro\\data\\venn")
v1 <- data.frame(read.table(list.files()[1], header = F)[, 1])
na <- strsplit(list.files()[1], split = "_")
colnames(v1) <- paste(na[[1]][c(1)], na[[1]][c(2)], na[[1]][c(3)], na[[1]][c(4)], na[[1]][c(2)], sep = "_")
v2 <- data.frame(read.table(list.files()[2], header = F)[, 1])
na <- strsplit(list.files()[2], split = "_")
rownames(v2) <- paste(na[[1]][c(1)], na[[1]][c(2)], na[[1]][c(3)], na[[1]][c(4)], na[[1]][c(2)], sep = "_")
v3 <- data.frame(read.table(list.files()[3], header = F)[, 1])
na <- strsplit(list.files()[3], split = "_")
colnames(v3) <- paste(na[[1]][c(1)], na[[1]][c(2)], na[[1]][c(3)], na[[1]][c(4)], na[[1]][c(2)], sep = "_")
v4 <- data.frame(read.table(list.files()[4], header = F)[, 1])
na <- strsplit(list.files()[4], split = "_")
colnames(v4) <- paste(na[[1]][c(1)], na[[1]][c(2)], na[[1]][c(3)], na[[1]][c(4)], na[[1]][c(2)], sep = "_")
v5 <- data.frame(read.table(list.files()[5], header = F)[, 1])
na <- strsplit(list.files()[5], split = ".")
colnames(v5) <- paste(na[[1]][c(1)], na[[1]][c(2)], na[[1]][c(3)], na[[1]][c(4)], na[[1]][c(2)], sep = "_")
v6 <- data.frame(read.table(list.files()[6], header = F)[, 1])
na <- strsplit(list.files()[6], split = "_")
colnames(v6) <- paste(na[[1]][c(1)], na[[1]][c(2)], na[[1]][c(3)], na[[1]][c(4)], na[[1]][c(2)], sep = "_")
v7 <- data.frame(read.table(list.files()[7], header = F)[, 1])
na <- strsplit(list.files()[7], split = "_")
colnames(v7) <- paste(na[[1]][c(1)], na[[1]][c(2)], na[[1]][c(3)], na[[1]][c(4)], na[[1]][c(2)], sep = "_")
v8 <- data.frame(read.table(list.files()[8], header = F)[, 1])
na <- strsplit(list.files()[8], split = "_")
colnames(v8) <- paste(na[[1]][c(1)], na[[1]][c(2)], na[[1]][c(3)], na[[1]][c(4)], na[[1]][c(2)], sep = "_")




lt <- list(
  v1 = v1, v3 = v3, v5 = v5, v6 = v6, v7 = v7, v8 = v8
)
# names(lt)<- c(colnames(v1), colnames(v3), colnames(v5), colnames(v6),
#              colnames(v7), colnames(v8))
names(lt) <- c(
  "C_24h_vs_DC_24h", "C_6h_vs_DC_6h", "MOCK_vs_C_24h", "MOCK_vs_C_6h",
  "VC_24h_vs_AZ_24h", "VC_6h_vs_AZ_6h"
)
m <- make_comb_mat(lt)
UpSet(m)


lt <- list(
  v1 = v1, v3 = v3, v5 = v5, v7 = v7, v8 = v8
)

names(lt) <- c(
  "C_24h_vs_DC_24h", "C_6h_vs_DC_6h", "MOCK_vs_C_24h",
  "VC_24h_vs_AZ_24h", "VC_6h_vs_AZ_6h"
)
library(VennDiagram)
venn.diagram(
  x = lt,
  category.names = names(lt),
  filename = "6_venn_diagramm.png",
  output = TRUE,
  imagetype = "png",
  height = 1000,
  width = 800,
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  lty = "blank",
  fill = c("yellow", "purple", "green", "blue", "red"),
  cex = 0.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer", total.population = 6,

  # cat.fontfamily = "sans",
  # rotation = 1
)




lt <- list(
  v1 = v1, v3 = v3, v5 = v5, v6 = v6, v7 = v7, v8 = v8
)
# names(lt)<- c(colnames(v1), colnames(v3), colnames(v5), colnames(v6),
#              colnames(v7), colnames(v8))
names(lt) <- c(
  "C_24h_vs_DC_24h", "C_6h_vs_DC_6h", "MOCK_vs_C_24h", "MOCK_vs_C_6h",
  "VC_24h_vs_AZ_24h", "VC_6h_vs_AZ_6h"
)

myV4 <- plotVenn(lt,
  nCycles = 2000,
  setColors = c("yellow", "purple", "green", "blue", "red", "black"),
  labelRegions = F, fontScale = 2, opacity = 0.2,
  borderWidth = 2, outFile = "mnVR.svg"
)

showSVG(myV4, opacity = 0.8, systemShow = T)
