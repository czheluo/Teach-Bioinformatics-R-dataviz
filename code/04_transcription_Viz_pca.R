#############################################
#  http://www.majorbio.com/
#  Copyright (C) 2019 in  Majorbio workshop
#  contact: meng.luo@majorbio.com
#############################################

set.seed(0614)

# require packages
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)
color <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = T)]
rcolor <- color[sample(1:length(color), length(color))]

# real dataset

pca <- read.csv("human.matrix.csv", header = T)
head(pca)

gro <- read.csv("human.group.csv", header = F)
head(gro)
# Principle component analysis

human.pca <- prcomp(as.matrix(t(pca[, c(2:27)])), scale = T)
human.pca.out <- as.data.frame(human.pca$x)
human.pca.out$group <- gro[, 1]

head(human.pca.out)
# SET THEM
theme <- theme(
  axis.text = element_text(size = 16),
  axis.title = element_text(size = 16, face = "bold"),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  panel.background = element_blank(),
  panel.border = element_rect(fill = NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  axis.text.x = element_text(colour = "black"),
  axis.text.y = element_text(colour = "black"),
  axis.ticks = element_line(colour = "black"),
  plot.margin = unit(c(1, 1, 1, 1), "line")
)

p <- ggplot(human.pca.out, aes(x = PC1, y = PC2, color = group))
p <- p + geom_point(size = 6) + theme
p

# add label text
p <- ggplot(human.pca.out, aes(x = PC1, y = PC2, color = group, label = row.names(human.pca.out)))
p <- p + geom_point(size = 4) + 
  geom_text(size = 3,position = position_nudge(y=2,x=-2)) + theme ## or  use vjust and hjust 
p
# explation

percentage <- round(human.pca$sdev / sum(human.pca$sdev) * 100, 2)
percentage <- paste(
  colnames(human.pca.out),
  "(", paste(as.character(percentage), "%", ")", sep = "")
)

p <- ggplot(human.pca.out, aes(x = PC1, y = PC2, color = group))
p <- p + geom_point(size = 6) + theme + xlab(percentage[1]) + ylab(percentage[2])
p

# change color

human.pca.out$group <- factor(human.pca.out$group,
  levels = c(
    "V_24h", "V_6h", "A_24h", "AV_24h",
    "C_6h", "MOCK", "A_6h", "C_24h", "AV_6h"
  )
)
png(paste("pca", ".png", sep = ""), width = 1000, height = 800)
p <- ggplot(human.pca.out, aes(x = PC1, y = PC2, color = group))
p <- p + geom_point(size = 6) +
  xlab(percentage[1]) + ylab(percentage[2]) +
  scale_color_manual(values = rcolor[sample(1:100)[c(1:9)]]) + theme
p
dev.off()

# save image

pdf("human.pca.pdf", width = 10, height = 10)

yy <- grid.arrange(p, nrow = 1)
op <- par(no.readonly = TRUE)
par(op)

dev.off()
# or

pdf("huaman.pca.p.pdf", width = 10, height = 10)
p
dev.off()




# CLEAN DATA
rm(list = ls(all.names = T))

# clean graphy

graphics.off()
