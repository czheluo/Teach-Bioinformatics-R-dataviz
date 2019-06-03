#############################################
#  http://www.majorbio.com/
#  Copyright (C) 2019 in  Majorbio workshop
#  contact: meng.luo@majorbio.com
#############################################

# require packages

library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)
color <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = T)]
rcolor <- color[sample(1:length(color), length(color))]

# example
# setosa (cimangyegucao) versicolor (yuanwei) (yanjiaocao)
head(iris)
df <- iris
df <- as.data.frame(iris)
row.names(df) <- paste(df$Species, row.names(df), sep = "_")
# df$Species <- NULL
head(df)

df.pca <- prcomp(df)
plot(df.pca$x[, 1], df.pca$x[, 2])
df.out <- as.data.frame(df.pca$x)
df.out$group <- sapply(strsplit(as.character(row.names(df)), "_"), "[[", 1)
head(df_out)

# plot

p <- ggplot(df.out, aes(x = PC1, y = PC2, color = group))
p <- p + geom_point()
p

theme <- theme(
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

p <- ggplot(df.out, aes(x = PC1, y = PC2, color = group))
p <- p + geom_point(size = 6) + theme
p
# add label text
p <- ggplot(df.out, aes(x = PC1, y = PC2, color = group, label = row.names(df)))
p <- p + geom_point() + geom_text(size = 3) + theme
p

percentage <- round(df.pca$sdev / sum(df.pca$sdev) * 100, 2)
percentage <- paste(
  colnames(df.out),
  "(", paste(as.character(percentage), "%", ")", sep = "")
)

p <- ggplot(df.out, aes(x = PC1, y = PC2, color = group))
p <- p + geom_point() + theme + xlab(percentage[1]) + ylab(percentage[2])
p


# change color

df.out$group <- factor(df.out$group, levels = c("virginica", "setosa", "versicolor"))

p <- ggplot(df.out, aes(x = PC1, y = PC2, color = group))
p <- p + geom_point() + theme + xlab(percentage[1]) + ylab(percentage[2]) +
  scale_color_manual(values = c("#FFFF00", "#00FFFF", "#FF00FF"))
p

# save image
pdf("pca.pdf", width = 10, height = 10)

yy <- grid.arrange(p, nrow = 1)
op <- par(no.readonly = TRUE)
par(op)

dev.off()
# or
pdf("pca.p.pdf", width = 10, height = 10)
p
dev.off()

# real dataset

pca <- read.csv("human.matrix.csv", header = T)
head(pca)

gro <- read.csv("group.csv", header = F)
head(gro)
human.pca <- prcomp(as.matrix(t(pca[, c(2:27)])), scale = T)
human.pca.out <- as.data.frame(human.pca$x)

human.pca.out$group <- gro[, 1]

head(human.pca.out)

p <- ggplot(human.pca.out, aes(x = PC1, y = PC2, color = group))
p <- p + geom_point()
p

p <- ggplot(human.pca.out, aes(x = PC1, y = PC2, color = group))
p <- p + geom_point() + theme
p
# add label text
p <- ggplot(human.pca.out, aes(x = PC1, y = PC2, color = group, label = row.names(human.pca.out)))
p <- p + geom_point() + geom_text(size = 3) + theme
p

percentage <- round(human.pca$sdev / sum(human.pca$sdev) * 100, 2)
percentage <- paste(
  colnames(human.pca.out),
  "(", paste(as.character(percentage), "%", ")", sep = "")
)

p <- ggplot(human.pca.out, aes(x = PC1, y = PC2, color = group))
p <- p + geom_point() + theme + xlab(percentage[1]) + ylab(percentage[2])
p


# change color

human.pca.out$group <- factor(human.pca.out$group,
  levels = c(
    "V_24h", "V_6h", "A_24h", "AV_24h",
    "C_6h", "MOCK", "A_6h", "C_24h", "AV_6h"
  )
)
png(paste("pca",".png",sep=""),width=1000, height=800)
p <- ggplot(human.pca.out, aes(x = PC1, y = PC2, color = group))
p <- p + geom_point(size=6)+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=16),
        #legend.box.background=element_rect("white"),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        plot.margin = unit(c(1, 1, 1, 1), "line")) + 
  xlab(percentage[1]) + ylab(percentage[2]) +
  scale_color_manual(values = rcolor[sample(1:100)[c(1:9)]])
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
