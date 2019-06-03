##########################################
#  http://www.majorbio.com/
#  Copyright (C) 2019 in Majorbio workshop
#  contact: meng.luo@majorbio.com
##########################################

# https://github.com/jokergoo/circlize

# circlize
load("mlm.rds")
mlm <- is.na(mlm)
circos.par("track.height" = 0.2)
circos.initialize(factors = mlm$SNP, x = mlm$BP / 1000000)
circos.track(
  factors = mlm$SNP, y = mlm$effect,
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"),
      CELL_META$sector.index
    )
    circos.axis(labels.cex = 0.6)
  }
)

# color
install.packages('RColorBrewer')
library(RColorBrewer)
color <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = T)]
rcolor <- color[sample(1:length(color), length(color))]

circos.trackPoints(mlm$SNP, mlm$BP / 1000000, mlm$effect, col = rcolor[1:7], pch = 19, cex = 0.5)


## heatmap

circos.track(ylim = c(-7.998938, 3.869004), panel.fun = function(x, y) {
  xlim <- CELL_META$xlim
  ylim <- CELL_META$ylim
  breaks <- seq(xlim[1], xlim[2], by = 1)
  n_breaks <- length(breaks)
  circos.rect(breaks[-n_breaks], rep(ylim[1], n_breaks - 1),
    breaks[-1], rep(ylim[2], n_breaks - 1),
    col = rand_color(n_breaks), border = NA
  )
})

## link (three link type)
chr1 <- mlm[which(mlm$SNP %in% "chr1"), 3][1:7]
chr6 <- mlm[which(mlm$SNP %in% "chr3"), 3][1:7]
for (i in 1:7) {
  circos.link("chr1", chr1[i], "chr6", chr6[i], h = 0.4)
}
## p value

circos.clear()
circos.par("track.height" = 0.2)
circos.initialize(factors = mlm$SNP, x = mlm$BP / 1000000)
circos.track(
  factors = mlm$SNP, y = -log10(mlm$HORVU7Hr1G116180.MLM),
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"),
      CELL_META$sector.index
    )
    circos.axis(labels.cex = 0.6)
  }
)

circos.trackPoints(mlm$SNP, mlm$BP / 1000000, -log10(mlm$HORVU7Hr1G116180.MLM), col = rcolor[1:7], pch = 19, cex = 0.5)

## heatmap
circos.track(ylim = c(min(-log10(mlm$HORVU7Hr1G116180.MLM)), max(-log10(mlm$HORVU7Hr1G116180.MLM))), panel.fun = function(x, y) {
  xlim <- CELL_META$xlim
  ylim <- CELL_META$ylim
  breaks <- seq(xlim[1], xlim[2], by = 0.1)
  n_breaks <- length(breaks)
  circos.rect(breaks[-n_breaks], rep(ylim[1], n_breaks - 1),
    breaks[-1], rep(ylim[2], n_breaks - 1),
    col = rand_color(n_breaks), border = NA
  )
})
## link (three link type)
chr1 <- mlm[which(mlm$SNP %in% "chr1"), 3][1:20]
chr6 <- mlm[which(mlm$SNP %in% "chr6"), 3][1:20]
chr2 <- mlm[which(mlm$SNP %in% "chr2"), 3][1:20]
chr7 <- mlm[which(mlm$SNP %in% "chr7"), 3][1:20]
chr4 <- mlm[which(mlm$SNP %in% "chr4"), 3][1:20]
chr5 <- mlm[which(mlm$SNP %in% "chr5"), 3][1:20]
chr3 <- mlm[which(mlm$SNP %in% "chr3"), 3][1:20]
for (i in 1:7) {
  circos.link("chr1", chr1[i], "chr6", chr6[i],
    col = rcolor[which(unique(mlm$SNP) %in% "chr6")],
    h = 0.4
  )
  circos.link("chr2", chr2[i], "chr7", chr7[i],
    col = rcolor[which(unique(mlm$SNP) %in% "chr7")],
    h = 0.5
  )
  circos.link("chr4", chr4[i], "chr5", chr5[i],
    col = rcolor[which(unique(mlm$SNP) %in% "chr5")],
    h = 0.5
  )
  circos.link("chr3", chr3[i], "chr4", chr4[i],
    col = rcolor[which(unique(mlm$SNP) %in% "chr4")],
    h = 0.5
  )
  circos.link("chr2", chr2[i], "chr5", chr4[i],
    col = rcolor[which(unique(mlm$SNP) %in% "chr4")],
    h = 0.5
  )
}
# save.image("circoslize.RData")
