##########################################
#  http://www.majorbio.com/
#  Copyright (C) 2019 in  Majorbio workshop
#  contact: meng.luo@majorbio.com
##########################################


set.seed(0614)

# require packages
library(ggplot2)
library(scales)
library(gridExtra)
library(RColorBrewer)
color <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = T)]
rcolor <- color[sample(1:length(color), length(color))]
library(EnhancedVolcano)
# setting data director

expres <- read.table("VC_24h_vs_AZ_24h.deseq2.xls", sep = "\t", header = T)

cols <- c("yes" = "red", "no" = "orange", "nosig" = "darkgrey", "up" = "#00B2FF", "down" = "#00B2FF")
# scatter plot
vol <- ggplot(expres, aes(
  x = log10(VC_24h_tpm),
  y = log10(AZ_24h_tpm),
  color = regulate
))
vol + geom_point(size = 2.5) +
  # xlim(c(0,5))+ylim(c(0,5))
  scale_colour_manual(values = c("#00B2FF", "darkgrey", "orange")) +
  theme_bw(base_size = 14) + # change overall theme
  theme(legend.position = "right") + # change the legend
  xlab(expression("Log10(VC_24h)")) + # Change X-Axis label
  ylab(expression("log10(AZ_24h)")) # Change Y-Axis label

# VOLCANO
vol <- ggplot(expres, aes(x = log2fc, y = -log10(padjust)))

vol + geom_point()

vol + geom_point(color = "red")

vol + geom_point(aes(color = "red"))

vol + geom_point(aes(color = regulate))

volcano <- vol + geom_point(aes(color = regulate)) + xlim(-4, 4) + ylim(0, 310)

volcano + labs(title = "Volcanoplot", x = "log2(FC)")

volcano + labs(
  title = "Volcanoplot", x = expression(log[2](FC)),
  y = expression(-log[10](FDR))
)

volcanop <- volcano + labs(
  title = "Volcanoplot", x = expression(log[2](FC)),
  y = expression(-log[10](FDR))
)

volcanop + scale_color_manual(values = c("green", "black", "red"))

volcanor <- volcanop + scale_color_manual(values = c("#00ba38", "#619cff", "#f8766d"))
volcanor + geom_hline(yintercept = 13) + geom_vline(xintercept = c(-1, 1))

volcanor + geom_hline(yintercept = 12, linetype = 4) + geom_vline(xintercept = c(-1, 1), linetype = 4)

ggsave("volcano.png")
ggsave("volcano8.png", volcano, width = 8, height = 8)

## R packages

png(paste("huamn.vol", ".png", sep = ""), width = 1000, height = 800)
EnhancedVolcano(expres,
  lab = rownames(expres),
  x = "log2fc",
  y = "padjust",
  # xlim = c(-6, 6),
  selectLab = "",
  title = "",
  pCutoff = 10e-12,
  FCcutoff = 1,
  xlab = "Log2FC",
  transcriptPointSize = 1.5,
  transcriptLabSize = 3.0,
  colAlpha = 1,
  cutoffLineType = "blank",
  cutoffLineCol = "black",
  cutoffLineWidth = 1,
  legendLabSize = 14,
  legendIconSize = 4,
  # transcriptPointSize =2, #c(1.6,1.6,1.6,1.6),
  # transcriptLabSize = 3, #c(4,4,4,4),
  # hline = c(10e-12, 10e-36, 10e-60, 10e-84),
  # hlineCol = c('grey0', 'grey25','grey50','grey75'),
  # hlineType = 'longdash',
  # hlineWidth = 0.8,
  gridlines.major = FALSE,
  gridlines.minor = FALSE
)
dev.off()
pdf(paste("human.vol", ".pdf", sep = ""), height = 9, width = 16)
EnhancedVolcano(expres,
  lab = rownames(expres),
  x = "log2fc",
  y = "padjust",
  # xlim = c(-6, 6),
  selectLab = "",
  title = "",
  pCutoff = 10e-12,
  FCcutoff = 1,
  xlab = "Log2FC",
  transcriptPointSize = 1.5,
  transcriptLabSize = 3.0,
  colAlpha = 1,
  cutoffLineType = "blank",
  cutoffLineCol = "black",
  cutoffLineWidth = 1,
  legendLabSize = 14,
  legendIconSize = 4,
  # hline = c(10e-12, 10e-36, 10e-60, 10e-84),
  # hlineCol = c('grey0', 'grey25','grey50','grey75'),
  # hlineType = 'longdash',
  # hlineWidth = 0.8,
  gridlines.major = FALSE,
  gridlines.minor = FALSE
)
dev.off()
