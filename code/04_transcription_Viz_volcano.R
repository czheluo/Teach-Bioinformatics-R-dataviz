##########################################
#  http://www.majorbio.com/
#  Copyright (C) 2019 in  Majorbio workshop
#  contact: meng.luo@majorbio.com
##########################################

# require packages
library(ggplot2)
library(scales)
library(gridExtra)
library(RColorBrewer)
color <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = T)]
rcolor <- color[sample(1:length(color), length(color))]

# setting data director

setwd("G:\\MAJORBIO\\R Course\\Intro\\data")

expres <- read.table("VC_24h_vs_AZ_24h.deseq2.xls",sep="\t",header = T)
# load("expres.rds")
cols <- c("yes" = "red", "no" = "orange", "nonsignificant" = "darkgrey", "up" = "#00B2FF", "down" = "#00B2FF")

vol <- ggplot(expres, aes(x = log2fc, y = -log10(padjust),fill = significant, labels=seq_id))

p <- vol + ggtitle(label = "Volcano Plot") +
  #scale_colour_manual(values = cols) +
  geom_point(size = 2.5, alpha = 1, na.rm = T, shape = 21, colour = "black") +
  theme_bw(base_size = 14) + # change overall theme
  theme(legend.position = "right") + # change the legend
  xlab(expression("Log2FC")) + # Change X-Axis label
  ylab(expression(-log[10]("FDR"))) + # Change Y-Axis label
  #geom_hline(yintercept = 1, colour="#990000", linetype="dashed") +
  geom_vline(xintercept = 2, colour="#990000", linetype="dashed") +
  geom_vline(xintercept = -2, colour="#990000", linetype="dashed") #+ # Add p-adj value cutoff line
  #scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values

#vol+scale_fill_npg()
#p+scale_color_aaas()
  
  expres$v24<-log10(rowMeans(expres[,c(8:10)]))
  expres$av24<-log10(rowMeans(expres[,c(11:13)]))#significant
  vol <- ggplot(expres, aes(x =log10(VC_24h_tpm), y = log10(AZ_24h_tpm),fill = regulate, labels=seq_id))
  
  vol + #ggtitle(label = "Volcano Plot", subtitle = "Colored by fold-change direction") +
    scale_colour_manual(values = cols) +
    geom_point(size = 2.5, alpha = 1, na.rm = T,shape = 21, colour = "white") +
    theme_bw(base_size = 14) + # change overall theme
    theme(legend.position = "right") + # change the legend
    xlab(expression("Log10(VC_24h)")) + # Change X-Axis label
    ylab(expression("log10(AZ_24h)")) + # Change Y-Axis label
    xlim(c(0,5))+ylim(c(0,5))+
    #geom_hline(yintercept = 1, colour="#990000", linetype="dashed") +
    #geom_vline(xintercept = 1, colour="#990000", linetype="dashed") +
    #geom_vline(xintercept = -1, colour="#990000", linetype="dashed") #+ # Add p-adj value cutoff line
  scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values
  


if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')
BiocManager::install('airway')
library(EnhancedVolcano)
library(airway)
library(magrittr)

data('airway')
airway$dex %<>% relevel('untrt')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library('DESeq2')

dds <- DESeqDataSet(airway, design = ~ cell + dex)
dds <- DESeq(dds, betaPrior=FALSE)
res1 <- results(dds,
                contrast = c('dex','trt','untrt'))
res1 <- lfcShrink(dds,
                  contrast = c('dex','trt','untrt'), res=res1)
res2 <- results(dds,
                contrast = c('cell', 'N061011', 'N61311'))
res2 <- lfcShrink(dds,
                  contrast = c('cell', 'N061011', 'N61311'), res=res2)
EnhancedVolcano(res1,
                lab = rownames(res1),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-5, 8))
head(res1)
EnhancedVolcano(expres,lab = rownames(expres),
                x= 'log2fc',
                y= 'padjust',
)
png(paste("huamn.vol",".png",sep=""),width=1000, height=800)
EnhancedVolcano(expres,lab = rownames(expres),
                x= 'log2fc',
                y= 'padjust',
               # xlim = c(-6, 6),
                selectLab="",
                title = "",
                pCutoff = 10e-12,
                FCcutoff = 1,
                xlab = "Log2FC",
                transcriptPointSize = 1.5,
                transcriptLabSize = 3.0,
                colAlpha = 1,
                cutoffLineType = 'blank',
                cutoffLineCol = 'black',
                cutoffLineWidth =1,
                legendLabSize = 14, 
                legendIconSize = 4,
                transcriptPointSize =2, #c(1.6,1.6,1.6,1.6), 
                transcriptLabSize = 3, #c(4,4,4,4), 
                #hline = c(10e-12, 10e-36, 10e-60, 10e-84),
                #hlineCol = c('grey0', 'grey25','grey50','grey75'),
                #hlineType = 'longdash',
                #hlineWidth = 0.8,
                gridlines.major = FALSE,
                gridlines.minor = FALSE)
dev.off()
pdf(paste("human.pc",".pdf",sep=""),height=9,width=16)
EnhancedVolcano(expres,lab = rownames(expres),
                x= 'log2fc',
                y= 'padjust',
                # xlim = c(-6, 6),
                selectLab="",
                title = "",
                pCutoff = 10e-12,
                FCcutoff = 1,
                xlab = "Log2FC",
                transcriptPointSize = 1.5,
                transcriptLabSize = 3.0,
                colAlpha = 1,
                cutoffLineType = 'blank',
                cutoffLineCol = 'black',
                cutoffLineWidth =1,
                legendLabSize = 14, 
                legendIconSize = 4,
                #hline = c(10e-12, 10e-36, 10e-60, 10e-84),
                #hlineCol = c('grey0', 'grey25','grey50','grey75'),
                #hlineType = 'longdash',
                #hlineWidth = 0.8,
                gridlines.major = FALSE,
                gridlines.minor = FALSE)
dev.off()
