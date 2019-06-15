##########################################
#  http://www.majorbio.com/
#  Copyright (C) 2019 in  Majorbio workshop
#  contact: meng.luo@majorbio.com
##########################################
## wb:https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html


#require
library(tidyverse)
library(DESeq2)
library(RColorBrewer)
color <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = T)]
rcolor <- color[sample(1:length(color), length(color))]

myWO<- read.csv("WT_KO_count.csv",header = T,row.names = "ensgene")
metWO<- read.csv("WT_KO.csv",header = T,row.names = "name")

all(rownames(metWO) %in% colnames(myWO))

all(rownames(metWO) == colnames(myWO))

dds <- DESeqDataSetFromMatrix(countData = myWO,
                              colData = metWO,
                              design = ~ dex)
dds

dds <- DESeq(dds)

sizeFactors(dds)
dispersions(dds)
results(dds)

res <- results(dds, tidy=TRUE)
res <- tbl_df(res)
res

res %>% 
  filter(padj<0.05) %>% 
  write_csv("sigresults.csv")

plotCounts(dds, gene="ENSMUSG00000002489", intgroup="dex")

# Return the data
plotCounts(dds, gene="ENSMUSG00000002489", intgroup="dex", returnData=TRUE)

# Plot it
plotCounts(dds, gene="ENSMUSG00000002489", intgroup="dex", returnData=TRUE) %>% 
  ggplot(aes(dex, count)) + 
  geom_boxplot(aes(fill=dex)) + scale_y_log10() 

# SET THEME
theme <- theme(
  axis.text = element_text(size = 16),
  axis.title = element_text(size = 16, face = "bold"),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  #panel.background = element_blank(),
  #panel.border = element_rect(fill = NA),
  #panel.grid.major = element_blank(),
  #panel.grid.minor = element_blank(),
  #strip.background = element_blank(),
  axis.text.x = element_text(colour = "black"),
  axis.text.y = element_text(colour = "black"),
  axis.ticks = element_line(colour = "black"),
  plot.margin = unit(c(1, 1, 1, 1), "line")
)

# Create the new column
res <- res %>% mutate(sig=padj<0.05)

# How many of each?
res %>% 
  group_by(sig) %>% 
  summarize(n=n())

res %>% 
  filter(!is.na(log2FoldChange)) %>% 
  ggplot(aes(baseMean, log2FoldChange, col=sig)) + 
  geom_point() + 
  scale_fill_manual(values = rcolor[c(1:3)])+
  scale_color_manual(values = rcolor[c(1:3)])+
  scale_x_log10() + 
  ggtitle("MA plot")

res %>% 
  filter(!is.na(log2FoldChange) & !is.na(pvalue)) %>% 
  ggplot(aes(log2FoldChange, -1*log10(pvalue), col=sig)) + 
  geom_point() + theme+
  scale_color_manual(values = rcolor[sample(1:100)[c(1:3)]])+
  ggtitle("Volcano plot")

### LOAD result

expres <- read.csv("WT_vs_KO.deseq2.xls",sep="\t",header = T)
expres <- expres %>% mutate(sig=padj<0.05)

head(expres)

png(paste("maplot", ".png", sep = ""), width = 1000, height = 800)

expres %>% 
  filter(!is.na(log2fc)) %>% 
  ggplot(aes(-1*log10(pvalue), log2fc, col=significant)) + 
  geom_point() + 
  scale_color_manual(values = rcolor[sample(1:100)[c(1:2)]])+
  #scale_x_log10() + 
  geom_hline(yintercept = c(-1,1), linetype = 4)+
  ggtitle("MA plot") +theme

dev.off()

pdf("Volcano.pdf", width = 10, height = 10)
expres %>% 
  filter(!is.na(log2fc) & !is.na(padjust)) %>% 
  ggplot(aes(log2fc, -1*log10(padjust), col=significant)) + 
  geom_point() + 
  ggtitle("Volcano plot") +theme

dev.off()












