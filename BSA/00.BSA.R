##########################################
#  http://www.majorbio.com/
#  Copyright (C) 2020 in Majorbio workshop
#  contact: meng.luo@majorbio.com
##########################################

################## required packages ###############################

#devtools::install_github("bmansfeld/QTLseqr")
#install.packages("vcfR")
#install.packages("tidyr")
#install.packages("dplyr")
#install.packages("CMplot")

library(CMplot)
library(QTLseqr)
library(vcfR)
library(tidyr)
library(dplyr)

################## set workdir ###############################

#setwd("D:\\BSA\\BSA")

vcf <- read.vcfR("pop.no.vcf.gz")

ltable <- data.frame(CHROM = getCHROM(vcf),
                     POS = getPOS(vcf),
                     REF = getREF(vcf),
                     ALT = getALT(vcf)
                     )

ad <- as.data.frame(extract.gt(vcf, "AD"))

ltable<-cbind(ltable,ad)

gt <- as.data.frame(extract.gt(vcf, "GT"))

#LowBulk
LB<-"YC-bulk"
LP<-"ZH30"

ltable<-ltable %>%
    separate(col = LB,
            into = c(paste("AD_REF.",LB,sep = ""),
                    paste("AD_ALT.",LB,sep = "")),
            sep = ",",convert = TRUE)

#HighBulk
HB<-"GC-bulk"
HP<-"JY102"
ltable<-ltable %>%
    separate(col = HB,
                    into = c(paste("AD_REF.",HB,sep = ""),
                             paste("AD_ALT.",HB,sep = "")),
                    sep = ",",convert = TRUE)#tidyr::



##filter markers
PS1 <- which(gt[,HP] == "0/0" &  gt[,LP] == "1/1")
#PS2 <- which(gt[,HP] == "1/1" &  gt[,LP] == "0/0")

ltable <- ltable[PS1,]


ltable <-ltable %>%
	mutate(REF_FRQ = (`AD_REF.GC-bulk` + `AD_REF.YC-bulk`) / (`AD_ALT.GC-bulk` + `AD_ALT.YC-bulk`
	+`AD_REF.GC-bulk` + `AD_REF.YC-bulk`))



table_filter<- filter(ltable, (`AD_REF.GC-bulk` +`AD_ALT.GC-bulk`+
                                          `AD_REF.YC-bulk`+`AD_ALT.YC-bulk`) >=
                                30)#dplyr::


#table_filter <- filter(table_filter, table_filter$REF_FRQ < 1 - refAlleleFreq &
#                table_filter$REF_FRQ > refAlleleFreq)

table_filter <- filter(table_filter, table_filter$REF_FRQ > 0.3)


load("tablefilter.RData")

write.table(table_filter, file = "soybeanF2.txt", sep = "\t", row.names = F, quote = F)

table <- importFromTable("soybeanF2.txt",
                      highBulk = HB,
                      lowBulk = LB,
                      chromList = unique(table_filter$CHROM),
                      sep = "\t")


table <- subset(table, !is.na(SNPindex.LOW) & !is.na(SNPindex.HIGH))



## G' stats

QTLseq_G <- runGprimeAnalysis(SNPset = table,
                        windowSize = 1e6,
                        outlierFilter = "deltaSNP")

### plot result ###
plotQTLStats(
    SNPset = QTLseq_G ,
    var = "Gprime",
    plotThreshold = TRUE,
    q = 0.001
)

getQTLTable(SNPset = QTLseq_G, alpha = 0.01, export = TRUE, fileName = "Gstat_QTL.csv")

### QTL-seq ###

QTLseq_G <- runQTLseqAnalysis(SNPset = QTLseq_G,
                        windowSize = 1e6,
                        popStruc = "F2",
                        bulkSize = c(30,30))


getQTLTable(SNPset = QTLseq_G, interval = 99, export = TRUE, fileName = "QTLseq_QTL.csv")

### plot result ###
plotQTLStats(
    SNPset = QTLseq_G,
    var = "deltaSNP",
    plotIntervals = TRUE)


##### Results VZ ######

library(CMplot)

SNP <- paste("SNP",c(1:length(QTLseq_G$CHROM)),sep="")
ch <- matrix(0,length(QTLseq_G$CHROM))
CHR <- unique(QTLseq_G$CHROM)

for (i in 1:length(CHR)){
    ch[which(QTLseq_G$CHROM %in% CHR[i])]<-rep(paste(i),length(QTLseq_G$CHROM==CHR))
}

Chr <- as.numeric(ch)

result <- data.frame(SNP,Chr,QTLseq_G$POS,QTLseq_G$Gprime,
                     QTLseq_G$tricubeDeltaSNP,QTLseq_G$CI_99)


getFDRThreshold <- function(pvalues, alpha = 0.01,method="BH"){
    sortedPvals <- sort(pvalues, decreasing = FALSE)
    pAdj <- p.adjust(sortedPvals, method = method)
    if (!any(pAdj < alpha)) {
        fdrThreshold <- NA
    } else {
        fdrThreshold <- sortedPvals[max(which(pAdj < alpha))]
    }
    return(fdrThreshold)
}


fdrT <- getFDRThreshold(QTLseq_G$pvalue, alpha = 0.01)

GprimeT <- QTLseq_G[which(QTLseq_G$pvalue == fdrT),"Gprime"]



CMplot(result[,c(1,2,3,4)], plot.type="m", LOG10=FALSE, cex.lab=2,ylab = NULL,
       ylim=c(round(min(result[,4]),digits=1),round(max(result[,4]),digits=1)),
       amplify=FALSE,main="G' value",chr.labels=NULL,
       threshold=GprimeT,threshold.lty=2,threshold.lwd=2, threshold.col="black",
       bin.size=1e6,type="p",multracks=FALSE,file.output=TRUE)#type="p"



CMplot(result[,c(1,2,3,4)], plot.type="m", LOG10=FALSE, cex.lab=2,ylab = NULL,
       ylim=c(round(min(result[,4]),digits=1),round(max(result[,4]),digits=1)),
       amplify=FALSE,main="G' value",chr.labels=NULL,
       threshold=GprimeT,threshold.lty=2,threshold.lwd=2, threshold.col="black",
       bin.size=1e6,type="p",multracks=FALSE,file.output=TRUE,
       chr.den.col=c("darkgreen", "yellow", "red"),
       )#type="p"


CMplot(result[,c(1,2,3,5)], plot.type="m", LOG10=FALSE, cex.lab=2,ylab = NULL,
       ylim=c(1,-1), amplify=FALSE,main="DeltaSNP",chr.labels=NULL,
       #threshold=c(min(result[,6]),-min(result[,6])),threshold.lty=c(2,2),
       #threshold.lwd=c(2,2), threshold.col=c("red","red"),
       bin.size=1e6,type="p",multracks=FALSE,file.output=TRUE,)#type="h"

CMplot(result[,c(1,2,3,5)], plot.type="m", LOG10=FALSE, cex.lab=2,ylab = NULL,
       ylim=c(1,-1), amplify=FALSE,main="DeltaSNP",chr.labels=NULL,
       threshold=c(min(result[,6]),-min(result[,6])),threshold.lty=c(2,2),
       threshold.lwd=c(2,2), threshold.col=c("red","red"),
       bin.size=1e6,type="p",multracks=FALSE,file.output=TRUE,
       chr.den.col=c("darkgreen", "yellow", "red"),
       signal.col=c("red","green"),signal.cex=c(1.5,1.5),
       signal.pch=c(19,19),file="jpg",
       )#type="h"


CMplot(result[,c(1:5)],type="p",plot.type="c",chr.labels=paste("Chr",c(1:20),sep=""),
       r=0.4,cir.legend=TRUE,outward=FALSE,cir.legend.col="black",LOG10=FALSE,
       cir.chr.h=1.3,chr.den.col="black",bin.size=1e6,ylab = NULL,file="jpg",
       memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10)


CMplot(result[,c(1:5)],type="p",plot.type="c",chr.labels=paste("Chr",c(1:20),sep=""),
       r=0.4,cir.legend=TRUE,outward=FALSE,cir.legend.col="black",LOG10=FALSE,
       cir.chr.h=1.3,ylab = NULL,file="jpg",#"jpg", "pdf", "tiff"
       memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10,
       chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6)




#png(paste("G",".png",sep=""),width=1000, height=800)
#dev.off()

#save.image("BSA.Rdata")

