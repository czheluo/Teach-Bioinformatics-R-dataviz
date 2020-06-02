##########################################
#  http://www.majorbio.com/
#  Copyright (C) 2019 in Majorbio workshop
#  contact: meng.luo@majorbio.com
##########################################

################## required packages ###############################


#devtools::install_github("bmansfeld/QTLseqr")
#install.packages("vcfR")
#install.packages("tidyr")

library(QTLseqr)
library(vcfR)
library(tidyr)
library(dplyr)
################## set workdir ###############################

setwd("D:\\R\\BSA")
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
    tidyr::separate(col = LB,
                    into = c(paste("AD_REF.",LB,sep = ""),
                             paste("AD_ALT.",LB,sep = "")),
                    sep = ",",convert = TRUE)

#HighBulk
HB<-"GC-bulk"
HP<-"JY102"
ltable<-ltable %>% 
    tidyr::separate(col = HB,
                    into = c(paste("AD_REF.",HB,sep = ""),
                             paste("AD_ALT.",HB,sep = "")),
                    sep = ",",convert = TRUE)


##filter markers

table_filter<- dplyr::filter(ltable, (`AD_REF.GC-bulk` +`AD_ALT.GC-bulk`+
                                          `AD_REF.YC-bulk`+`AD_ALT.YC-bulk`) >= 
                                 10)

PS1 <- which(gt[,HP] == "0/0" &  gt[,LP] == "1/1")
#PS2 <- which(gt[,HP] == "1/1" &  gt[,LP] == "0/0")

table <- ltable[PS1,]



write.table(table, file = "soybean.txt", sep = "\t", row.names = F, quote = F)

table <- importFromTable("soybean.txt",
                      highBulk = HB,
                      lowBulk = LB,
                      chromList = unique(ltable$CHROM),
                      sep = "\t")


table <- subset(table, !is.na(SNPindex.LOW) & !is.na(SNPindex.HIGH))



## G' stats

table <- runGprimeAnalysis(SNPset = table,
                        windowSize = 1e6,
                        outlierFilter = "deltaSNP")

### plot result ###
plotQTLStats(
    SNPset = table ,
    var = "Gprime",
    plotThreshold = TRUE,
    q = 0.01
)

### QTL-seq ###

table <- runQTLseqAnalysis(SNPset = table,
                        windowSize = 1e6,
                        popStruc = "F2",
                        bulkSize = c(30,30))



### plot result ###
plotQTLStats(
    SNPset = table,
    var = "deltaSNP",
    plotIntervals = TRUE)


##### Result VZ ######

library(CMplot)





