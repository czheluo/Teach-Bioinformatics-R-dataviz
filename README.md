###  WORKSHOP Overview <a href=""><img src="Fig/outline.png" align="right" alt="logo" height="157" width="200" /></a>

### Big Data and Data Science using R for researchers 
 This is a short course for Rstat and Rdataviz. Mainly is designed to provide a good opportunity for researchers to learn R (An interactive approach to statistical computing). Specifically designed for data analysis and graphics (ggplot2) and the visual analysis of the results related to transcriptome analysis (volcano plot, bubble plot, complexheatmap (go enrichment Set analysis, kegg analysis (pathsway analysis)) and venn graph). Statistical genomics (biometric models) to be included in later courses (eg: GWAS (EWAS and TWAS)). Looking forward to meetting you at majorbio in shanghai. 
 
 ### Location
 		<br>
		<iframe src="https://www.google.com/maps/embed?pb=!1m18!1m12!1m3!1d3416.2558021437867!2d121.62565031548901!3d31.102621574761713!2m3!1f0!2f0!3f0!3m2!1i1024!2i768!4f13.1!3m3!1m2!1s0x35b278a7530eb1dd%3A0xb40b78f5079ac68b!2sCentury+Medicine+Park!5e0!3m2!1sen!2sus!4v1531624183346" width="600" height="450" frameborder="0" style="border:0" allowfullscreen></iframe>
		<br>
		<br>

### Schedule

| Topic | Time | Day|
| :---: | :---: | :---: |
| Introduction to Linux for bioinformatics | 08:30 - 12:00 | 06.13
| DataScience For Epigenetics | 13:30 - 17:30 | 06.13
| Transcriptome analysis procedure I (QC, Assembly, Mapping and RSEM) | 08:30 - 12:00 | 06.14
| Transcriptome analysis procedure II (DE(edgeR), PCA AND (GO, KEGG) set analysis)  | 13:30 - 17:30 | 06.14

### Contact Me:
> [meng.luo@majorbio.com](Meng Luo) OR [czheluo@gmail.com](Chenzhe Luo) 

# Installation
 
**Required packages** can be installed on Windows, Linux and MacOS(NO TRAYING) with following steps:

**installation for newest R packages**
Required packages can be installed with following R codes for Windows:  
```r
>source(install_packages.R)
## Offline installation
>sourceInstall(path = setwd()) # the PATH WAS THE PACKAGES LOACATION
## online installation
>sourceInstall(win = T, packages = T) # Recommended
```
Required packages can be installed with following R codes for Linux:  
```r
>source(install_packages.R)
## Offline installation
>sourceInstall(path = setwd()) # the PATH WAS THE PACKAGES LOACATION
## online installation
>sourceInstall(linux = T, packages = T) # Recommended
```
