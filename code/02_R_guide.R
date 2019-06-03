##########################################
#  http://www.majorbio.com/
#  Copyright (C) 2019 in Majorbio workshop
#  contact: meng.luo@majorbio.com
##########################################


# Workspace
getwd()
setwd("dir") # setwd("G:\\MAJORBIO\\project\\name\\")
save.image("name.Rdata")
savehistory("name.Rhistory")
dir.create(file.path(getwd(),'dirname'), showWarnings = FALSE)
setwd(file.path(getwd(), 'dirname'))

# Clean
rm(Typename)
rm(list = ls(all.names = TRUE)) # will clear all objects includes hidden objects
gc() # free up memrory and report the memory usage
graphics.off() 

