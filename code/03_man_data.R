##########################################
#  http://www.majorbio.com/
#  Copyright (C) 2019 in Majorbio workshop
#  contact: meng.luo@majorbio.com
##########################################

# Import Data
# SET WORK DIRECTOR
getwd()
setwd()
# readtable
df1 <- read.table("group.txt",
  header = FALSE
)

# Read in csv files
df2 <- read.table("group.csv",
  header = FALSE, sep = ","
)

df3 <- read.csv("group.txt",
  header = FALSE, sep = ""
)

df3 <- read.csv("group.csv",
  header = FALSE, sep = ""
)

# save RDA or RData Files into R
save(df1, file = "df1.RData")
save.image(df2, file = "df2.RData")


# Read RDA or RData Files into R
load("df2.RData")
load("df2.RData")


# write .csv file
write.csv(df1, file = "df1.csv", quote = F)


# write  .CSV or txt files
write.table(df1, file = "df1.csv", sep = ",", quote = F)
write.table(df1, file = "df1.txt", sep = ",", quote = F)

# remove variables

rm(df1)
rm(list = ls(all.names = TRUE)) # will clear all objects includes hidden objects
