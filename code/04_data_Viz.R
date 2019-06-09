##########################################
#  http://www.majorbio.com/
#  Copyright (C) 2019 in Majorbio workshop
#  contact: meng.luo@majorbio.com
##########################################

# load packages
library(ggplot2)
library(datasets)

# set data dir

setwd()

data(iris)
str(iris)

# set reproduced 

set.seed(0613)

# quick plot and ggplot

qplot(Sepal.Width, Sepal.Length, data = iris)

ggplot(iris,aes(x=Sepal.Width, y=Sepal.Length)) +geom_point()

qplot(Sepal.Width, Sepal.Length, data = iris, 
      color = Species, shape = Species) + 
  theme_bw() # theme_classic() 

ggplot(iris,aes(x=Sepal.Width, y=Sepal.Length, color = Species)) +
  geom_point() + 
  theme_bw()

qplot(Sepal.Width, Sepal.Length, data = iris, 
      color = Species, shape = Species, 
      geom = c("point", "smooth"), 
      method = lm, se = FALSE) + theme_bw()

ggplot(iris,aes(x=Sepal.Width, y=Sepal.Length,  color = Species)) +
  geom_point() + geom_smooth(method = lm)+
  theme_bw()


#facet

ggplot(iris,aes(x=Sepal.Width, y=Sepal.Length,  color = Species)) +
  geom_point() + geom_smooth(method = lm)+
  theme_bw()+facet_wrap(~Species,dir = "h")

ggplot(iris,aes(x=Sepal.Width, y=Sepal.Length,  color = Species)) +
  geom_point() + geom_smooth(method = lm)+
  theme_bw()+facet_grid(~Species)



# Viz large datasets
# read dataset
train <- read.csv('Big_Mart_Dataset.csv',header = T)

# plot
ggplot(train, aes(Item_Visibility, Item_MRP)) + geom_point(aes(color = Item_Type)) + 
  scale_x_continuous("Item Visibility", breaks = seq(0,0.35,0.05))+
  scale_y_continuous("Item MRP", breaks = seq(0,270,by = 30))+ 
  theme_bw() + labs(title="Scatterplot") + facet_wrap( ~ Item_Type)


## boxplot 
ggplot(train, aes(Outlet_Identifier, Item_Outlet_Sales)) + geom_boxplot(fill = "red")+
  scale_y_continuous("Item Outlet Sales", breaks= seq(0,15000, by=500))+
  labs(title = "Box Plot", x = "Outlet Identifier")



# Open a pdf file
pdf("rplot.pdf") 
# Create a plot
ggplot(train, aes(Outlet_Identifier, Item_Outlet_Sales)) + geom_boxplot(fill = "red")+
  scale_y_continuous("Item Outlet Sales", breaks= seq(0,15000, by=500))+
  labs(title = "Box Plot", x = "Outlet Identifier")
# Close the pdf file
dev.off() 



# Open a pdf file
png("rplot.png") 
# Create a plot
ggplot(train, aes(Outlet_Identifier, Item_Outlet_Sales)) + geom_boxplot(fill = "red")+
  scale_y_continuous("Item Outlet Sales", breaks= seq(0,15000, by=500))+
  labs(title = "Box Plot", x = "Outlet Identifier")
# Close the pdf file
dev.off() 

# Clean data

rm(list = ls(all.names = T))
