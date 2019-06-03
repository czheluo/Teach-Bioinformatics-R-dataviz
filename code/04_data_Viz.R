##########################################
#  http://www.majorbio.com/
#  Copyright (C) 2019 in Majorbio workshop
#  contact: meng.luo@majorbio.com
##########################################

# Scatterplots (plot(...),scatter(...))

library(datasets)
data(iris)
str(iris)
plot(iris$Sepal.Width, 
     iris$Sepal.Length)

# compare

library(ggplot2)
plot(iris$Sepal.Width,  iris$Sepal.Length)
qplot(Sepal.Width, Sepal.Length, data = iris)
qplot(Sepal.Width, Sepal.Length, data = iris, color = Species, shape = Species) + theme_bw() # theme_classic() 
qplot(Sepal.Width, Sepal.Length, data = iris, color = Species, shape = Species, geom = c("point", "smooth"), method = lm, se = FALSE) + theme_bw()


# A basic scatterplot
plot(x = 1:10,
     y = 1:10,
     xlab = "X Axis label",
     ylab = "Y Axis label",
     main = "Main Title")

# 3D Scatterplot with Coloring and Vertical Lines and Regression Plane
library(scatterplot3d)
attach(mtcars)
s3d <-scatterplot3d(wt,disp,mpg, pch=16, highlight.3d=TRUE,
                    type="h", main="3D Scatterplot")
fit <- lm(mpg ~ wt+disp)
s3d$plane3d(fit)


# Lines (plot(...), lines(...))
x <- c(1:5); y <- x # create some data
par(pch=22, col="red") # plotting symbol and color
par(mfrow=c(2,4)) # all plots on one page
opts = c("p","l","o","b","c","s","S","h")
for(i in 1:length(opts)){
  heading = paste("type=",opts[i])
  plot(x, y, type="n", main=heading)
  lines(x, y, type=opts[i])
}

# barplot (barplot(.))
#barplot(height, width = 1, space = NULL,
#        names.arg = NULL, legend.text = NULL,
#        horiz = FALSE, density = NULL, angle = 45,
#        col = NULL, border = par("fg"),
#        main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
#        xlim = NULL, ylim = NULL, xpd = TRUE, log = "",
#        axes = TRUE, axisnames = TRUE,
#        cex.axis = par("cex.axis"), cex.names = par("cex.axis"),
#        plot = TRUE, axis.lty = 0, .)

# Simple Bar Plot

counts <- table(mtcars$gear)
barplot(counts, main="Car Distribution",
        xlab="Number of Gears")

# Stacked Bar Plot with Colors and Legend
counts <- table(mtcars$vs, mtcars$gear)
barplot(counts, main="Car Distribution by Gears and VS",
        xlab="Number of Gears", col=c("darkblue","red"),
        legend = rownames(counts)) 


# Grouped Bar Plot
counts <- table(mtcars$vs, mtcars$gear)
barplot(counts, main="Car Distribution by Gears and VS",
        xlab="Number of Gears", col=c("darkblue","red"),
        legend = rownames(counts), beside=TRUE)

# Histograms(hist(.))
#hist(x, breaks = "Sturges",
#     freq = NULL, probability = !freq,density = NULL, col = NULL, border = NULL,
#     main = paste("Histogram of" , xname), xlim = range(breaks), ylim = NULL,
#     xlab = xname, ylab,
#     axes = TRUE, plot = TRUE, labels = FALSE, .)

# Simple Histogram
hist(mtcars$mpg)

# Colored Histogram with Different Number of Bins
hist(mtcars$mpg, breaks=12, col="red") 


# Add a Normal Curve
x <- mtcars$mpg
h<-hist(x, breaks=10, col="red", xlab="Miles Per Gallon",
        main="Histogram with Normal Curve")
xfit<-seq(min(x),max(x),length=40)
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=2) 

# Density Plot (density(.))
# Kernel Density Plot
d <- density(mtcars$mpg) # returns the density data
plot(d) # plots the results 

# Filled Density Plot
d <- density(mtcars$mpg)
plot(d, main="Kernel Density of Miles Per Gallon")
polygon(d, col="red", border="blue") 

# boxplot (boxplot(.))
# S3 method for formula
#boxplot(formula, data = NULL, ., subset, na.action = NULL,
#        xlab = paste(names(mf)[-response], collapse = " : "),
#        ylab = names(mf)[ response],
#        add = FALSE, ann = !add,
#        drop = FALSE, sep = ".", lex.order = FALSE)

# Boxplot of MPG by Car Cylinders
boxplot(mpg~cyl,data=mtcars, main="Car Milage Data",
        xlab="Number of Cylinders", ylab="Miles Per Gallon") 


# Notched Boxplot of Tooth Growth Against 2 Crossed Factors
# boxes colored for ease of interpretation
boxplot(len~supp*dose, data=ToothGrowth, notch=TRUE,
        col=(c("gold","darkgreen")),
        main="Tooth Growth", xlab="Suppliment and Dose") 


# Violin Plots
# install.packages('vioplot')
library(vioplot)
x1 <- mtcars$mpg[mtcars$cyl==4]
x2 <- mtcars$mpg[mtcars$cyl==6]
x3 <- mtcars$mpg[mtcars$cyl==8]
vioplot(x1, x2, x3, names=c("4 cyl", "6 cyl", "8 cyl"),
        col="gold")
title("Violin Plots of Miles Per Gallon")


# Viz large datasets
library (ggplot2)
train <- read.csv('Big_Mart_Dataset.csv',header = T)
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


