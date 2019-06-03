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

# Load data
data("mtcars")
df <- mtcars

# Convert cyl as a grouping variable
df$cyl <- as.factor(df$cyl)

# Inspect the data
head(df[, c("wt", "mpg", "cyl", "qsec")], 4)

ggplot(df, aes(x = wt, y = mpg)) +
  geom_point(aes(color = cyl, size = qsec), alpha = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  scale_size(range = c(0.5, 12)) # Adjust the range of points size

# Color and shape depend on factor (categorical variable)
ggplot(iris, aes(
  x = Sepal.Length, y = Sepal.Width, color = Petal.Length,
  size = Petal.Length
)) +
  geom_point(alpha = 0.6)

# setting data director

go <- read.csv("goenrich.csv", header = T)

head(go)

ggplot(go, aes(
  x = richfactor, y = Description, color = Petal.Length,
  size = Petal.Length
)) +
  geom_point(alpha = 0.6)
tw <- go[c(1:20), ]
ggplot(tw, aes(x = richfactor, y = Description)) +
  geom_point(aes(color = Pvalue_corrected, size = Number), alpha = 0.5) # +
# scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
# scale_size(range = c(0.5, 12))  # Adjust the range of points size

# for finding the font in our computer
windowsFonts()
## BUDDLE PLOT
ggplot(tw, aes(x = richfactor, y = Description)) +
  geom_point(aes(color = Pvalue_corrected, size = Number), alpha = 1) +
  scale_color_gradient(low = "black", high = "red") +
  theme(
    legend.position = "right",
    axis.text = element_text(
      family = "serif",
      size = 10, face = "bold"
    ),
    legend.text = element_text(
      family = "serif",
      size = 10, face = "bold"
    ),
    legend.title = element_text(
      family = "serif",
      size = 10, face = "bold"
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  labs(size = "Number of Genes", colour = "FDR") +
  ylab(NULL) + xlab("Rich Factor")

#guides(fill = guide_legend(title = c("FDR", "number genes")))
#theme(
#  plot.title = element_text(color = "red", size = 14, face = "bold.italic"),
#  axis.title.x = element_text(color = "blue", size = 14, face = "bold"),
#  axis.title.y = element_text(color = "#993333", size = 14, face = "bold")
#)

# color
# Default version: just say color=your numeric column
ggplot(mtcars, aes(x = wt, y = mpg, color = disp)) + geom_point(size = 4)

# 1 - Scale_fill_gradient
ggplot(mtcars, aes(x = wt, y = mpg, color = disp)) + geom_point(size = 5) +
  scale_color_gradient(low = "black", high = "red")

# 2 - Three colours gradient
ggplot(mtcars, aes(x = wt, y = mpg, color = disp)) + geom_point(size = 5) +
  scale_fill_gradient2(midpoint = 3)

# 3 - scale_colour_gradientn() and scale_fill_gradientn(): a custom n-colour gradient.
# Also works with cm.colors, heat.colors, and the colors of the package "colorspace"
ggplot(mtcars, aes(x = wt, y = mpg, color = disp)) + geom_point(size = 5) +
  scale_color_gradientn(colours = terrain.colors(7))

# 4 - Using R colorbrewer
ggplot(mtcars, aes(x = wt, y = mpg, color = disp)) + geom_point(size = 5) +
  scale_color_distiller(palette = "RdPu")
