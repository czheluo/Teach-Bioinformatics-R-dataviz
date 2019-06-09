##########################################
#  http://www.majorbio.com/
#  Copyright (C) 2019 in  Majorbio workshop
#  contact: meng.luo@majorbio.com
##########################################

set.seed(0613)

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
# setwd()

go <- read.csv("goenrich.csv", header = T)

head(go)

tw <- go[sample(c(1:25)), ]

# for finding the font in our computer
# windowsFonts()

## BUDDLE PLOT
p <- ggplot(tw, aes(x = richfactor, y = Description, color = Pvalue_corrected, size = Number)) +
  geom_point(alpha = 1) +
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
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    # panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  labs(size = "Number of Genes", colour = "FDR") +
  ylab(NULL) + xlab("Rich Factor")

# facet
facet_grid(~Term.Type, scales = "free_y")

# save image
png(paste("bubble", ".png", sep = ""), width = 1000, height = 800)
print(p)
dev.off()

p <- ggplot(tw, aes(x = richfactor, y = Description, color = Pvalue_corrected, size = Number)) +
  geom_point(alpha = 1) +
  # scale_color_gradient(low = "black", high = "red") +
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

pdf("bubble.pdf", width = 10, height = 10)

print(p)

dev.off()

# CLEAN DATA
rm(list = ls(all.names = T))

# clean image

graphics.off()
