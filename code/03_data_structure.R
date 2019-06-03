##########################################
#  http://www.majorbio.com/
#  Copyright (C) 2019 in Majorbio workshop
#  contact: meng.luo@majorbio.com
##########################################

# Data Type

# character
v <- c('a' , "good", "TRUE", '23.4')
print(class(v))

# numeric ()  real or decimal
v <- c(23.4, 36.6, 63.3)
print(class(v))

# integer
v <- 2L
print(class(v))

# missing value
i <- NA # (NULL)

# logical
v <- TRUE
print(class(v))

# complex
v <- 3 + 2i
print(class(v))

# raw
v <- charToRaw("Hello")
print(class(v))


# Data Structure

# vector
vec <- c(1, 2, 3, 10, 100)
apple <- c('red','green',"yellow")
print(vec)
print(apple)

# matrix
mat1 <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2)
mat2 <- matrix(c('a','a','b','c','b','a'), nrow = 2, ncol = 3, byrow = TRUE)
print(mat1)
print(mat2)
                
# math
vec + 4 
vec * 4
vec * vec # element-wise multiplication

vec %*% vec # matrix multiplication (inner)

vec %o% vec # matrix multiplication (outer)

t(mat1) # transpose

# array
a <- array(c('green','yellow'),dim = c(3,3,2))
a <- array(c(1,2,3),dim = c(3,3,2))

print(a)


# data frames
BMI <- data.frame(
  gender = c("Male", "Male","Female"), 
  height = c(152, 171.5, 165), 
  weight = c(81,93, 78),
  Age = c(42,38,26)
)
print(BMI)


# factor
# Create a vector
apple.colors <- c('green', 'green', 'yellow', 'red', 'red ','red', 'green')

# Create a factor object
factor.apple <- factor(apple.colors)

# Print the factor
print(factor.apple)
print(nlevels(factor.apple))

# list
# Create a list.
list1 <- list(c(2,5,3), 21.3, sin)
# Print the list.
print(list1)

# Subsetting Data 

# Three types
#integers
apple_colors <- c('green', 'green', 'yellow', 'red', 'red ','red', 'green')
apple.colors[c(1,2)]
apple.colors[-c(1,2)]
#name
names(apple.colors)<-c(letters[c(1:7)])
print(apple.colors)
apple.colors[c('a','b')]

# logicals
apple.colors[c(TRUE,FALSE,TRUE,FALSE,TRUE,FALSE,TRUE)]

# Import Data 

#readtable
df1 <- read.table("https://s3.amazonaws.com/assets.datacamp.com/blog_assets/test.txt", 
                  header = FALSE)
# Read in csv files
df2 <- read.table("https://s3.amazonaws.com/assets.datacamp.com/blog_assets/test.csv", 
                  header = FALSE,
                  sep = ",")

df3 <- read.csv("https://s3.amazonaws.com/assets.datacamp.com/blog_assets/test.csv",
                header = FALSE)

df4 <- read.csv2("https://s3.amazonaws.com/assets.datacamp.com/blog_assets/test.csv", 
                 header= FALSE)

# Read a delimited file
df1 <- read.delim("https://s3.amazonaws.com/assets.datacamp.com/blog_assets/test_delim.txt", sep="$") 
df2 <- read.delim2("https://s3.amazonaws.com/assets.datacamp.com/blog_assets/test_delim.txt", sep="$")

# Read a .xls or .xlsx files
library(readxl)
df1 <- read_excel("<name and extension of your file>")

# Read RDA or RData Files into R
load('name.RData')
load('name.RDA')


#write .csv file
write.csv(df1, file= 'df1.csv',quote=F)


#write  .CSV or txt files
write.table(df1,file='df1.csv',sep=',',quote=F)
write.table(df1,file='df1.txt',sep=',',quote=F)


#save RDA or RData Files into R
save(df1,file='name.RData')
save.image(df1,file='name.RData')




