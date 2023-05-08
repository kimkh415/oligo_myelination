# convert slide-seq raw data to matrix mtx format



library(DropletUtils)
library(Matrix)
library(data.table)
library(dplyr)


exp_file = "Puck_200306_02.digital_expression.txt"

count = as.data.frame(fread(exp_file, header=T))
rownames(count) = count[,1]
count = count[,-1]
count = as.matrix(count)
x = Matrix(count, sparse=T)

write10xCounts(path=".", x=x, overwrite=T)
#system("gzip data/*")

print('done')

