setwd("C:\\Users\\li.pan\\Desktop\\售后汇总\\26年3月\\DZOE2024031499")
data <- read.delim('fpkm.xls',header = T,
                   sep = '\t',
                   check.names = F)
data <- data[,-1]
data <- log2(data+1)
gens <- apply(data,1,sd)
for ( i in seq(0,1,0.1)){
  m <- i
  print(paste0('gens>',m))
  print(sum(gens>m))
}








