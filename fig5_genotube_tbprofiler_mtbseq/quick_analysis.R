{
  rm(list = ls())
  gc()
  graphics.off()
}


deletion <- read.table("deletion_distribution.txt")
median(sapply(unique(deletion$V1), function(x) sum(deletion$V2[deletion$V1 == x])))

insertion <- read.table("insertion_distribution.txt")
summary(sapply(unique(insertion$V1), function(x) sum(insertion$V2[insertion$V1 == x])))


readLines("test")
