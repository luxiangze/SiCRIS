setwd("./")
library(ggpubr)

# Create some data format
reads <- read.table("zhangt.gene.raw.count.txt",header=TRUE)

# Basic density plot
# Add mean line and marginal rug
pdf("gene.pdf")
gghistogram(reads, x = "zhangt",
  fill = "lightgray", # 设置填充色
  add = "mean", # 添加均值线
  rug = TRUE # 添加轴须线
)
dev.off()
