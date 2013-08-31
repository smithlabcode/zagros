library(grid)
args = commandArgs()
s = substr(args[3],1,nchar(args[3])-4)
pdf(paste(s,".pdf",sep = ""),width=14,height=4, pagecentre=TRUE)
data = read.table(args[3], header=FALSE)$V1
hist(data,probability="true", xlim=c(0,160), breaks=160, col="lightblue", main="locations", ylab="Density")
dev.off()

