library(grid)
args = commandArgs()
s = substr(args[3],1,nchar(args[3])-4)
pdf(paste(s,".pdf",sep = ""),width=14,height=4, pagecentre=TRUE)
data = read.table(args[3], header=FALSE)$V1
ksize = 10
toDrop = 0
plot(data, type="l", col="lightblue", lwd=4, main="single stranded preference", ylab="Probability of being unpaired")
dev.off()

