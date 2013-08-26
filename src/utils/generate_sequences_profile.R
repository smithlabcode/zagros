library(lattice)
ltheme <- canonical.theme(color = FALSE) 
ltheme$strip.background$col <- "transparent" 
lattice.options(default.theme = ltheme) 
args = commandArgs()
s = substr(args[3],1,nchar(args[3])-4)
pdf(paste(s,".pdf",sep = ""), width=14, height=7*as.integer(as.numeric(args[5])/30))
x1 <- read.table(args[3])$V1
x2 <- read.table(args[3])$V2
x3 <- read.table(args[3])$V3
y <- matrix(x1, as.numeric(args[6]), as.numeric(args[5]), dimnames=list(c(1:as.numeric(args[6])), c(as.numeric(args[5]):1)))
text <- matrix(x2, as.numeric(args[6]), as.numeric(args[5]), dimnames=list(c(1:as.numeric(args[6])), c(as.numeric(args[5]):1)))
color <- matrix(x3, as.numeric(args[6]), as.numeric(args[5]), dimnames=list(c(1:as.numeric(args[6])), c(as.numeric(args[5]):1)))
z <- levelplot(y, aspect="fill", xlab="", ylab="", colorkey = FALSE, par.settings = list(axis.line = list(col = 0)),scales=list( x = list(draw = FALSE)),
               panel=function(...) {
               arg <- list(...)
               panel.levelplot(...)
               panel.text(arg$x, arg$y, text, col=color)})
print(z)
dev.off()

