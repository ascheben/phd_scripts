#!/usr/bin/Rscript

# Make linkage disequilibrium (LD) plot
# Input is plink LD result table with R2 values

setwd("/path/to/my/genoTable")
#args <- commandArgs(trailingOnly = TRUE)
#stopifnot(length(args) == 2)

tab <- read.table("distance.ld", sep="\t")

distance <- tab$V1
LD.data <- tab$V2
n<-96
HW.st<-c(C=0.1)
bnap.HW.nonlinear<-nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=100))
# Then to plot this for (for example) positions 1 to 10000:
options(bitmapType='cairo')
png(file="/path/to/output/chrA01.png")
newdist <- seq(1,10000, 10)
plot(newdist, predict(bnap.HW.nonlinear,list(distance=newdist)), type="l", ylim=c(0.0, 0.6), xlab="Distance",ylab="r²")
#legend("topright", legend=c("Brachy","Sub"), lty=c(1,1), col=c("black","red"))
dev.off()
x <- c(1000, 5000, 10000, 500000, 1000000)
y.bnap <- predict(bnap.HW.nonlinear, list(distance=x))
df <- as.data.frame(t(data.frame( y.bnap)))
names(df) <- x
rownames(df) <- c("chrA01bnap")
write.csv(df, "/path/to/output/chrA01.csv")








































