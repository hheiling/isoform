
library(MASS)
library(isoform)

# --------------------------------------------------------------
# we simulate a gene of 5 exons with three isoforms:
# isoform 1: 1-2-3-4-5
# isoform 2: 1-2-5
# isoform 3: 2-3-5
# suppose the effective length of each exon is 100bp
# and the effecitve length of each exon junction is 50bp

# include three isoforms that are not expressed 
# isoform 1-4-5
# isoform 1-2-3-4
# isoform 1-2-4-5
# --------------------------------------------------------------

# for each isoform, there are 15 observations
# exons 1-5, exon juctions 1-2, 2-3, 3-4, 4-5, 
# and 1-3, 1-4, 1-5, 2-4, 2-5, and 3-5

x1 = c(rep(100,5), rep(50,4), rep(0,6))

x2 = c(100, 100, 0, 0, 100, 50, 0, 0, 0, 0, 0, 0, 0, 50, 0)
x3 = c(0, 100, 100, 0, 100, 0, 50, 0, 0, 0, 0, 0, 0, 0, 50)
x4 = c(100, 0, 0, 100, 100, 0, 0, 0, 50, 0, 50, 0, 0, 0, 0)
x5 = c(100, 100, 100, 100, 0, 50, 50, 50, 0, 0, 0, 0, 0, 0, 0)
x6 = c(100, 100, 0, 100, 100, 50, 0, 0, 50, 0, 0, 0, 50, 0, 0)

X  = cbind(x1, x2, x3, x4, x5, x6)
rnames = c(1:5, "1-2", "2-3", "3-4", "4-5", "1-3", "1-4", "1-5", "2-4", "2-5", "3-5")
rownames(X) = c(paste("exon", rnames, sep=""))
X

mu = 2*rowSums(X[,1:3])
set.seed(2002)
y  = rnegbin(mu, theta=5)

is = isoSelect(y, X)
is

mu = X[,is$w2use] %*% is$b2use
plot(mu, y)
abline(0,1)

simuE5I6 = data.frame(y=y, X=X)
names(simuE5I6) = c("y", "x1", "x2", "x3", "x4", "x5", "x6")

save("simuE5I6", file="simuE5I6.rda")
