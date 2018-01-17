#
# Example from OECD 2006
#
# The data set 'trout' contains
# dose levels of 10, 25, 60, 150 and 1000 ppm,
# and a control as well as
# the response in trout weight in mg
# 

data(trout)
attach(trout)
xmean <- tapply(Y, DOSE, mean)
xn <- tapply(Y, DOSE, length)
xsd <- tapply(Y, DOSE, sd)
xse <- xsd / sqrt(xn)

ans <- data.frame(MEAN = round(xmean,1), SE = round(xse, 3), n = xn)
rownames(ans) <- levels(DOSE)

# 
ans

#
# Check for normality  
#
fit <- aov(Y ~ DOSE - 1)
shapiro.test(residuals(fit))

##
## Check for homogeneous variances
##
require(car)
leveneTest(fit)

##
## Perform Tanhame-Dunnett test
##
td.out <- tamhaneDunnettTest(Y, DOSE, alternative = 'less')
summary(td.out)

##
## Perform Dunnett test
##
summary(dunnettTest(Y, DOSE, alternative = 'less'))
require(multcomp)
summary(glht(fit, linfct=mcp(DOSE = "Dunnett"), alternative = 'less'))

#
# Test for monotonicity
# normalised ranks
#
k <- length(xn)
mat <- contr.poly(k)
mat <- mat[,1:2]
colnames(mat) <- c("Linear trend", "Quadratic trend")
mat <- t(mat)

m <- length(Y)
Rij <- rank(Y)
Rn <- Rij / m
Rfit <- aov(Rn ~ DOSE)
summary(glht(Rfit, linfct = mcp(DOSE = mat)))

##
## Perform step-down Jonckheere test
##
pval <- rep(NA, k)
z <- rep(NA, k)
for (i in k:2){
  YY <- Y[as.numeric(DOSE) <= i]
  DDOSE <- DOSE[as.numeric(DOSE) <= i]
  js.out <- jonckheereTest(YY, DDOSE, alternative = 'less')
  z[i] <- js.out$statistic
  pval[i] <- js.out$p.value
}

out <- data.frame(MEAN = round(xmean, 1), z = round(z, 3), p = round(pval, 3))
rownames(out) <- levels(DOSE)
out

#
# Employ pairwise Wilcox test with Holm adjustment
#
w.out <- pairwise.wilcox.test(Y, DOSE, alternative = 'less', p.adj= 'none')
pval <- as.vector(w.out$p.value[,1])
padj <- p.adjust(pval, 'holm')
lev <- levels(DOSE)
H0 <- rep(NA, k-1)
for (i in 2:k){
  H0[i-1] <- paste(lev[i], "-", lev[1], ">= 0")
}
wb.out <- data.frame(p = pval, padj = padj)
rownames(wb.out) <- H0
wb.out

