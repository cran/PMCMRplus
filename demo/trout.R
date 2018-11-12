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
H0 <- rep(NA, k)
z <- rep(NA, k)
for (i in k:2){
  YY <- Y[as.numeric(DOSE) <= i]
  DDOSE <- DOSE[as.numeric(DOSE) <= i]
  js.out <- jonckheereTest(YY, DDOSE, alternative = 'less')
  z[i] <- js.out$statistic
  pval[i] <- js.out$p.value

  H0[i] <- ifelse(i == 2,
                  paste0("1 == 2"),
                  paste0("1 == ... == ", i))

  ##H0[i] <- paste0(levels(DOSE)[1:i], collapse = " = ")
}

symp <- symnum(pval[2:k], corr=FALSE,
               cutpoints = c(0,  .001,.01,.05, .1, 1),
               symbols = c("***","**","*","."," "))

out <- data.frame(z = round(z[2:k], 3),
                  p  = format.pval(pval[2:k]),
                  symp)
rownames(out) <- H0[2:k]
names(out) <- c("z", "Pr(>z)", "")

print.ME <- function(x, ...) {
  cat("\n\tStep-down Jonckheere trend test\n\n")
  cat("data: Y and DOSE\n")
  cat("alternative hypothesis: less\n")
  cat("H0\n")
  print(out)
  invisible(x)
}
print.ME(out)
#
# Perform pairwise Wilcox test with Holm adjustment
#
summary(manyOneUTest(Y ~ DOSE, alternative = "less", p.adjust = "holm"))

#
# Perform Williams trend test
#
williamsTest(Y ~ DOSE, alternative = "less")

#
# Perform Shirley's test
#
shirleyWilliamsTest(Y ~ DOSE, alternative = "less")

detach(trout)
