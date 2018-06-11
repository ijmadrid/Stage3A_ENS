res <- read.csv("./Stage 3A ENS/DNA_Religation/projects/RCLPolymer_TimesBelowThreshold/results/2018_06_10_04_08_trackDSB_fixedDSBs.csv")

res[is.na(res)] <- Inf

library(ggplot2)

res$genomicDistance <- res$B1 - res$A1 - 1

{p <- ggplot(data = res, aes(x = genomicDistance, y = Nc, fill = a1.b2_ComebackRate)) + geom_tile()
p}

res$genomicDistance <- res$B1 - res$A1 - 1
res$genomicDistance <- as.factor(res$genomicDistance)
ggplot(data = res, aes(x = Nc, y = a1.a2_ComebackRate, color = genomicDistance)) + geom_line(size=2)


ggplot(data = res, aes(x = Nc, y = a1.b2_meanComebackTime, color = genomicDistance)) + geom_line() +
  geom_ribbon(aes(x = Nc, 
                  ymin = a1.b2_meanComebackTime - a1.b2_meanComebackTime_dx,
                  ymax = a1.b2_meanComebackTime + a1.b2_meanComebackTime_dx), 
              alpha = 0.2)

res$misrepairRate <- (1/res$a1.b1_meanComebackTime + 1/res$a1.b2_meanComebackTime + 
                      1/res$a2.b1_meanComebackTime + 1/res$a2.b2_meanComebackTime)

res$repairRate <- 1/res$a1.a2_meanComebackTime + 1/res$b1.b2_meanComebackTime

res$repairComeback <- (res$a1.a2_meanComebackTime + res$b1.b2_meanComebackTime)/2

res$encounterRate <- res$misrepairRate + res$repairRate

ggplot(data = res, aes(x = Nc, y = repairComeback,  color = genomicDistance)) + geom_line()
  
ggplot(data = res, aes(x = Nc, y = a1.b1_ComebackRate +
                                   a1.b2_ComebackRate + 
                                   a2.b1_ComebackRate +
                                   a2.b2_ComebackRate,  color = genomicDistance)) + geom_line()

res$EstimatedProba <- (res$a1.a2_ComebackRate + res$b1.b2_ComebackRate)/
                      (res$a1.b1_ComebackRate + res$a1.b2_ComebackRate + 
                       res$a2.b1_ComebackRate + res$a2.b2_ComebackRate + 
                       res$a1.a2_ComebackRate + res$b1.b2_ComebackRate)

res$EstimatedProba <- (1/res$a1.a2_meanComebackTime + 1/res$b1.b2_meanComebackTime)/
                      (1/res$a1.b1_meanComebackTime + 1/res$a1.b2_meanComebackTime + 
                       1/res$a2.b1_meanComebackTime + 1/res$a2.b2_meanComebackTime + 
                       1/res$a1.a2_meanComebackTime + 1/res$b1.b2_meanComebackTime)

res$EstimatedProba <- ()


alphas <- res[c("Nc","polymerAlpha","centerAlpha","a1_Alpha","a2_Alpha","b1_Alpha","b2_Alpha")]
alphas.melted <- melt(alphas, id = "Nc")

ggplot(data = alphas.melted, aes(x = Nc, y = value, color = variable)) + geom_point() + geom_line()


variances <- res[c("Nc",
                  "a1.a2_interbreakDistance_meanOfVariances",
                  "a1.b1_interbreakDistance_meanOfVariances",
                  "a1.b2_interbreakDistance_meanOfVariances",
                  "a2.b1_interbreakDistance_meanOfVariances",
                  "a2.b2_interbreakDistance_meanOfVariances",
                  "b1.b2_interbreakDistance_meanOfVariances")]

variances.melted <- melt(variances, id = "Nc")

ggplot(data = variances.melted, aes(x = Nc, y = value, color = variable)) + geom_line(size = 1) + geom_point(size = 3)

ggplot(data = res, aes(x = Nc, y = a1.a2_interbreakDistance_totalVariance, color = genomicDistance)) + geom_point() + geom_line()


# Numeric vector from string lists
library(gsubfn)

msds <- res$a1_MSD
str2numvector <- function(strlist){
  z <- strapply(as.character(strlist), "(\\d*\\.*\\d+|\\d+\\.|\\d*\\.*\\d+e\\+\\d+|\\d*\\.*\\d+e\\-\\d+|nan)", as.numeric)
  z <- cbind(sapply(z, 2, FUN = unlist))
  return(z)
}

msds.num <- str2numvector(msds)
timeline <- c(1,2,3,29997,29998,29999)

msd.data <- data.frame(msd = msds.num[,1], time = timeline)
ggplot(data = msd.data, aes(x= time, y = msd)) + geom_line()
# Histograms

exponential <- function(t,A,l){
  return(A*exp(-l*t))
}

{df <- data.frame(realtime = seq(1, 2, 0.1))
df$y1 <- exponential(df$realtime, res$a1.a2_ComebackAmplitude[3] + res$b1.b2_ComebackAmplitude[3], 
                     res$a1.a2_ComebackRate[3] + res$b1.b2_ComebackRate[3] )
df$y2 <- exponential(df$realtime, res$a1.a2_ComebackAmplitude[6] + res$b1.b2_ComebackAmplitude[6], 
                     res$a1.a2_ComebackRate[6] + res$b1.b2_ComebackRate[6] )
df$y3 <- exponential(df$realtime, res$a1.a2_ComebackAmplitude[9] + res$b1.b2_ComebackAmplitude[9], 
                     res$a1.a2_ComebackRate[9] + res$b1.b2_ComebackRate[9] )
df$y4 <- exponential(df$realtime, res$a1.a2_ComebackAmplitude[12] + res$b1.b2_ComebackAmplitude[12], 
                     res$a1.a2_ComebackRate[12] + res$b1.b2_ComebackRate[12] )
df$y5 <- exponential(df$realtime, res$a1.a2_ComebackAmplitude[15] + res$b1.b2_ComebackAmplitude[15], 
                     res$a1.a2_ComebackRate[15] + res$b1.b2_ComebackRate[15] )
df$y6 <- exponential(df$realtime, res$a1.a2_ComebackAmplitude[18] + res$b1.b2_ComebackAmplitude[18], 
                     res$a1.a2_ComebackRate[18] + res$b1.b2_ComebackRate[18] )
df$y7 <- exponential(df$realtime, res$a1.a2_ComebackAmplitude[21] + res$b1.b2_ComebackAmplitude[21], 
                     res$a1.a2_ComebackRate[21] + res$b1.b2_ComebackRate[21] )
df$y8 <- exponential(df$realtime, res$a1.a2_ComebackAmplitude[24] + res$b1.b2_ComebackAmplitude[24], 
                     res$a1.a2_ComebackRate[24] + res$b1.b2_ComebackRate[24] )
}
library(reshape)
df.melted <- melt(df, id = c("realtime"))
ggplot(data = df.melted, aes(x = realtime, y = value, color = variable)) + geom_line(size = 1)
