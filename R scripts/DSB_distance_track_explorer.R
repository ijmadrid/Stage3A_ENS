res <- read.csv("./Stage 3A ENS/DNA_Religation/projects/RCLPolymer_TimesBelowThreshold/results/2018_06_11_18_40_trackDSB_fixedDSBs.csv")

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

amps <- res[c("Nc","a1_Amplitude","a2_Amplitude","b1_Amplitude","b2_Amplitude")]
amps.melted <- melt(amps, id = "Nc")

alphas <- res[c("Nc","a1_Alpha","a2_Alpha","b1_Alpha","b2_Alpha")]
alphas.melted <- melt(alphas, id = "Nc")

ggplot(data = alphas.melted, aes(x = Nc, y = value, color = variable, shape = variable)) + geom_point(size = 3) + geom_line(size = 1)


#############################
#############################

ggplot(data = res, aes(x = Nc, y = estimatedProba_byMeanTimes)) + geom_line()



#############################

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

exp.fit <- function(df, k){
  return(exponential(df$realtime, res$a1.a2_ComebackAmplitude[k] + res$b1.b2_ComebackAmplitude[k], 
              res$a1.a2_ComebackRate[k] + res$b1.b2_ComebackRate[k] ))
}

exp.fit <- function(df, k){
  min.rate = res$a1.a2_ComebackRate[k] + res$a1.b1_ComebackRate[k] +
    res$a1.b2_ComebackRate[k] + res$a2.b1_ComebackRate[k] +
    res$a2.b2_ComebackRate[k] + res$b1.b2_ComebackRate[k]
  
  min.ampl = min.rate * 
    res$a1.a2_ComebackAmplitude[k]/res$a1.a2_ComebackRate[k] *
    res$a1.b1_ComebackAmplitude[k]/res$a1.b1_ComebackRate[k] *
    res$a1.b2_ComebackAmplitude[k]/res$a1.b2_ComebackRate[k] *
    res$a2.b1_ComebackAmplitude[k]/res$a2.b1_ComebackRate[k] *
    res$a2.b2_ComebackAmplitude[k]/res$a2.b2_ComebackRate[k] *
    res$b1.b2_ComebackAmplitude[k]/res$b1.b2_ComebackRate[k]
  
  min.ampl = min.rate
  
  return(exponential(df$realtime, min.ampl, min.rate))
}

{df <- data.frame(realtime = seq(1, 2, 0.05))
df$"3"  <- exp.fit(df,1)
df$"4"  <- exp.fit(df,2)
df$"5"  <- exp.fit(df,3)
df$"6"  <- exp.fit(df,4)
df$"7"  <- exp.fit(df,5)
df$"8"  <- exp.fit(df,6)
df$"9"  <- exp.fit(df,7)
df$"10" <- exp.fit(df,8)
df$"11" <- exp.fit(df,9)
df$"12" <- exp.fit(df,10)
df$"13" <- exp.fit(df,11)
df$"14" <- exp.fit(df,12)
df$"15" <- exp.fit(df,13)
df$"16" <- exp.fit(df,14)
df$"17" <- exp.fit(df,15)
df$"18" <- exp.fit(df,16)
df$"19" <- exp.fit(df,17)
df$"20" <- exp.fit(df,18)
df$"21" <- exp.fit(df,19)
df$"22" <- exp.fit(df,20)
df$"23" <- exp.fit(df,21)
}

library(reshape)
df.melted <- melt(df, id = c("realtime"), variable_name = "nc")
df.melted$nc <- as.numeric(df.melted$nc)
ggplot(data = df.melted, aes(x = realtime, y = value, color = nc)) + geom_point() + scale_color_gradientn(colours = heat.colors(20))
