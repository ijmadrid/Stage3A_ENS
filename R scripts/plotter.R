res <- read.csv("./Stage 3A ENS/DNA_Religation/projects/RCLPolymer_CLs@DamageFoci/results/2018_05_16_17_52_proba-v_sigma.csv")

library(ggplot2)

res$excludedVolumeCutOff <- as.factor(res$excludedVolumeCutOff)

{
p <- ggplot(data = res, aes(x = excludedVolumeSpringConstant, y = repair_probability)) + geom_line(size = 1.5) + geom_ribbon(data = res, aes(x = excludedVolumeSpringConstant, 
                                  ymin = repair_probability - repair_CIhalflength,
                                  ymax = repair_probability + repair_CIhalflength), alpha = 0.2) + theme_bw()
p
}

{q <- ggplot(data = res, aes(x = excludedVolumeSpringConstant, 
                            y = meanFET)) + geom_line(size = 1.5) + geom_ribbon(aes(x = excludedVolumeSpringConstant, 
                                                                             ymin = meanFET - halfCI_FET,
                                                                             ymax = meanFET + halfCI_FET), alpha = 0.2) + theme_bw()
q+ scale_color_manual(values=c("turquoise1"))}

#### Curve fitting
exp.fit <- lm(data = res, formula = log(meanFET) ~ excludedVolumeSpringConstant)
x <- seq(0,2,0.05)
fit <- exp(predict(exp.fit,list(excludedVolumeSpringConstant=x)))
predicted_df <- data.frame(fit = fit, x=x)
q + geom_line(data = predicted_df, aes(x = x, y = fit), colour = 'red')

{k <- ggplot(data = res, aes(x = excludedVolumeSpringConstant, 
                            y = Ensemble_MSRG, 
                            color = keepCL)) + geom_line(size = 1.5)  + geom_ribbon(aes(x = excludedVolumeSpringConstant, 
                                                                                                     ymin = Ensemble_MSRG - Ensemble_MSRG_dx,
                                                                                                     ymax = Ensemble_MSRG + Ensemble_MSRG_dx), alpha = 0.2) + theme_bw()
k + scale_color_manual(values=c("turquoise1"))}


{
  q <- ggplot(data = res, aes(x = Ensemble_MSRG, 
                             y = repair_probability)) + geom_point(size = 1.5) + theme_bw()

  
  #### fitting
  fit <- lm(data = res, formula = repair_probability ~ poly(sqrt(Ensemble_MSRG),2))
  
  predicted_df = data.frame(x = seq(0.5,1.2,0.05),
                            y = predict(fit, list(Ensemble_MSRG = seq(0.5,1.2,0.05))))
  q + geom_line(data = predicted_df, aes(x = x, y = y), colour = 'red')
  
}





{
  q <- ggplot(data = res, aes(x = meanFET, 
                              y = repair_probability)) + geom_point(size = 1.5) + theme_bw()
  
  
  #### fitting
  fit <- lm(data = res, formula = repair_probability ~ poly(meanFET))
  x = seq(0, 12, 0.1)
  predicted_df = data.frame(x = x,
                            y = predict(fit, list(meanFET = x)))
  q + geom_line(data = predicted_df, aes(x = x, y = y), colour = 'red')
  
}



##### Violin plots : STATISTICAL PROPERTY vs. ENCOUNTER EVENT

# 1. String lists to numeric vectors

library(gsubfn)
str2numvector <- function(strlist){
  z <- strapply(as.character(strlist), "(\\d*\\.*\\d+|\\d+\\.|\\d*\\.*\\d+e\\+\\d+|\\d*\\.*\\d+e\\-\\d+)", as.numeric)
  #z <- strapply(as.character(strlist), "(-?(\\d*\\.*\\d+|\\d+\\.))", as.numeric)
  #z <- lapply(as.character(strlist), function(x) gsub("[^[:digit:].e+]","",x))
  z <- cbind(lapply(z, 2, FUN = unlist))
  return(z)
#return(unlist(lapply(z, mean)))
}

msrgs <- str2numvector(res$MSRG_atEncounter)
msrgs
plot(res$Ensemble_MSRG, unlist(lapply(msrgs, mean)))
FETs <- str2numvector(res$FETs)

# 2.



f <- function(x, a, l){
  return(a*exp(-l*x))
}

{
  k <- ggplot(data = res, aes(x = ))
}

library(gsubfn)
msrgs <- res$MSRG_atEncounter
z <- strapply(as.character(msrgs), "(-?(\\d*\\.*\\d+|\\d+\\.))", as.numeric)
z <- cbind(lapply(z, 2, FUN = unlist))
means <- unlist(lapply(z, mean))
{k <- ggplot(data = res, aes(x = excludedVolumeCutOff, 
                             y = means, 
                             color = keepCL)) + geom_line(size = 1.5) + theme_bw()
k}
