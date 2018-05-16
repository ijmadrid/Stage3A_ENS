res <- read.csv("./Stage 3A ENS/DNA_Religation/projects/RCLPolymer_CLs@DamageFoci/results/2018_05_16_16_09_proba-v_sigma.csv")

library(ggplot2)

res$excludedVolumeCutOff <- as.factor(res$excludedVolumeCutOff)

{
p <- ggplot(data = res, aes(x = excludedVolumeCutOff, y = repair_probability, color = keepCL)) + geom_line(size = 1.5) + geom_ribbon(data = res, aes(x = excludedVolumeCutOff, 
                                  ymin = repair_probability - repair_CIhalflength,
                                  ymax = repair_probability + repair_CIhalflength), alpha = 0.2) + theme_bw()
p
}

{q <- ggplot(data = res, aes(x = excludedVolumeCutOff, 
                            y = meanFET, 
                            color = keepCL)) + geom_line(size = 1.5) + geom_ribbon(aes(x = excludedVolumeCutOff, 
                                                                             ymin = meanFET - halfCI_FET,
                                                                             ymax = meanFET + halfCI_FET), alpha = 0.2) + theme_bw()
q}


{k <- ggplot(data = res, aes(x = excludedVolumeCutOff, 
                            y = Ensemble_MSRG, 
                            color = keepCL)) + geom_line(size = 1.5)  + theme_bw()
k}

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
