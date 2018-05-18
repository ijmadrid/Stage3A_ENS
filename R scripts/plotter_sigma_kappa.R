res <- read.csv("./Stage 3A ENS/DNA_Religation/projects/RCLPolymer_CLs@DamageFoci/results/2018_05_17_18_10_proba-v_kappa.csv")
  
library(ggplot2)

res$excludedVolumeCutOff <- as.factor(res$excludedVolumeCutOff)

{
p <- ggplot(data = res, aes(x = excludedVolumeSpringConstant, y = repair_probability, color = excludedVolumeCutOff)) + geom_line(size = 1.5) + geom_ribbon(data = res, aes(x = excludedVolumeSpringConstant, 
                                  ymin = repair_probability - repair_CIhalflength,
                                  ymax = repair_probability + repair_CIhalflength), alpha = 0.2) + theme_bw()
p
}



####
{
  p <- ggplot(data = res, aes(x = excludedVolumeCutOff, y = excludedVolumeSpringConstant)) + geom_raster(aes(fill = repair_probability), interpolate = T)
  
  p <- p + scale_fill_gradientn(colours= heat.colors(100)) +
    geom_hline(yintercept = 0.6,linetype = 2,size=1) +
    geom_vline(xintercept = 0.05,linetype = 2,size=1) +
    geom_vline(xintercept = 0.2,linetype = 2,size=1) 
  
  p
}
####


{q <- ggplot(data = res, aes(x = excludedVolumeSpringConstant, 
                             y = meanFET, 
                             color = excludedVolumeCutOff)) + 
    geom_line(size = 1.5) + 
    geom_ribbon(aes(x = excludedVolumeSpringConstant, 
                    ymin = meanFET - halfCI_FET,
                    ymax = meanFET + halfCI_FET), 
                alpha = 0.2) + 
    theme_bw()
q}

####
{
  q2 <- ggplot(data = res, aes(x = excludedVolumeCutOff, y = excludedVolumeSpringConstant)) + geom_raster(aes(fill = meanFET), interpolate = T)
  
  q2 <- q2 + scale_fill_gradientn(colours= heat.colors(100)) +
    geom_hline(yintercept = 0.6,linetype = 2,size=1) +
    geom_vline(xintercept = 0.05,linetype = 2,size=1) +
    geom_vline(xintercept = 0.2,linetype = 2,size=1)
  
  q2
}
####


{k <- ggplot(data = res, aes(x = excludedVolumeSpringConstant, 
                             y = Ensemble_MSRG, 
                             color = excludedVolumeCutOff)) + 
    geom_line(size = 1.5) + 
    geom_ribbon(aes(x = excludedVolumeSpringConstant, 
                    ymin = Ensemble_MSRG - Ensemble_MSRG_dx,
                    ymax = Ensemble_MSRG + Ensemble_MSRG_dx), 
                alpha = 0.2) + 
    theme_bw()
k}

####
{
  k2 <- ggplot(data = res, aes(x = excludedVolumeCutOff, y = excludedVolumeSpringConstant)) + geom_raster(aes(fill = Ensemble_MSRG), interpolate = T)
  
  k2 <- k2 + scale_fill_gradientn(colours= heat.colors(100)) +
    geom_hline(yintercept = 0.6,linetype = 2,size=1) +
    geom_vline(xintercept = 0.05,linetype = 2,size=1) +
    geom_vline(xintercept = 0.2,linetype = 2,size=1)
  
  k2
}
####








#### Curve fitting
exp.fit <- lm(data = res, formula = log(meanFET) ~ excludedVolumeSpringConstant)
x <- seq(0,2,0.05)
fit <- exp(predict(exp.fit,list(excludedVolumeSpringConstant=x)))
predicted_df <- data.frame(fit = fit, x=x)
q + geom_line(data = predicted_df, aes(x = x, y = fit), colour = 'red')

{k <- ggplot(data = res, aes(x = excludedVolumeSpringConstant, 
                            y = Ensemble_MSRG)) + geom_line(size = 1.5)  + geom_ribbon(aes(x = excludedVolumeSpringConstant, 
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
  z <- strapply(as.character(strlist), "(\\d*\\.*\\d+|\\d+\\.|\\d*\\.*\\d+e\\+\\d+|\\d*\\.*\\d+e\\-\\d+|nan)", as.numeric)
  z <- cbind(sapply(z, 2, FUN = unlist))
  return(z)
}

msrgs <- str2numvector(res$MSRG_atEncounter)
msrgs[is.nan(msrgs)] <- NA
plot(res$Ensemble_MSRG, colMeans(msrgs,na.rm = T))
FETs <- str2numvector(res$FETs)
FETs[is.nan(FETs)] <- NA

str2eventlist <- function(strlist){
  z <- strapply(as.character(strlist), "(\\w+|\\w*\\-\\w+)")
  z <- cbind(sapply(z, 2, FUN = unlist))
  return(z)
}
events <- str2eventlist(res$events)
events[events == "NA"] <- NA

df <- data.frame(events = events, fets = FETs, msrgs = msrgs)

# 2. Violin plots
library(manipulate)

manipulate(
  ggplot(data = df <- subset(data.frame(events = events[,column], 
                                        fets = FETs[,column], 
                                        msrgs = msrgs[,column]), 
                             !is.na(events)), 
         aes(x = events, y = df[,factor], fill = events, color = events)) +
    geom_violin() + geom_boxplot(width=0.1, fill = "white",color="black") +
    coord_flip() + ylab(factor) +
    scale_fill_manual(values = c("red","coral","orange","yellow","green2")) + 
    scale_color_manual(values = c("red","coral","orange","yellow","green2")) +
    theme_bw() + theme(legend.position="none"),
  column = slider(1,30),
  factor = picker('fets','msrgs')
)

manipulate(
ggplot(data = df, aes(x = events.10, y = df[,factor], fill = events.10)) + geom_violin() + geom_boxplot(width=0.1, fill = "white") + theme(legend.position="none") + coord_flip(),
factor = picker("msrgs.10","fets.10")
)



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
