res <- read.csv("../DNA_Religation/projects/RCLPolymer_CLs@DamageFoci/results/2018_05_14_14_37_proba_vs_genomicDistance_Ncs_.csv")

library(ggplot2)

{
p <- ggplot(data = res, aes(x = genomicDistance, y = repair_probability, color = keepCL)) + geom_line(size = 1.5)

p + geom_errorbar(data = res, aes(x = genomicDistance, 
                                  ymin = repair_probability - repair_CIhalflength,
                                  ymax = repair_probability + repair_CIhalflength), size = 1, width = 0.5) + theme_bw()

}

{p <- ggplot(data = res, aes(x = genomicDistance, 
                            y = meanFET, 
                            color = keepCL)) + geom_line(size = 1.5) + geom_errorbar(aes(x = genomicDistance, 
                                                                             ymin = meanFET - halfCI_FET,
                                                                             ymax = meanFET + halfCI_FET), size = 1, width = 1) + theme_bw()
p}


{q <- ggplot(data = res, aes(x = genomicDistance, 
                            y = Ensemble_MSRG, 
                            color = keepCL)) + geom_line(size = 1.5) + theme_bw()
q}

f <- function(x, a, l){
  return(a*exp(-l*x))
}

{
  h <- ggplot(data = res, aes(x = ))
}
