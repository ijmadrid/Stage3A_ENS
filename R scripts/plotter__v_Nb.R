res <- read.csv("./Stage 3A ENS/DNA_Religation/projects/RCLPolymer_RepairProbaVsNumberOfBreaks/results/2018_06_01_17_19_StatsVsNb.csv")

library(ggplot2)


res$genomicDistance <- as.factor(res$genomicDistance)

{
  p <- ggplot(data = res, aes(x = Nb, y = repair_probability,color = genomicDistance)) + geom_line(size = 1.5) + geom_ribbon(data = res, aes(x = Nb, 
                                                                                                                     ymin = repair_probability - repair_CIhalflength,
                                                                                                                     ymax = repair_probability + repair_CIhalflength), alpha = 0.2) + theme_bw()
  p
}

mixed <- data.frame(Nb = seq(2,9))
mixed$proba <- 1/(2*mixed$Nb - 1)
q <- ggplot(data = mixed, aes(x = Nb, y = proba)) + geom_line(color = 'red')
p + geom_line(aes(x = Nb, y = 1/(2*Nb - 1)), size = 1.5, color = 'red')


{
  p <- ggplot(data = res, aes(x = Nb, y = meanFET,color = genomicDistance)) + geom_line(size = 1.5) + geom_ribbon(data = res, aes(x = Nb, 
                                                                                                                                  ymin = meanFET - halfCI_FET,
                                                                                                                                  ymax = meanFET + halfCI_FET), alpha = 0.2) + theme_bw()
  p
}
