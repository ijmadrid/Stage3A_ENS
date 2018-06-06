res <- read.csv("./Stage 3A ENS/DNA_Religation/projects/RCLPolymer_TimesBelowThreshold/results/2018_06_06_16_38_trackDSB.csv")


library(ggplot2)

{p <- ggplot(data = res, aes(x = genomicDistance, y = Nc, fill = a1.b2_ComebackRate)) + geom_tile()
p}

res$genomicDistance <- as.factor(res$genomicDistance)
ggplot(data = res, aes(x = Nc, y = a1.a2_ComebackRate, color = genomicDistance)) + geom_line(size=2)




ggplot(data = res, aes(x = Nc, y = a1.a2_meanComebackTime, color = genomicDistance)) + geom_line() +
  geom_ribbon(aes(x = Nc, 
                  ymin = a1.a2_meanComebackTime - a1.a2_meanComebackTime_dx,
                  ymax = a1.a2_meanComebackTime + a1.a2_meanComebackTime_dx), 
              alpha = 0.2)

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

ggplot(data = res, aes(x = Nc, y = EstimatedProba, color = genomicDistance)) + geom_line()
