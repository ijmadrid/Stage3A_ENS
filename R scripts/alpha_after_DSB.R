res_01keep <- read.csv("./Stage 3A ENS/DNA_Religation/projects/AnomalousExponentAfterDSBs/results/2018_06_22_17_39_trackAlphaAfterDSB_vsNc.csv")
res_01remo <- read.csv("./Stage 3A ENS/DNA_Religation/projects/AnomalousExponentAfterDSBs/results/2018_06_27_13_08_trackAlphaAfterDSB_vsNc.csv")


get_alpha <- function(res){
alphas <- res[, grepl("alpha.", names(res))]
alphas$nc <- res$Nc
alphas$alpha_monomers <- NULL
return(alphas)
}

alphas.k <- get_alpha(res_01keep)
alphas.r <- get_alpha(res_01remo)
  
library(ggplot2)
library(reshape2)

#alphas.melted <- melt(alphas, id.vars = "nc")

alphas.k.melted <- melt(alphas.k, id.vars = "nc")
alphas.r.melted <- melt(alphas.r, id.vars = "nc")

alphas.avg <- data.frame(n = seq(100),
                         keep = sapply(alphas.k[,!names(alphas.k) %in% ("nc")], mean),
                         remo = sapply(alphas.r[,!names(alphas.r) %in% ("nc")], mean))
 
ggplot(data = alphas.r.melted, aes(x = variable, y = value, color = nc)) + 
  geom_point() + 
  geom_vline(xintercept = c(46,47,54,55)) +
  scale_color_gradient(low = "yellow", high = "red")


alphas.avg.melted <- melt(alphas.avg, id.vars = c("n"))

ggplot(data = alphas.avg.melted, aes(x= n, y = value, color = variable)) + 
  geom_point() +
  geom_vline(xintercept = c(46,47,54,55))

library(manipulate)

manipulate(
  {ggplot(data = alphas.melted[alphas.melted$nc == xnc,], aes(x = variable, y = value)) + 
      geom_point() + 
      geom_vline(xintercept = c(31,32,69,70))},
  xnc = slider(2,48,step=2)
)


ggplot(data = res, aes(x = Nc, y = alpha.54)) + geom_line()
