res.k <- read.csv("../DNA_Religation/projects/RCLPolymer_CLs@DamageFoci/results/2018_06_27_17_38_proba-v_gNc_withVE.csv")
res.r <- read.csv("../DNA_Religation/projects/RCLPolymer_CLs@DamageFoci/results/2018_06_27_19_15_proba-v_gNc_withVE.csv")


library(ggplot2)
library(RColorBrewer)

res.k$excludedVolumeCutOff <- as.factor(res.k$excludedVolumeCutOff)

res.k$genomicDistance <- as.factor(res.k$genomicDistance)
{
  p <- ggplot(data = res.k, aes(x = Nc, y = repair_probability, color = genomicDistance)) + geom_line()
  p
}


res.r$excludedVolumeCutOff <- as.factor(res.r$excludedVolumeCutOff)

res.r$genomicDistance <- as.factor(res.r$genomicDistance)
{
  p <- ggplot(data = res.r, aes(x = Nc, y = repair_probability, color = genomicDistance)) + geom_line()
  p
}


#### Simulataneous plot

res.k$Scenario <- "Keeping CLs"
res.r$Scenario <- "Removing CLs"

res <- rbind(res.k, res.r)

{
p <- ggplot(data = res, mapping = aes(x = Nc, y = repair_probability, color = genomicDistance)) + 
  facet_grid(.~panel, scale="free") + 
  geom_line(data = res.k, stat = "identity") + 
  geom_line(data = res.r, stat = "identity")
p
}



{
  q <- ggplot(data = res, mapping = aes(x = Nc, y = Ensemble_MSRG, color = Scenario, shape = Scenario)) + 
    facet_grid(.~genomicDistance, scale="free") + 
    geom_line(data = res.k, stat = "identity", size = 2) + geom_point(data = res.k, stat = "identity", size = 4) +
    geom_line(data = res.r, stat = "identity", size = 2) + geom_point(data = res.r, stat = "identity", size = 4)
  q 
}

####
{
  w = c(rep(rep(1,5),3),rep(9,5),rep(rep(4,5),10))
  p <- ggplot(data = res, aes(x = genomicDistance, y = Nc, width = w)) + 
   geom_tile(aes(fill=repair_probability))
  
  p +  scale_fill_gradientn(colours= heat.colors(100))
}
####


##### VERIFICATION   ADAPTIVE-EPSILON VS MRG
res$genomicDistance <- as.factor(res$genomicDistance)
{
  p <- ggplot(data = res, aes(x = Nc, y = encounterDistance)) + geom_point(size = 2) + geom_line(size = 1)
  p + theme_bw() #+ geom_line(aes(x = Nc, y = Ensemble_MSRG, color = genomicDistance), size = 1) 
}

{
  p <- ggplot(data = res, aes(x = encounterDistance, y = sqrt(Ensemble_MSRG), color = genomicDistance)) + geom_point()
  p
}

relation <- lm(data = res, formula = excludedVolumeCutOff ~ sqrt(Ensemble_MSRG))
summary(relation)

##########################################

library(plotly)

#######
pp <- ggplotly(p)
ggplotly(p) %>% add_surface(z = pp$z)
######

########
{proba <- matrix(res$repair_probability, ncol = 14)
mfet <- matrix(res$meanFET, ncol = 14)
gs <- c(seq(2,4),seq(9,49,4))
nc <- seq(2,19,4)
py <- plot_ly(x = gs, y = nc, z = proba) %>% add_surface()
py}
#########


########
{
  proba <- matrix(res$repair_probability, ncol = 18)
  gs <- seq(2,19)
  nc <- seq(2,18,2)
  py <- plot_ly(data = res, x = gs, y = nc, z = proba) %>% add_surface()
  py
}



########


########
proba.REMOVING <- matrix(res$repair_probability, ncol = 2)
proba.KEEPING <- matrix(res$repair_probability, ncol = 2)
gs = c(4,20)
nc = seq(20,59,4)
py <- plot_ly(x = gs, y = nc, z = proba.REMOVING)
py
py <- plot_ly(x = gs, y = nc) %>% add_surface(z = proba.KEEPING, cmin = min(proba.KEEPING), cmax = max(proba.KEEPING), colorscale = list(c(0,1),c("rgb(255,112,184)","rgb(128,0,64)"))) %>% 
  add_surface(z = proba.REMOVING,cmin = min(proba.REMOVING), cmax = max(proba.REMOVING), colorscale = list(c(0,1),c("rgb(107,184,214)","rgb(0,90,124)")))
py
#########