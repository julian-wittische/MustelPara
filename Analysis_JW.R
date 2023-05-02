setwd("T:/Julian/")

library("readxl")
library(MCMCglmm)


df <- read_excel("mesurements_S.nasicola_host-species_max.xlsx", "Sheet5")
df$sex.of.worm <- as.factor(df$sex.of.worm)
df$host.species <- as.factor(df$host.species)
df$sex.of.host <- as.factor(df$sex.of.host)
df$which.side <- as.factor(df$which.side)
df$id <- as.factor(df$id)


m.full <- MCMCglmm(lenth.of.worm ~ sex.of.worm*number.of.worms + host.species * sex.of.host + 
                     condylobasal.length.of.host.species + which.side
                     , random = ~id, data = df, family="gaussian")

summary(m.full)

par(mfrow=c(8,2), mar=c(2,2,1,0))
plot(m.full$Sol, auto.layout=F)


