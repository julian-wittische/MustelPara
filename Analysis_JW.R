################################################################################
############### ANALYSIS OF MUSTELID PARASITE LENGTH AND NUMBER ################
################################################################################

# Author: Julian Wittische
# Data: November 2023

# /!\ WARNING /!\: when there is no random effects in the Gaussian model, performance::r2()
# does not work (called by tab_model()). Use glmmTMB::r.squaredGLMM() instead
# It is being fix by the performance team

library("readxl")
df <- read_excel("mesurements_S.nasicola_host-species_max.xlsx", "Sheet5")
colnames(df)[1] <- "length.of.worm"
df$sex.of.worm <- as.factor(df$sex.of.worm)
df$host.species <- as.factor(df$host.species)
df$sex.of.host <- as.factor(df$sex.of.host)
df$which.side <- as.factor(df$which.side)
df$id <- as.factor(df$id)
levels(df$host.species)[levels(df$host.species)=="M.vision"] <- "N.vison"
df$host.species
################################################################################
library(glmmTMB)
library(performance)
library(MuMIn)
library(sjPlot)
library(ggplot2)
library(ggeffects)

### Check collinearity

# Condylobasal length of host and host species are highly correlated
mod_nointer <- glmmTMB( length.of.worm ~ sex.of.worm + number.of.worms + host.species + sex.of.host + 
                          condylobasal.length.of.host.species + which.side + (1|id),
                        data = df, na.action = na.fail)
check_collinearity(mod_nointer)

mod_nointer <- glmmTMB( number.of.worms  ~ sex.of.worm + length.of.worm + host.species + sex.of.host + 
                          condylobasal.length.of.host.species + which.side + (1|id),
                        data = df, na.action = na.fail)
check_collinearity(mod_nointer)

# Length of worms ~ condylobasal length
mod_nointer <- glmmTMB( length.of.worm ~ sex.of.worm + number.of.worms + sex.of.host + 
                          condylobasal.length.of.host.species + which.side + (1|id),
                        data = df, na.action = na.fail)
check_collinearity(mod_nointer)

# Length of worms ~ host species
mod_nointer <- glmmTMB( length.of.worm ~ sex.of.worm + number.of.worms + host.species + sex.of.host + 
                          which.side + (1|id),
                        data = df, na.action = na.fail)
check_collinearity(mod_nointer)

# Number of worms ~ condylobasal length
mod_nointer <- glmmTMB( number.of.worms ~ sex.of.worm + length.of.worm + sex.of.host + 
                          condylobasal.length.of.host.species + which.side + (1|id),
                        data = df, na.action = na.fail,
                        family=poisson(link="log"))
check_collinearity(mod_nointer)

# Number of worms ~ host species
mod_nointer <- glmmTMB( number.of.worms ~ sex.of.worm + length.of.worm + host.species + sex.of.host + 
                          which.side + (1|id),
                        data = df, na.action = na.fail,
                        family=poisson(link="log"))
check_collinearity(mod_nointer)
#NOTE: as expected, host species and its condylobasal length is highly correlated
#ACTION: we have to get rid of one -> do with one then the other

################################################################################
################################################################################
######################## LENGTH OF WORM ########################################
################################################################################
################################################################################

################################################################################
######################## WITH CONDYLOBASAL LENGTH ##############################
################################################################################

### Let's put id as random effects and use REML (Zuur et al 2009)

mod_norand <- glmmTMB( length.of.worm ~ (sex.of.worm + number.of.worms + sex.of.host + 
                          condylobasal.length.of.host.species + which.side)^2,
                        data = df, na.action = na.fail, REML = TRUE)

mod_rand <- glmmTMB( length.of.worm ~ (sex.of.worm + number.of.worms + sex.of.host + 
                          condylobasal.length.of.host.species + which.side)^2 + (1|id),
                     data = df, na.action = na.fail, REML = TRUE)

summary(mod_norand)$AIC[1]
summary(mod_rand)$AIC[1]
#NOTE: apparently no need to include id even though many worms come from same individuals

### Full model summary with first level interactions

mod_full <- glmmTMB(length.of.worm ~ (sex.of.worm + number.of.worms + sex.of.host + 
                         condylobasal.length.of.host.species + which.side)^2,
                    data = df, na.action = na.fail, REML = FALSE)

summary(mod_full)
summary(mod_full)$AIC[1]

### Model selection using AIC (fixed effects)
### Use ML not REML before running model averaging (Ben Bolker on Cross-Validated)
### Using AICc is even more controversial than AIC for GLMM -> use AIC to rank
### abund_T_randland_dred <- MuMIn::dredge(abund_T_randland, rank = "AIC")

### Let's susbet to models with an AIC within 2 of the lowest.
dred <- dredge(mod_full, rank = "AIC")
dreddelta2 <- subset(dred, delta<2)
beepr::beep(3)
### Most parsimonious model
mod_parsi <- get.models(dreddelta2, 11)[[1]]
summary(mod_parsi)
tab_model(mod_parsi)

mean(df$length.of.worm[df$sex.of.worm=="F"])
mean(df$length.of.worm[df$sex.of.worm=="M"])

### The above results hint that there is not one "true" model, so let's try model averaging
mod_full_avg <- model.avg(dreddelta2, rank="AIC", fit=TRUE)
mod_full_avg_summary <- summary(mod_full_avg)
mod_full_avg_summary
mod_full_avg_summary$coefmat.full[,1]
mod_full_avg_summary$coef.nmod

### Plotting
best_plot <- plot_model(mod_parsi, type = "pred", terms = c("condylobasal.length.of.host.species","sex.of.worm"),
                               title = '',
                               axis.title = c("Condylobasal length of host", expression("Predicted length of worm")),
                               legend.title = "Sex of worm",
                               ci.lvl = 0.95, se=TRUE) + scale_x_continuous(limits = c(0, 800))

best_plot
jpeg("Figure1.jpg", width=14, height=8, units="cm", res=600)

best_plot
dev.off()

intercept_plot <- plot_model(mod_parsi, type = "pred", terms = c("condylobasal.length.of.host.species [0,800]","sex.of.worm"),
                        title = '',
                        axis.title = c("Condylobasal length of host", expression("Predicted length of worm")),
                        legend.title = "Sex of worm",
                        ci.lvl = 0.95, se=TRUE)
intercept_plot

################################################################################
######################## WITH HOST SPECIES #####################################
################################################################################

### Let's put id as random effects and use REML (Zuur et al 2009)

mod_norandHS <- glmmTMB( length.of.worm ~ (sex.of.worm + number.of.worms + sex.of.host + 
                                          host.species + which.side)^2,
                       data = df, na.action = na.fail, REML = TRUE)

mod_randHS <- glmmTMB( length.of.worm ~ (sex.of.worm + number.of.worms + sex.of.host + 
                                        host.species + which.side)^2 + (1|id),
                     data = df, na.action = na.fail, REML = TRUE)

summary(mod_norandHS)$AIC[1]
summary(mod_randHS)$AIC[1]
#NOTE: apparently no need to include id even though many worms come from same individuals

### Full model summary with first level interactions

mod_fullHS <- glmmTMB(length.of.worm ~ (sex.of.worm + number.of.worms + sex.of.host + 
                                       host.species + which.side)^2,
                    data = df, na.action = na.fail, REML = FALSE)

summary(mod_fullHS)
summary(mod_fullHS)$AIC[1]

### Model selection using AIC (fixed effects)
### Use ML not REML before running model averaging (Ben Bolker on Cross-Validated)
### Using AICc is even more controversial than AIC for GLMM -> use AIC to rank
### abund_T_randland_dred <- MuMIn::dredge(abund_T_randland, rank = "AIC")

### Let's susbet to models with an AIC within 2 of the lowest.
dredHS <- dredge(mod_fullHS, rank = "AIC")
dreddelta2HS <- subset(dredHS, delta<2)
beepr::beep(3)
### Most parsimonious model
mod_parsiHS <- get.models(dreddelta2HS, 1)[[1]]
summary(mod_parsiHS)
tab_model(mod_parsiHS)

mean(df$length.of.worm[df$sex.of.worm=="F"])
mean(df$length.of.worm[df$sex.of.worm=="M"])

### The above results hint that there is not one "true" model, so let's try model averaging
mod_full_avgHS <- model.avg(dreddelta2HS, rank="AIC", fit=TRUE)
mod_full_avg_summaryHS <- summary(mod_full_avgHS)
mod_full_avg_summaryHS
mod_full_avg_summaryHS$coefmat.full[,1]
mod_full_avg_summaryHS$coef.nmod

### Plotting
best_plot <- plot_model(mod_parsiHS, type = "pred", terms = c("host.species","sex.of.worm"),
                        title = '',
                        axis.title = c("Host species", expression("Predicted length of worm")),
                        legend.title = "Sex of worm",
                        ci.lvl = 0.95, se=TRUE)
jpeg("Figure2.jpg", width=14, height=8, units="cm", res=600)
best_plot
dev.off()

################################################################################
################################################################################
######################## NUMBERS OF WORMS ######################################
################################################################################
################################################################################

### Check collinearity

mod_nointer <- glmmTMB( number.of.worms ~ sex.of.worm + length.of.worm + host.species + sex.of.host + 
                          condylobasal.length.of.host.species + which.side + (1|id),
                        data = df, na.action = na.fail,
                        family=poisson(link="log"))
check_collinearity(mod_nointer)
#NOTE: as expected, host species and its condylobasal length is highly correlated
#ACTION: we have to get rid of one -> do with one then the other

################################################################################
######################## WITH CONDYLOBASAL LENGTH ##############################
################################################################################

### Let's put id as random effects and use REML (Zuur et al 2009)

mod_norandNUM <- glmmTMB( number.of.worms ~ (sex.of.worm + length.of.worm + sex.of.host + 
                                          condylobasal.length.of.host.species + which.side)^2,
                       data = df, na.action = na.fail, REML = TRUE,
                       family=poisson(link="log"))

mod_randNUM <- glmmTMB( number.of.worms ~ (sex.of.worm + length.of.worm + sex.of.host + 
                                          condylobasal.length.of.host.species + which.side)^2 + (1|id),
                     data = df, na.action = na.fail, REML = TRUE,
                     family=poisson(link="log"))

summary(mod_norandNUM)$AIC[1]
summary(mod_randNUM)$AIC[1]
#NOTE: apparently we need to include id

### Full model summary with first level interactions

mod_fullNUM <- glmmTMB(number.of.worms ~ (sex.of.worm + length.of.worm + sex.of.host + 
                                         condylobasal.length.of.host.species + which.side)^2+ (1|id),
                    data = df, na.action = na.fail, REML = FALSE,
                    family=poisson(link="log"))

summary(mod_fullNUM)
summary(mod_fullNUM)$AIC[1]
check_overdispersion(mod_fullNUM)

### Model selection using AIC (fixed effects)
### Use ML not REML before running model averaging (Ben Bolker on Cross-Validated)
### Using AICc is even more controversial than AIC for GLMM -> use AIC to rank
### abund_T_randland_dred <- MuMIn::dredge(abund_T_randland, rank = "AIC")

### Let's susbet to models with an AIC within 2 of the lowest.
dredNUM <- dredge(mod_fullNUM, rank = "AIC")
dreddelta2NUM <- subset(dredNUM, delta<2)
beepr::beep(3)
### Most parsimonious model
mod_parsiNUM <- get.models(dreddelta2NUM, 1)[[1]]
summary(mod_parsiNUM)
tab_model(mod_parsiNUM)

mean(df$number.of.worms[df$sex.of.worm=="F"])
mean(df$number.of.worms[df$sex.of.worm=="M"])

### The above results hint that there is not one "true" model, so let's try model averaging
mod_full_avgNUM <- model.avg(dreddelta2NUM, rank="AIC", fit=TRUE)
mod_full_avg_summaryNUM <- summary(mod_full_avgNUM)
mod_full_avg_summaryNUM
mod_full_avg_summaryNUM$coefmat.full[,1]
mod_full_avg_summaryNUM$coef.nmod

### Plotting
best_plotNUM <- plot_model(mod_parsiNUM, type = "pred", terms = c("sex.of.host","sex.of.worm", "which.side"),
                        title = '',
                        axis.title = c("Sex of host", expression("Predicted number of worms")),
                        legend.title = "Sex of worm",
                        ci.lvl = 0.95, se=TRUE)
best_plotNUM
jpeg("Figure3.jpg", width=14, height=8, units="cm", res=600)
best_plotNUM
dev.off()

################################################################################
######################## WITH HOST SPECIES #####################################
################################################################################

### Let's put id as random effects and use REML (Zuur et al 2009)

mod_norandNUMHS <- glmmTMB( number.of.worms ~ (sex.of.worm + length.of.worm + sex.of.host + 
                                               host.species + which.side)^2,
                          data = df, na.action = na.fail, REML = TRUE,
                          family=poisson(link="log"))

mod_randNUMHS <- glmmTMB( number.of.worms ~ (sex.of.worm + length.of.worm + sex.of.host + 
                                             host.species + which.side)^2 + (1|id),
                        data = df, na.action = na.fail, REML = TRUE,
                        family=poisson(link="log"))

summary(mod_norandNUMHS)$AIC[1]
summary(mod_randNUMHS)$AIC[1]
#NOTE: apparently we need to include id

### Full model summary with first level interactions

mod_fullNUMHS <- glmmTMB(number.of.worms ~ (sex.of.worm + length.of.worm + sex.of.host + 
                                           host.species + which.side)^2+ (1|id),
                       data = df, na.action = na.fail, REML = FALSE,
                       family=poisson(link="log"))

summary(mod_fullNUMHS)
summary(mod_fullNUMHS)$AIC[1]

### Model selection using AIC (fixed effects)
### Use ML not REML before running model averaging (Ben Bolker on Cross-Validated)
### Using AICc is even more controversial than AIC for GLMM -> use AIC to rank
### abund_T_randland_dred <- MuMIn::dredge(abund_T_randland, rank = "AIC")

### Let's susbet to models with an AIC within 2 of the lowest.
dredNUMHS <- dredge(mod_fullNUMHS, rank = "AIC")
dreddelta2NUMHS <- subset(dredNUMHS, delta<2)
beepr::beep(3)
### Most parsimonious model
mod_parsiNUMHS <- get.models(dreddelta2NUMHS, 1)[[1]]
summary(mod_parsiNUMHS)
tab_model(mod_parsiNUMHS)

mean(df$number.of.worms[df$sex.of.worm=="F"])
mean(df$number.of.worms[df$sex.of.worm=="M"])

### The above results hint that there is not one "true" model, so let's try model averaging
mod_full_avgNUMHS <- model.avg(dreddelta2NUMHS, rank="AIC", fit=TRUE)
mod_full_avg_summaryNUMHS <- summary(mod_full_avgNUMHS)
mod_full_avg_summaryNUMHS
mod_full_avg_summaryNUMHS$coefmat.full[,1]
mod_full_avg_summaryNUMHS$coef.nmod

### Plotting
best_plotNUMHS <- plot_model(mod_parsiNUMHS, type = "pred", terms = c("host.species","sex.of.host", "which.side"),
                           title = '',
                           axis.title = c("Host species", expression("Predicted number of worms")),
                           legend.title = "Sex of host",
                           ci.lvl = 0.95, se=TRUE)
best_plotNUMHS
jpeg("Figure4.jpg", width=14, height=8, units="cm", res=600)
best_plotNUMHS
dev.off()


################################################################################
######################## WITH HOST SPECIES #####################################
################################################################################

mod_ermina <- glmmTMB(length.of.worm ~ (sex.of.worm + condylobasal.length.of.host.species)^2,
                    data = df[df$host.species=="M.ermina",], na.action = na.fail, REML = FALSE)

tab_model(mod_ermina)

plot_ermina <- plot_model(mod_ermina, type = "pred", terms = c("condylobasal.length.of.host.species","sex.of.worm"),
                          title = '',
                          axis.title = c("Condylobasal length of host", expression("Predicted length of worms")),
                          legend.title = "",
                          ci.lvl = 0.95, se=TRUE)
plot_ermina

mod_nivalis <- glmmTMB(length.of.worm ~ (sex.of.worm + condylobasal.length.of.host.species)^2,
                      data = df[df$host.species=="M.nivalis",], na.action = na.fail, REML = FALSE)

tab_model(mod_nivalis)

plot_nivalis <- plot_model(mod_nivalis, type = "pred", terms = c("condylobasal.length.of.host.species","sex.of.worm"),
                           title = '',
                           axis.title = c("Condylobasal length of host", expression("Predicted length of worms")),
                           legend.title = "",
                           ci.lvl = 0.95, se=TRUE)
plot_nivalis


mod_putorius <- glmmTMB(length.of.worm ~ (sex.of.worm + condylobasal.length.of.host.species)^2,
                      data = df[df$host.species=="M.putorius",], na.action = na.fail, REML = FALSE)

tab_model(mod_putorius)

plot_putorius <- plot_model(mod_putorius, type = "pred", terms = c("condylobasal.length.of.host.species","sex.of.worm"),
                           title = '',
                           axis.title = c("Condylobasal length of host", expression("Predicted length of worms")),
                           legend.title = "",
                           ci.lvl = 0.95, se=TRUE)
plot_putorius

mod_vison <- glmmTMB(length.of.worm ~ (sex.of.worm + condylobasal.length.of.host.species)^2,
                      data = df[df$host.species=="N.vison",], na.action = na.fail, REML = FALSE)

tab_model(mod_vison)

plot_vison <- plot_model(mod_vison, type = "pred", terms = c("condylobasal.length.of.host.species","sex.of.worm"),
                            title = '',
                            axis.title = c("Condylobasal length of host", expression("Predicted length of worms")),
                            legend.title = "",
                            ci.lvl = 0.95, se=TRUE)
plot_vison

################################################################################
################### REDO MODELS WITHIN SPECIES #################################
################################################################################

### M. ermina
ermina <- df[df$host.species=="M.ermina",]

mod_norand_ermina <- glmmTMB( length.of.worm ~ (sex.of.worm + number.of.worms + sex.of.host + 
                                           condylobasal.length.of.host.species + which.side)^2,
                       data = ermina, na.action = na.fail, REML = TRUE)

mod_rand_ermina <- glmmTMB( length.of.worm ~ (sex.of.worm + number.of.worms + sex.of.host + 
                                         condylobasal.length.of.host.species + which.side)^2 + (1|id),
                     data = ermina, na.action = na.fail, REML = TRUE)

summary(mod_norand_ermina)$AIC[1]
summary(mod_rand_ermina)$AIC[1]

# No need to include random effects

mod_ermina_full <- glmmTMB(length.of.worm ~ (sex.of.worm + number.of.worms + sex.of.host + 
                                          condylobasal.length.of.host.species + which.side)^2,
                      data = ermina, na.action = na.fail, REML = FALSE)

### Let's susbet to models with an AIC within 2 of the lowest.
dredermina <- dredge(mod_ermina_full, rank = "AIC")
dreddelta2ermina <- subset(dredermina, delta<2)
beepr::beep(3)
### Most parsimonious model
mod_parsi_ermina <- get.models(dreddelta2ermina, 8)[[1]]
summary(mod_parsi_ermina)
tab_model(mod_parsi_ermina)

### M. nivalis
nivalis <- df[df$host.species=="M.nivalis",]

mod_norand_nivalis <- glmmTMB( length.of.worm ~ (sex.of.worm + number.of.worms + sex.of.host + 
                                                  condylobasal.length.of.host.species + which.side)^2,
                              data = nivalis, na.action = na.fail, REML = TRUE)

mod_rand_nivalis <- glmmTMB( length.of.worm ~ (sex.of.worm + number.of.worms + sex.of.host + 
                                                condylobasal.length.of.host.species + which.side)^2 + (1|id),
                            data = nivalis, na.action = na.fail, REML = TRUE)

summary(mod_norand_nivalis)$AIC[1]
summary(mod_rand_nivalis)$AIC[1]

# No need to include random effects

mod_nivalis_full <- glmmTMB(length.of.worm ~ (sex.of.worm + number.of.worms + sex.of.host + 
                                               condylobasal.length.of.host.species + which.side)^2,
                           data = nivalis, na.action = na.fail, REML = FALSE)

### Let's susbet to models with an AIC within 2 of the lowest.
drednivalis <- dredge(mod_nivalis_full, rank = "AIC")
dreddelta2nivalis <- subset(drednivalis, delta<2)
beepr::beep(3)

### Most parsimonious model (THERE ARE TO MODELS WITH THE SAME DEGREES OF FREEDOM)
mod_parsi_nivalis <- get.models(dreddelta2nivalis, 13)[[1]]
summary(mod_parsi_nivalis)
tab_model(mod_parsi_nivalis)

mod_parsi_nivalis2 <- get.models(dreddelta2nivalis, 35)[[1]]
summary(mod_parsi_nivalis2)
tab_model(mod_parsi_nivalis2)

### M. putorius
putorius <- df[df$host.species=="M.putorius",]

mod_norand_putorius <- glmmTMB( length.of.worm ~ (sex.of.worm + number.of.worms + sex.of.host + 
                                                  condylobasal.length.of.host.species + which.side)^2,
                              data = putorius, na.action = na.fail, REML = TRUE)

mod_rand_putorius <- glmmTMB( length.of.worm ~ (sex.of.worm + number.of.worms + sex.of.host + 
                                                condylobasal.length.of.host.species + which.side)^2 + (1|id),
                            data = putorius, na.action = na.fail, REML = TRUE)

summary(mod_norand_putorius)$AIC[1]
summary(mod_rand_putorius)$AIC[1]

# No need to include random effects

mod_putorius_full <- glmmTMB(length.of.worm ~ (sex.of.worm + number.of.worms + sex.of.host + 
                                               condylobasal.length.of.host.species + which.side)^2,
                           data = putorius, na.action = na.fail, REML = FALSE)

### Let's susbet to models with an AIC within 2 of the lowest.
dredputorius <- dredge(mod_putorius_full, rank = "AIC")
dreddelta2putorius <- subset(dredputorius, delta<2)
beepr::beep(3)
### Most parsimonious model
mod_parsi_putorius <- get.models(dreddelta2putorius, 1)[[1]]
summary(mod_parsi_putorius)
tab_model(mod_parsi_putorius)    

### N. vison
vison <- df[df$host.species=="N.vison",]

mod_norand_vison <- glmmTMB( length.of.worm ~ (sex.of.worm + number.of.worms + sex.of.host + 
                                                    condylobasal.length.of.host.species + which.side)^2,
                                data = vison, na.action = na.fail, REML = TRUE)

mod_rand_vison <- glmmTMB( length.of.worm ~ (sex.of.worm + number.of.worms + sex.of.host + 
                                                  condylobasal.length.of.host.species + which.side)^2 + (1|id),
                              data = vison, na.action = na.fail, REML = TRUE)

summary(mod_norand_vison)$AIC[1]
summary(mod_rand_vison)$AIC[1]

# No need to include random effects

mod_vison_full <- glmmTMB(length.of.worm ~ (sex.of.worm + number.of.worms + sex.of.host + 
                                                 condylobasal.length.of.host.species + which.side)^2,
                             data = vison, na.action = na.fail, REML = FALSE)

### Let's susbet to models with an AIC within 2 of the lowest.
dredvison <- dredge(mod_vison_full, rank = "AIC")
dreddelta2vison <- subset(dredvison, delta<2)
beepr::beep(3)
### Most parsimonious model
mod_parsi_vison <- get.models(dreddelta2vison, 1)[[1]]
summary(mod_parsi_vison)
tab_model(mod_parsi_vison)
beepr::beep(1)


###############################################################################################################
###############################################################################################################
### CODE FOR ALAIN: FOR LOOP TO EXTRACT R2
R2df <- data.frame(marginal=numeric(), conditional=numeric())
for (i in 1: nrow(dreddelta2)){
  R2df[i,] <- r.squaredGLMM(get.models(dreddelta2, i)[[1]])
}

### CODE FOR ALAIN: plotting model averaging results
# Coefficients plots
plot_models(mod_full_avg)

library(parameters)
library(see)

mod_full_avg <- model.avg(dreddelta2, rank="AIC", fit=TRUE)
mod_full_avg_summary
mod <- model_parameters(mod_full_avg)
plot(mod)

# NOT TRIED:
# library(papeR)
# prettify()

# Prediction plots

plot_model(mod_full_avg_summary, type = "pred", terms = c("condylobasal.length.of.host.species","sex.of.worm"),
           title = '',
           axis.title = c("Condylobasal length of host", expression("Predicted length of worms")),
           legend.title = "",
           ci.lvl = 0.95, se=TRUE) #No CI, known issue

plot(ggpredict(mod_full_avg, terms = c("condylobasal.length.of.host.species","sex.of.worm"))

# NOT TRIED:
# library(tidyverse)
# library(AICcmodavg)
# # create new data with only CBL varying?
# modavg <- modavgPred(dredelta2, newdata = .)
# ggplot(modavg,aes(x=CBL, y = mod.avg.pred, group = type,
#            color = type, ymin = lower.CL, ymax=upper.CL,
#            fill = type)) +
# geom_line()
#ggeffects


