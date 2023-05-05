################################################################################
######
################################################################################

library("readxl")
df <- read_excel("mesurements_S.nasicola_host-species_max.xlsx", "Sheet5")
df$sex.of.worm <- as.factor(df$sex.of.worm)
df$host.species <- as.factor(df$host.species)
df$sex.of.host <- as.factor(df$sex.of.host)
df$which.side <- as.factor(df$which.side)
df$id <- as.factor(df$id)

################################################################################
library(glmmTMB)
library(performance)
library(MuMIn)
library(sjPlot)
library(ggplot2)
library(ggeffects)

################################################################################
################################################################################
######################## LENGTH OF WORM ########################################
################################################################################
################################################################################

### Check collinearity

mod_nointer <- glmmTMB( lenth.of.worm ~ sex.of.worm + number.of.worms + host.species + sex.of.host + 
                   condylobasal.length.of.host.species + which.side + (1|id),
                   data = df, na.action = na.fail)
check_collinearity(mod_nointer)
#NOTE: as expected, host species and its condylobasal length is highly correlated
#ACTION: we have to get rid of one -> do with one then the other

################################################################################
######################## WITH CONDYLOBASAL LENGTH ##############################
################################################################################

### Let's put id as random effects and use REML (Zuur et al 2009)

mod_norand <- glmmTMB( lenth.of.worm ~ (sex.of.worm + number.of.worms + sex.of.host + 
                          condylobasal.length.of.host.species + which.side)^2,
                        data = df, na.action = na.fail, REML = TRUE)

mod_rand <- glmmTMB( lenth.of.worm ~ (sex.of.worm + number.of.worms + sex.of.host + 
                          condylobasal.length.of.host.species + which.side)^2 + (1|id),
                     data = df, na.action = na.fail, REML = TRUE)

summary(mod_norand)$AIC[1]
summary(mod_rand)$AIC[1]
#NOTE: apparently no need to include id even though many worms come from same individuals

### Full model summary with first level interactions

mod_full <- glmmTMB(lenth.of.worm ~ (sex.of.worm + number.of.worms + sex.of.host + 
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

### Most parsimonious model
mod_parsi <- get.models(dreddelta2, 11)[[1]]
summary(mod_parsi)
tab_model(mod_parsi)

mean(df$lenth.of.worm[df$sex.of.worm=="F"])
mean(df$lenth.of.worm[df$sex.of.worm=="M"])

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
                               legend.title = "",
                               ci.lvl = 0.95, se=TRUE)
best_plot

################################################################################
######################## WITH HOST SPECIES #####################################
################################################################################

### Let's put id as random effects and use REML (Zuur et al 2009)

mod_norandHS <- glmmTMB( lenth.of.worm ~ (sex.of.worm + number.of.worms + sex.of.host + 
                                          host.species + which.side)^2,
                       data = df, na.action = na.fail, REML = TRUE)

mod_randHS <- glmmTMB( lenth.of.worm ~ (sex.of.worm + number.of.worms + sex.of.host + 
                                        host.species + which.side)^2 + (1|id),
                     data = df, na.action = na.fail, REML = TRUE)

summary(mod_norandHS)$AIC[1]
summary(mod_randHS)$AIC[1]
#NOTE: apparently no need to include id even though many worms come from same individuals

### Full model summary with first level interactions

mod_fullHS <- glmmTMB(lenth.of.worm ~ (sex.of.worm + number.of.worms + sex.of.host + 
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

### Most parsimonious model
mod_parsiHS <- get.models(dreddelta2HS, 1)[[1]]
summary(mod_parsiHS)
tab_model(mod_parsiHS)

mean(df$lenth.of.worm[df$sex.of.worm=="F"])
mean(df$lenth.of.worm[df$sex.of.worm=="M"])

### The above results hint that there is not one "true" model, so let's try model averaging
mod_full_avgHS <- model.avg(dreddelta2HS, rank="AIC", fit=TRUE)
mod_full_avg_summaryHS <- summary(mod_full_avgHS)
mod_full_avg_summaryHS
mod_full_avg_summary$coefmat.full[,1]
mod_full_avg_summary$coef.nmod

### Plotting
best_plot <- plot_model(mod_parsiHS, type = "pred", terms = c("host.species","sex.of.worm"),
                        title = '',
                        axis.title = c("Host species", expression("Predicted length of worm")),
                        legend.title = "",
                        ci.lvl = 0.95, se=TRUE)
best_plot

################################################################################
################################################################################
######################## NUMBERS OF WORMS ######################################
################################################################################
################################################################################

### Check collinearity

mod_nointer <- glmmTMB( number.of.worms ~ sex.of.worm + lenth.of.worm + host.species + sex.of.host + 
                          condylobasal.length.of.host.species + which.side + (1|id),
                        data = df, na.action = na.fail)
check_collinearity(mod_nointer)
#NOTE: as expected, host species and its condylobasal length is highly correlated
#ACTION: we have to get rid of one -> do with one then the other

################################################################################
######################## WITH CONDYLOBASAL LENGTH ##############################
################################################################################

### Let's put id as random effects and use REML (Zuur et al 2009)

mod_norandNUM <- glmmTMB( number.of.worms ~ (sex.of.worm + lenth.of.worm + sex.of.host + 
                                          condylobasal.length.of.host.species + which.side)^2,
                       data = df, na.action = na.fail, REML = TRUE)

mod_randNUM <- glmmTMB( number.of.worms ~ (sex.of.worm + lenth.of.worm + sex.of.host + 
                                          condylobasal.length.of.host.species + which.side)^2 + (1|id),
                     data = df, na.action = na.fail, REML = TRUE)

summary(mod_norandNUM)$AIC[1]
summary(mod_randNUM)$AIC[1]
#NOTE: apparently we need to include id

### Full model summary with first level interactions

mod_fullNUMHS <- glmmTMB(number.of.worms ~ (sex.of.worm + lenth.of.worm + sex.of.host + 
                                         condylobasal.length.of.host.species + which.side)^2+ (1|id),
                    data = df, na.action = na.fail, REML = FALSE)

summary(mod_fullNUMHS)
summary(mod_fullNUMHS)$AIC[1]

### Model selection using AIC (fixed effects)
### Use ML not REML before running model averaging (Ben Bolker on Cross-Validated)
### Using AICc is even more controversial than AIC for GLMM -> use AIC to rank
### abund_T_randland_dred <- MuMIn::dredge(abund_T_randland, rank = "AIC")

### Let's susbet to models with an AIC within 2 of the lowest.
dredNUM <- dredge(mod_fullNUMHS, rank = "AIC")
dreddelta2NUM <- subset(dredNUM, delta<2)

### Most parsimonious model
mod_parsiNUM <- get.models(dreddelta2NUM, 10)[[1]]
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
                        legend.title = "",
                        ci.lvl = 0.95, se=TRUE)
best_plotNUM

################################################################################
######################## WITH HOST SPECIES #####################################
################################################################################

### Let's put id as random effects and use REML (Zuur et al 2009)

mod_norandNUMHS <- glmmTMB( number.of.worms ~ (sex.of.worm + lenth.of.worm + sex.of.host + 
                                               host.species + which.side)^2,
                          data = df, na.action = na.fail, REML = TRUE)

mod_randNUMHS <- glmmTMB( number.of.worms ~ (sex.of.worm + lenth.of.worm + sex.of.host + 
                                             host.species + which.side)^2 + (1|id),
                        data = df, na.action = na.fail, REML = TRUE)

summary(mod_norandNUMHS)$AIC[1]
summary(mod_randNUMHS)$AIC[1]
#NOTE: apparently we need to include id

### Full model summary with first level interactions

mod_fullNUMHS <- glmmTMB(number.of.worms ~ (sex.of.worm + lenth.of.worm + sex.of.host + 
                                           host.species + which.side)^2+ (1|id),
                       data = df, na.action = na.fail, REML = FALSE)

summary(mod_fullNUMHS)
summary(mod_fullNUMHS)$AIC[1]

### Model selection using AIC (fixed effects)
### Use ML not REML before running model averaging (Ben Bolker on Cross-Validated)
### Using AICc is even more controversial than AIC for GLMM -> use AIC to rank
### abund_T_randland_dred <- MuMIn::dredge(abund_T_randland, rank = "AIC")

### Let's susbet to models with an AIC within 2 of the lowest.
dredNUMHS <- dredge(mod_fullNUMHS, rank = "AIC")
dreddelta2NUMHS <- subset(dredNUMHS, delta<2)

### Most parsimonious model
mod_parsiNUMHS <- get.models(dreddelta2NUMHS, 3)[[1]]
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
best_plotNUMHS <- plot_model(mod_parsiNUMHS, type = "pred", terms = c("host.species","sex.of.host"),
                           title = '',
                           axis.title = c("Host species", expression("Predicted number of worms")),
                           legend.title = "",
                           ci.lvl = 0.95, se=TRUE)
best_plotNUMHS
