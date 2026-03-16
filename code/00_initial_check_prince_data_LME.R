

library(tidyverse)
library(here)
library(MuMIn)
library(effects)
library(ggeffects)
library(car)
library(lme4)
library(bbmle)
library(patchwork)




project_dir <- here::here()

raw_data <- file.path(paste(project_dir, "0_data", "Prince_2023_MK_precise.csv", sep = "/"))
mk <- read_csv(raw_data) #read in raw data
mk_filtered <- mk %>% filter(standardised.M.K<5 & standardised.M.K>0.3) %>% 
  filter(standardised.K<2) %>% 
  filter(Hoenig.M<2) %>% 
  drop_na(standardised.M.K, latitude, species) %>% 
  mutate(slat = as.vector(scale(latitude, scale = T)), sM.K = as.vector(scale(standardised.M.K, scale = T))) %>% 
  group_by(species) %>% 
  mutate(M_weight= 1/(sd(standardised.M.K)/sqrt(n()))) %>% 
  ungroup() %>% 
  drop_na(M_weight) %>% 
  filter(M_weight!="Inf")


#random slope and intercept
mod1 <- lmer(log(standardised.M.K) ~ slat*species+
                         (1+slat|species),
                       mk_filtered, REML = F, 
             #weights = M_weight, 
             control=lmerControl(optCtrl=list(maxfun=1e6, optimizer = "bobyqa")))

#random intercept only - use this model
mod0 <- lmer(log(standardised.M.K) ~ slat*species+
               (1|species),
             mk_filtered, REML = T, 
             #weights = M_weight, 
             control=lmerControl(optCtrl=list(maxfun=1e6, optimizer = "bobyqa")))

AICctab(mod1,mod0)
plot(mod0)

options(na.action=na.fail)
moddredge<-dredge(mod1,trace=2)

subset(moddredge, delta<3)

mod2 <- lmer(log(standardised.M.K) ~ slat+
               (1+slat|species),
             mk_filtered, REML = T, weights = M_weight, control=lmerControl(optCtrl=list(maxfun=1e6, optimizer = "bobyqa")))

summary(mod0)
Anova(mod0, type='III')
plot(mod2)
plot(allEffects(mod0))

moddf<- as.data.frame(Effect(c('slat'),mod2, xlevels=list(slat=c(seq(min(mk_filtered$slat), 
                                                                           max(mk_filtered$slat), by=0.1)))))

moddf <- moddf %>% mutate(transfit = exp(fit), 
                          transupper = exp(upper), 
                          translower = exp(lower))



ggplot(moddf, aes(slat, transfit))+
  geom_line()+
  geom_ribbon(aes(ymin = translower, ymax = transupper), alpha = .25)+
  theme_classic()


# Model to include in paper ####

project_dir <- here::here()

raw_data <- file.path(paste(project_dir, "0_data", "Prince_2023_MK_precise.csv", sep = "/"))
mk <- read_csv(raw_data) #read in raw data

mk_filtered <- mk %>% filter(standardised.M.K<=5 & standardised.M.K>=0.3) %>% 
  filter(standardised.K<3) %>% 
  filter(Hoenig.M<3) %>% 
  drop_na(standardised.M.K, latitude, species) %>% 
  mutate(slat = as.vector(scale(latitude, scale = T)), sM.K = as.vector(scale(standardised.M.K, scale = T))) %>% 
  group_by(species) %>% 
  mutate(n = n()) %>% 
  filter(n >1) %>% #only include species with > 1 M/K
  ungroup()

#random intercept only 
mod <- lmer(log(standardised.M.K) ~ slat+
               (1|species),
             mk_filtered, REML = T,  
             control=lmerControl(optCtrl=list(maxfun=1e6)))

plot(mod)
plot(allEffects(mod))
summary(mod)

#extract model data
moddf<- as.data.frame(Effect(c('slat'),mod, xlevels=list(slat=c(seq(min(mk_filtered$slat), 
                                                                     max(mk_filtered$slat), by=0.1)))))
#change back from log
moddf <- moddf %>% mutate(transfit = exp(fit), 
                          transupper = exp(upper), 
                          translower = exp(lower))



modplot <- ggplot(moddf, aes(slat, transfit))+
  geom_line(linewidth=1)+
  geom_ribbon(aes(ymin = translower, ymax = transupper), alpha = .25)+
  theme_classic()+
  theme(axis.title = element_text(size=13),
        axis.title.x = element_text(vjust = -0.2), 
        axis.title.y = element_text(vjust = 1.5),
        axis.text = element_text(size=10),
        plot.margin = margin(t=25, r=25, l=25, b=25))+
  labs(y="M/K ratio", x="Latitude (scaled)")+
  geom_rug(data = mk_filtered, aes(x = slat), sides = "b", inherit.aes = FALSE, linewidth=0.2)

#save plot 
plot_path <- file.path(paste(project_dir, "02_pipeline", 
                             "08_updated_results_steps_2024",
                             "output plots combined steps", sep = "/"))

ggsave("MK_Lat_model.png",
       plot = modplot, 
       path = plot_path, 
       width = 15, height = 15, units = "cm")


modplot2 <- ggplot(moddf, aes(slat, transfit))+
  geom_point(data=mk_filtered, aes(slat, standardised.M.K), shape=16, size= 1, colour="skyblue", alpha=0.35)+
  geom_line(linewidth=0.5)+
  geom_ribbon(aes(ymin = translower, ymax = transupper), alpha = .25)+
  theme_classic()+
  theme(axis.title = element_text(size=13),
        axis.title.x = element_text(vjust = -0.2), 
        axis.title.y = element_text(vjust = 1.5),
        axis.text = element_text(size=10),
        axis.line = element_line(linewidth = 0.3),
        plot.margin = margin(t=25, r=25, l=25, b=25))+
  labs(y="M/K ratio", x="Latitude (scaled)")


modplot2

#save plot 
plot_path <- file.path(paste(project_dir, "02_pipeline", 
                             "08_updated_results_steps_2024",
                             "output plots combined steps", sep = "/"))

ggsave("MK_Lat_model_dot.png",
       plot = modplot2, 
       path = plot_path, 
       width = 15, height = 15, units = "cm")

### model on K to include in paper #### 

project_dir <- here::here()

raw_data <- file.path(paste(project_dir, "0_data", "Prince_2023_MK_precise.csv", sep = "/"))
mk <- read_csv(raw_data) #read in raw data

mk_filtered <- mk %>% filter(standardised.M.K<=5 & standardised.M.K>=0.3) %>% 
  filter(standardised.K<3) %>% 
  filter(Hoenig.M<3) %>% 
  drop_na(standardised.K, latitude, species) %>% 
  mutate(slat = as.vector(scale(latitude, scale = T)), sK = as.vector(scale(standardised.K, scale = T))) %>% 
  group_by(species) %>% 
  mutate(n = n()) %>% 
  filter(n >1) %>% #only include species with > 1 M/K
  ungroup()

#random intercept only 
mod <- lmer(log(standardised.K) ~ slat+
              (1|species),
            mk_filtered, REML = T,  
            control=lmerControl(optCtrl=list(maxfun=1e6)))

plot(mod)
plot(allEffects(mod))
summary(mod)

#extract model data
moddf<- as.data.frame(Effect(c('slat'),mod, xlevels=list(slat=c(seq(min(mk_filtered$slat), 
                                                                    max(mk_filtered$slat), by=0.1)))))
#change back from log
moddf <- moddf %>% mutate(transfit = exp(fit), 
                          transupper = exp(upper), 
                          translower = exp(lower))



modplot <- ggplot(moddf, aes(slat, transfit))+
  geom_line(linewidth=1)+
  geom_ribbon(aes(ymin = translower, ymax = transupper), alpha = .25)+
  theme_classic()+
  theme(axis.title = element_text(size=13),
        axis.title.x = element_text(vjust = -0.2), 
        axis.title.y = element_text(vjust = 1.5),
        axis.text = element_text(size=10),
        plot.margin = margin(t=25, r=25, l=25, b=25))+
  labs(y=expression("K (year"^-1*")"), x="Latitude (scaled)")+
  geom_rug(data = mk_filtered, aes(x = slat), sides = "b", inherit.aes = FALSE, linewidth=0.2)

#save plot 
plot_path <- file.path(paste(project_dir, "02_pipeline", 
                             "08_updated_results_steps_2024",
                             "output plots combined steps", sep = "/"))

ggsave("K_Lat_model.png",
       plot = modplot, 
       path = plot_path, 
       width = 15, height = 15, units = "cm")


modplot2 <- ggplot(moddf, aes(slat, transfit))+
  geom_point(data=mk_filtered, aes(slat, standardised.K), shape=16, size= 1, colour="skyblue", alpha=0.35)+
  geom_line(linewidth=0.5)+
  geom_ribbon(aes(ymin = translower, ymax = transupper), alpha = .25)+
  theme_classic()+
  theme(axis.title = element_text(size=13),
        axis.title.x = element_text(vjust = -0.2), 
        axis.title.y = element_text(vjust = 1.5),
        axis.text = element_text(size=10),
        axis.line = element_line(linewidth = 0.3),
        plot.margin = margin(t=25, r=25, l=25, b=25))+
  labs(y=expression("K (year"^-1*")"), x="Latitude (scaled)")


modplot2

#save plot 
plot_path <- file.path(paste(project_dir, "02_pipeline", 
                             "08_updated_results_steps_2024",
                             "output plots combined steps", sep = "/"))

ggsave("K_Lat_model_dot.png",
       plot = modplot2, 
       path = plot_path, 
       width = 15, height = 15, units = "cm")
                      
