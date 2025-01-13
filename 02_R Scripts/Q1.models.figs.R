library(msm)
library(base)
library(modeest)
library(reshape)
library(dplyr)
library(tidyverse)
library(data.table)

setwd("/Users/lanie/Documents/UM/AR_Prod-BioDiv")
Q1.data<-read.csv("/Users/lanie/Documents/UM/AR_Prod-BioDiv/biodiversity_metadata3.csv") %>%
  select(-c(X))%>%
  rename(site = block)

###Q1 LMEs##################
library(lme4)
library(MuMIn)
library(partR2)
library(performance)
library(nortest)

stand <- function(X) { (X-mean(X,na.rm=T))/(2*sd(X,na.rm=T)) }#function
stand.data<-Q1.data %>% mutate(stand.rich = stand(sp_richness)) %>% 
  mutate(stand.div = stand(sp_evenness^2)) %>%
  mutate(stand.meanLmax = stand(log(Lmax + 1))) %>%
  mutate(stand.skew = stand(log(skewness + 10))) %>%
  mutate(stand.tp = stand(log(TP + 1))) %>%
  mutate(stand.FFG = stand(log(ffg_evenness +1))) %>%
  mutate(stand.biomass = stand(log(biomass +10)))

model1<-lmer(log(SP)~ stand.rich + stand.div + stand.meanLmax + stand.tp  + stand.skew + stand.FFG + (1|site:date), data = stand.data, REML = FALSE)
summary(model1)
rm=resid(model1)
print(shapiro.test(rm))
print(lillie.test(rm))
check_heteroscedasticity(model1)
options(na.action=na.fail)
d.model1 <- as.data.frame(subset(dredge(model1,extra = c("AIC","r.squaredGLMM")),delta < 2))
part.fit<-partR2(model1, partvars = c("stand.rich", "stand.div", "stand.meanLmax", "stand.tp", "stand.skew", 
                                      "stand.FFG"), nboot = 100, R2_type = "conditional", max_level = 1)
summary(part.fit)

model2<-lmer(log(n_supply) ~ stand.rich + stand.div + stand.meanLmax + stand.tp + stand.FFG  + stand.skew + (1|site:date), data = stand.data, REML=FALSE)
summary(model2)
rm2=resid(model2)
print(shapiro.test(rm2))
print(lillie.test(rm2))
check_heteroscedasticity(model2)
d.model2 <- as.data.frame(subset(dredge(model2,extra = c("AIC","r.squaredGLMM")),delta < 2))
part.fit2<-partR2(model2, partvars = c("stand.rich", "stand.div", "stand.meanLmax", "stand.tp", "stand.skew", 
                                       "stand.FFG"), nboot = 100, R2_type = "conditional", max_level = 1)
summary(part.fit2)

model3<-lmer(log(p_supply)~ stand.rich + stand.div + stand.meanLmax + stand.tp + stand.skew + stand.FFG  + (1|site:date), data = stand.data, REML = FALSE)
summary(model3)
rm3=resid(model3)
print(shapiro.test(rm3))
print(lillie.test(rm3))
check_heteroscedasticity(model3)
d.model3 <- as.data.frame(subset(dredge(model3,extra = c("AIC","r.squaredGLMM")),delta < 2))
part.fit3<-partR2(model3, partvars = c("stand.rich", "stand.div", "stand.meanLmax", "stand.tp", "stand.skew", 
                                       "stand.FFG"), nboot = 100, R2_type = "conditional", max_level = 1)
summary(part.fit3)

model4<-lmer(log(n.turnover)~ stand.rich + stand.div  + stand.meanLmax + stand.tp + stand.skew + stand.FFG + (1|site:date), data = stand.data, REML = FALSE)
summary(model4)
rm4=resid(model4)
print(shapiro.test(rm4))
print(lillie.test(rm4))
check_heteroscedasticity(model4)
d.model4 <- as.data.frame(subset(dredge(model4,extra = c("AIC","r.squaredGLMM")),delta < 2))
part.fit4<-partR2(model4, partvars = c("stand.rich", "stand.div", "stand.meanLmax", "stand.tp", "stand.skew", 
                                       "stand.FFG"), nboot = 100, R2_type = "conditional", max_level = 1)
summary(part.fit4)

model5<-lmer(log(p.turnover)~ stand.rich + stand.div + stand.meanLmax + stand.tp + stand.skew + stand.FFG + (1|site:date), data = stand.data, REML = FALSE)
summary(model5)
rm5=resid(model5)
print(shapiro.test(rm5))
print(lillie.test(rm5))
check_heteroscedasticity(model5)
d.model5 <- as.data.frame(subset(dredge(model5,extra = c("AIC","r.squaredGLMM")),delta < 2))
part.fit5<-partR2(model5, partvars = c("stand.rich", "stand.div", "stand.meanLmax", "stand.tp", "stand.skew", 
                                       "stand.FFG"), nboot = 100, R2_type = "conditional", max_level = 1)
summary(part.fit5)

subset.data<-stand.data %>% filter(!row_number() %in% c(62,73, 146, 59, 233, 227))
model6<-lmer(biomass.turnover ~ stand.rich + stand.div + stand.FFG + stand.skew + stand.meanLmax + stand.tp + (1|site:date), data = subset.data, REML = FALSE)
summary(model6)
rm6=resid(model6)
print(shapiro.test(rm6))
print(lillie.test(rm6))
check_heteroscedasticity(model6)
d.model6 <- as.data.frame(subset(dredge(model6,extra = c("AIC","r.squaredGLMM")),delta < 2))
part.fit6<-partR2(model6, partvars = c("stand.rich", "stand.div", "stand.meanLmax", "stand.tp", "stand.skew", 
                                       "stand.FFG"), nboot = 100, R2_type = "conditional", max_level = 1)
summary(part.fit6)

####Model Loop#############
response1 = cbind(log(Q1.data$n_supply),log(Q1.data$p_supply),log(Q1.data$SP),
                  log(Q1.data$n.turnover),log(Q1.data$p.turnover),(Q1.data$biomass.turnover), log(Q1.data$biomass))

response.names = c("N supply","P supply","Secondary Production", "N supply turnover", "P supply turnover", "Biomass turnover", "Biomass")
coef = matrix(ncol = 6,nrow = 7)
output = matrix(ncol = 6,nrow = 7)
output2 = matrix(ncol = 3,nrow = 7)
fredhead2 = list()
confints = list()

for (i in 1:7){
  
  model1=lmer(response1[,i] ~ stand.rich + stand.div + stand.meanLmax + stand.tp + stand.FFG + stand.skew + (1|site:date), data = stand.data,REML = FALSE)
  
  options(na.action=na.fail)
  dpp <- as.data.frame(subset(dredge(model1,extra = c("AIC","r.squaredGLMM")),delta < 2))
  fred=sapply(dpp,as.numeric)
  fred[is.na(fred)] = 0
  fred = data.frame(fred)
  
  fred$exp = exp(-0.5* fred$delta)
  fred$w = fred$exp/sum(fred$exp)
  fredhead2[[i]]= fred[1:6,]
  confints[[i]] = confint(model1,method = "boot")
  
  coef[i,]= colSums(data.frame(FFG = ifelse(abs(fred$stand.FFG)>0,fred$stand.FFG*fred$w,0),
                               Diversity = ifelse(abs(fred$stand.div)>0,fred$stand.div*fred$w,0),
                               Lmax= ifelse(abs(fred$stand.meanLmax)>0,fred$stand.meanLmax*fred$w,0),
                               Richness = ifelse(abs(fred$stand.rich)>0,fred$stand.rich*fred$w,0),
                               s.l = ifelse(abs(fred$stand.skew)>0,fred$stand.skew*fred$w,0),
                               TP = ifelse(abs(fred$stand.tp)>0,fred$stand.tp*fred$w,0)))
  
  
  output[i,]= colMeans(data.frame(FFG = ifelse(abs(fred$stand.FFG)>0,fred$w,0),
                                  Diversity = ifelse(abs(fred$stand.div)>0,fred$w,0),
                                  Lmax = ifelse(abs(fred$stand.meanLmax)>0,fred$w,0),
                                  Richness = ifelse(abs(fred$stand.rich)>0,fred$w,0),
                                  s.l = ifelse(abs(fred$stand.skew)>0,fred$w,0),
                                  TP = ifelse(abs(fred$stand.tp)>0,fred$w,0)))
  
  
  output2[i,] =  c(round(fred$r.squaredGLMM1[1],2),round(fred$r.squaredGLMM2[1],2),round(fred$AICc[1],2))
  
}

output3 = cbind(output,output2)
coef2 = cbind(coef,output2)
colnames(output3) = c('TG evenness', 'Sp evenness', expression("L[max]"), "Richness", "Skewness", 'Trophic level', "r2fixed","r2cond","AICc")
rownames(output3) = response.names
output3
colnames(coef2) = c('TG evenness', 'Sp evenness','Max length', "Richness", "Skewness", 'Trophic level', "r2fixed","r2cond","AICc")
rownames(coef2) = response.names
coef2
names(fredhead2) = response.names
names(confints) = response.names

coef3<-data.frame(coef2[,1:7])
graph.data<-data.frame(coef2) %>% rownames_to_column(var = "Response") %>% 
  pivot_longer(c(TG.evenness:Trophic.level), names_to = "Predictor", values_to = "Parameter_estimates") %>% 
  #pivot_longer(c(r2fixed:r2cond), names_to = "r2", values_to = "fit") %>%
  #dplyr::group_by(Response) %>%
  mutate(abs.est.par = abs(Parameter_estimates)) %>%
  mutate(Predictor = case_when(
    Predictor == "Richness" ~ "Sp richness",
    Predictor == "Max.length" ~ "Lmax",
    Predictor == "Sp.evenness" ~ "Sp evenness",
    Predictor == "TG.evenness" ~ "FG evenness",
    Predictor == "Trophic.level" ~ "Trophic level",
    T ~ Predictor
  )) %>%
  arrange(Response,desc(abs(Parameter_estimates))) #, .by_group = TRUE)

####Repeat for biomass turnover with subset.data#####
response.names4 = c("N supply","P supply","Secondary Production", "N supply turnover", "P supply turnover", "Biomass turnover", "Biomass")
coef4 = matrix(ncol = 6,nrow = 7)
output4 = matrix(ncol = 6,nrow = 7)
output5 = matrix(ncol = 3,nrow = 7)
fredhead5 = list()
confints4 = list()

model4=lmer(biomass.turnover~ stand.rich + stand.div + stand.meanLmax + stand.tp + stand.skew + stand.FFG + (1|site:date), data = subset.data, REML = FALSE)

options(na.action=na.fail)
dpp4 <- as.data.frame(subset(dredge(model4,extra = c("AIC","r.squaredGLMM")),delta < 2))
fred4=sapply(dpp4,as.numeric)
fred4[is.na(fred4)] = 0
fred4 = data.frame(fred4)

fred4$exp = exp(-0.5* fred4$delta)
fred4$w = fred4$exp/sum(fred4$exp)
fredhead5= fred4[1:6,]
confints4 = confint(model4,method = "boot")

coef4= colSums(data.frame(FFG = ifelse(abs(fred4$stand.FFG)>0,fred4$stand.FFG*fred4$w,0),
                          Diversity = ifelse(abs(fred4$stand.div)>0,fred4$stand.div*fred4$w,0),
                          Lmax= ifelse(abs(fred4$stand.meanLmax)>0,fred4$stand.meanLmax*fred4$w,0),
                          Richness = ifelse(abs(fred4$stand.rich)>0,fred4$stand.rich*fred4$w,0),
                          s.l = ifelse(abs(fred4$stand.skew)>0,fred4$stand.skew*fred4$w,0),
                          TP = ifelse(abs(fred4$stand.tp)>0,fred4$stand.tp*fred4$w,0)))


output4= colMeans(data.frame(FFG = ifelse(abs(fred4$stand.FFG)>0,fred4$w,0),
                             Diversity = ifelse(abs(fred4$stand.div)>0,fred4$w,0),
                             Lmax = ifelse(abs(fred4$stand.meanLmax)>0,fred4$w,0),
                             Richness = ifelse(abs(fred4$stand.rich)>0,fred4$w,0),
                             s.l = ifelse(abs(fred4$stand.skew)>0,fred4$w,0),
                             TP = ifelse(abs(fred4$stand.tp)>0,fred4$w,0)))


output5 =  c(round(fred4$r.squaredGLMM1[1],2),round(fred4$r.squaredGLMM2[1],2),round(fred4$AICc[1],2))

output6 = cbind(output4,output5)
coef6 = cbind(coef4,output5)
colnames(output6) = c('TG evenness', 'Sp evenness', expression("L[max]"), "Richness", "Skewness", 'Trophic level', "r2fixed","r2cond","AICc")
rownames(output6) = response.names4
output6
colnames(coef5) = c('TG evenness', 'Sp evenness','Max length', "Richness", "Skewness", 'Trophic level', "r2fixed","r2cond","AICc")
rownames(coef5) = response.names4
coef5
names(fredhead6) = response.names4
names(confints4) = response.names4

red.d.model6<-data.frame(cbind(coef4, output5)) %>%
  rownames_to_column(var = "Predictor") 

bio.turnover.graph.data<-red.d.model6 %>%
  mutate(abs.est.par = abs(coef4)) %>%
  mutate(Predictor = case_when(
    Predictor == "Richness" ~ "Sp richness",
    Predictor == "Diversity" ~ "Sp evenness",
    Predictor == "FFG" ~ "FG evenness",
    Predictor == "TP" ~ "Trophic level",
    Predictor == "s.l" ~ "Skewness",
    T ~ Predictor
  )) %>%
  arrange(Predictor,desc(abs(abs.est.par))) #, .by_group = TRUE)
"#4B0055" "#422C70" "#185086" "#007094" "#008E98" "#00A890" "#00BE7D" "#6CD05E" "#BBDD38" "#FDE333"
######Q1 figure#####################
graph.data$Predictor<-factor(graph.data$Predictor, levels = c("Sp evenness", "Sp richness", "FG evenness", "Skewness", "Lmax", "Trophic level"))

(N.supp.plot<-ggplot(data = graph.data %>% filter(Response == 'N supply'), aes(x=fct_rev(fct_reorder(Predictor, abs.est.par)),
                                                                               y = abs.est.par)) +
    geom_col(aes(fill = Predictor), color = "darkgray", width = 1, alpha = 0.7)+
    ylab("Weighted Parameter Estimates")+
    xlab("N supply")+
    scale_y_continuous(expand = c(0,0), limits = c(0, 1.1), breaks = c(0, 0.5, 1))+
    scale_x_discrete(expand = c(0,0))+
    scale_fill_manual(values = c("#4B0055", "#008E98", "#BBDD38", "#185086","#6CD05E","#FDE333"))+
    annotate("text", x = 1, y =0.5, label = "+", size = 7) +
    annotate("text", x = 2, y =0.5, label = "+", size = 7) +
    annotate("text", x = 3, y = 0.5, label = "-", size = 7) +
    annotate("text", x= 4.8, y=1, label = "paste(italic(R) ^ 2,\"m = 0.57\")", size =7, parse = TRUE) +
    annotate("text", x= 4.9, y=0.9, label = "paste(italic(R) ^ 2,\"c = 0.83\")", size =7, parse = TRUE) +
    theme_bw(base_size = 22)+
    theme(panel.border = element_blank(), axis.line.y = element_line(color = "darkgray", linewidth = 1), axis.line.x = element_line(color = "darkgray"), 
          panel.grid.major = element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(size =18),
          panel.grid.minor = element_blank(), legend.position = "none", legend.title = element_blank())
)

(P.supp.plot<-ggplot(data = graph.data %>% filter(Response == 'P supply'), aes(x=fct_rev(fct_reorder(Predictor, abs.est.par)),
                                                                               y = abs.est.par)) +
    geom_col(aes(fill = Predictor), color = "darkgray", width = 1, alpha = 0.7)+
    ylab("Weighted")+
    xlab("P supply")+
    scale_y_continuous(expand = c(0,0), limits = c(0, 1.1), breaks = c(0, 0.5, 1))+
    scale_x_discrete(expand = c(0,0))+
    scale_fill_manual(values = c("#4B0055", "#008E98", "#BBDD38", "#185086","#6CD05E","#FDE333"))+
    annotate("text", x = 1, y =0.5, label = "+", size = 7) +
    annotate("text", x = 2, y =0.5, label = "+", size = 7) +
    annotate("text", x = 3, y =0.5, label = "-", size = 7) +
    annotate("text", x = 4, y =0.5, label = "+", size = 7) +
    annotate("text", x= 4.8, y=1, label = "paste(italic(R) ^ 2,\"m = 0.49\")", size =7, parse = TRUE) +
    annotate("text", x= 4.9, y=0.9, label = "paste(italic(R) ^ 2,\"c = 0.76\")", size =7, parse = TRUE) +
    theme_bw(base_size = 22)+
    theme(panel.border = element_blank(), axis.line.y = element_blank(), axis.line.x = element_line(color = "darkgray"), panel.grid.major = element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.y = element_text(color = "white", margin = margin(0, 28, 0, 0)),
          panel.grid.minor = element_blank(), legend.position = "none")
)

(Prod.plot<-ggplot(data = graph.data %>% filter(Response == 'Secondary Production'), aes(x=fct_rev(fct_reorder(Predictor, abs.est.par)),
                                                                                         y = abs.est.par)) +
    geom_col(aes(fill = Predictor), color = "darkgray", width = 1, alpha = 0.7)+
    scale_y_continuous(expand = c(0,0), limits = c(0, 1.1), breaks = c(0, 0.5, 1))+
    ylab("Weighted")+
    xlab("Secondary production")+
    scale_x_discrete(expand = c(0,0))+
    scale_fill_manual(values = c("#4B0055", "#008E98", "#BBDD38", "#185086","#6CD05E","#FDE333"))+
    annotate("text", x = 1, y =0.5, label = "+", size = 7) +
    annotate("text", x = 2, y =0.5, label = "+", size = 7) +
    annotate("text", x = 3, y =0.5, label = "-", size = 7) +
    annotate("text", x = 4, y =0.5, label = "+", size = 7) +
    annotate("text", x= 4.8, y=1, label = "paste(italic(R) ^ 2,\"m = 0.65\")", size =7, parse = TRUE) +
    annotate("text", x= 4.9, y=0.9, label = "paste(italic(R) ^ 2,\"c = 0.80\")", size =7, parse = TRUE) +
    theme_bw(base_size = 22)+
    theme(panel.border = element_blank(), axis.line.y = element_blank(), axis.line.x = element_line(color = "darkgray"), panel.grid.major = element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.y = element_text(color = "white", margin = margin(0, 17, 0, 0)),
          panel.grid.minor = element_blank(), legend.position = "none")
)

(N.turnover.plot<-ggplot(data = graph.data%>% filter(Response == 'N supply turnover'), aes(x=fct_rev(fct_reorder(Predictor, abs.est.par)),
                                                                                         y = abs.est.par)) +
    geom_col(aes(fill = Predictor), color = "darkgray", width = 1, alpha = 0.7)+
    scale_y_continuous(expand = c(0,0), limits = c(0, 1.1), breaks = c(0, 0.5, 1))+
    ylab("Weighted")+
    xlab("N supply turnover")+
    scale_x_discrete(expand = c(0,0))+
    scale_fill_manual(values = c("#4B0055", "#008E98", "#BBDD38", "#185086","#6CD05E","#FDE333"))+
    annotate("text", x = 1, y =0.5, label = "-", size = 7) +
    annotate("text", x = 2, y =0.5, label = "-", size = 7) +
    annotate("text", x = 3, y =0.5, label = "-", size = 7) +
    annotate("text", x = 4, y =0.5, label = "-", size = 7) +
    annotate("text", x= 4.8, y=1, label = "paste(italic(R) ^ 2,\"m = 0.61\")", size =7, parse = TRUE) +
    annotate("text", x= 4.9, y=0.9, label = "paste(italic(R) ^ 2,\"c = 0.63\")", size =7, parse = TRUE) +
    theme_bw(base_size = 22)+
    theme(panel.border = element_blank(), axis.line.y = element_blank(), axis.line.x = element_line(color = "darkgray"), panel.grid.major = element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.y = element_text(color = "white", margin = margin(0, 32, 0, 0)),
          panel.grid.minor = element_blank(), legend.position = "none")
)

(P.turnover.plot<-ggplot(data = graph.data%>% filter(Response == 'P supply turnover'), aes(x=fct_rev(fct_reorder(Predictor, abs.est.par)),
                                                                                         y = abs.est.par)) +
    geom_col(aes(fill = Predictor), color = "darkgray", width = 1, alpha = 0.7)+
    scale_y_continuous(expand = c(0,0), limits = c(0, 1.1), breaks = c(0, 0.5, 1))+
    ylab("Weighted")+
    xlab("P supply turnover")+
    scale_x_discrete(expand = c(0,0))+
    scale_fill_manual(values = c("#4B0055", "#008E98", "#BBDD38", "#185086","#6CD05E","#FDE333"))+
    annotate("text", x = 1, y =0.5, label = "-", size = 7) +
    annotate("text", x = 2, y =0.5, label = "-", size = 7) +
    annotate("text", x= 4.8, y=1, label = "paste(italic(R) ^ 2,\"m = 0.48\")", size =7, parse = TRUE) +
    annotate("text", x= 4.9, y=0.9, label = "paste(italic(R) ^ 2,\"c = 0.70\")", size =7, parse = TRUE) +
    theme_bw(base_size = 22)+
    theme(panel.border = element_blank(), axis.line.y = element_blank(), axis.line.x = element_line(color = "darkgray"), panel.grid.major = element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.y = element_text(color = "white", margin = margin(0, 43, 0, 0)),
          panel.grid.minor = element_blank(), legend.position = "none")
)

bio.turnover.graph.data$Predictor<-factor(bio.turnover.graph.data$Predictor, levels = c("Sp evenness", "Sp richness", "FG evenness", "Skewness", "Lmax", "Trophic level"))

(bio.turnover.plot<-ggplot(data = bio.turnover.graph.data, aes(x=fct_rev(fct_reorder(Predictor, abs.est.par)),
                                                       y = abs.est.par)) +
    geom_col(aes(fill = Predictor), color = "darkgray", width = 1, alpha = 0.7)+
    scale_y_continuous(expand = c(0,0), limits = c(0, 1.1), breaks = c(0, 0.5, 1))+
    ylab("Weighted")+
    xlab("Biomass turnover")+
    scale_x_discrete(expand = c(0,0))+
    scale_fill_manual(values = c("#4B0055", "#008E98", "#BBDD38", "#185086","#6CD05E","#FDE333"))+
    annotate("text", x = 1, y =0.5, label = "-", size = 7) +
    annotate("text", x = 2, y =0.5, label = "-", size = 7) +
    annotate("text", x = 3, y =0.5, label = "-", size = 7) +
    annotate("text", x= 4.8, y=1, label = "paste(italic(R) ^ 2,\"m = 0.47\")", size =7, parse = TRUE) +
    annotate("text", x= 4.9, y=0.9, label = "paste(italic(R) ^ 2,\"c = 0.50\")", size =7, parse = TRUE) +
    theme_bw(base_size = 22)+
    theme(panel.border = element_blank(), axis.line.y = element_blank(), axis.line.x = element_line(color = "darkgray"), panel.grid.major = element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.y = element_text(color = "white", margin = margin(0, 18, 0, 0)),
          panel.grid.minor = element_blank(), legend.position = "top", legend.title = element_blank(), legend.key.size = unit(1, "cm"), legend.text = element_text(size = 20))+
    guides(fill = guide_legend(nrow = 1))
)

(n.rich.plot2<-ggplot(data = Q1.data, aes(x=sp_richness, y=log(n_supply)))+
    ylab("N supply")+
    xlab("Sp richness")+
    geom_point(color = "#008E98", size = 3, alpha = 0.7)+
    geom_smooth(method = "lm", color = "darkgray", fill = "lightgray", lwd = 1)+
    scale_x_continuous(expand = c(0.01,0.01))+
    scale_y_continuous(expand = c(0.01, 0.01))+
    theme_bw(base_size = 22)+
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text = element_text(size = 18), axis.title.y = element_text(margin = margin(0, 10, 0, 0))))

(p.rich.plot2<-ggplot(data = Q1.data, aes(x=sp_richness, y=log(p_supply)))+
    ylab("P supply")+
    xlab("Sp richness")+
    geom_point(color = "#008E98", size = 3, alpha = 0.7)+
    geom_smooth(method = "lm", color = "darkgray", fill = "lightgray", lwd = 1)+
    scale_x_continuous(expand = c(0.01,0.01))+
    scale_y_continuous(expand = c(0.01, 0.01))+
    theme_bw(base_size = 22)+
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text = element_text(size = 18)))

(prod.rich.plot2<-ggplot(data = Q1.data, aes(x=sp_richness, y=log(SP)))+
    ylab("Secondary production")+
    xlab("Sp richness")+
    geom_point(color = "#008E98", size = 3, alpha = 0.7)+
    geom_smooth(method = "lm", color = "darkgray", fill = "lightgray", lwd = 1)+
    scale_x_continuous(expand = c(0.01,0.01))+
    scale_y_continuous(expand = c(0.01, 0.01))+
    theme_bw(base_size = 22)+
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text = element_text(size = 18)))

(n.turnover.lmax.plot2<-ggplot(data = Q1.data, aes(x=log(Lmax +1), y=log(n.turnover)))+
    ylab("N supply turnover")+
    xlab("Lmax") +
    geom_point(color = "#6CD05E", size = 3, alpha = 0.7)+
    geom_smooth(method = "lm", color = "darkgray", fill = "lightgray", lwd = 1)+
    scale_x_continuous(expand = c(0.01,0.01))+
    scale_y_continuous(expand = c(0.01, 0.01), breaks = c(-8, -7, -6))+
    theme_bw(base_size = 22)+
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text = element_text(size = 18)))

(p.turnover.lmax.plot2<-ggplot(data = Q1.data, aes(x=log(Lmax +1), y=log(p.turnover)))+
    ylab("P supply turnover")+
    xlab("Lmax") +
    geom_point(color = "#6CD05E", size = 3, alpha = 0.7)+
    geom_smooth(method = "lm", color = "darkgray", fill = "lightgray", lwd = 1)+
    scale_x_continuous(expand = c(0.01,0.01))+
    scale_y_continuous(expand = c(0.01, 0.01))+
    theme_bw(base_size = 22)+
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text = element_text(size = 18)))  

(bio.turnover.lmax.plot2<-ggplot(data = subset.data, aes(x=log(Lmax+1), y=biomass.turnover))+
    ylab("Biomass turnover")+
    xlab("Lmax") +
    geom_point(color = "#6CD05E", size = 3, alpha = 0.7)+
    geom_smooth(method = "lm", color = "darkgray", fill = "lightgray", lwd = 1)+
    scale_x_continuous(expand = c(0.01,0.01))+
    scale_y_continuous(expand = c(0.01, 0.01))+
    theme_bw(base_size = 22)+
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text = element_text(size = 18), axis.title.y = element_text(margin = margin(0, 0, 0, 0))))

library(cowplot)
Q1.cowplot.grid<-plot_grid(N.supp.plot, P.supp.plot, Prod.plot, 
                           N.turnover.plot,  P.turnover.plot, bio.turnover.plot+ theme(legend.position = "none"), 
                           n.rich.plot2, p.rich.plot2, prod.rich.plot2, n.turnover.lmax.plot2, p.turnover.lmax.plot2, bio.turnover.lmax.plot2, 
                           nrow = 2, rel_heights = c(0.65, 0.5), rel_widths = c(1.01, 1.01, 1, 1.01, 1.05, 0.99))
legend<-get_legend(bio.turnover.plot + guides(color = guide_legend(nrow = 1)) + theme(legend.position = "right", legend.justification = "center"))
Q1.cowplot.grid2<-plot_grid(legend, Q1.cowplot.grid, nrow = 2, rel_heights = c(0.1, 0.75))

library(showtext)

showtext_auto()
showtext_opts(dpi = 300)

ggsave("Q1.plot.horz5.tiff", plot = Q1.cowplot.grid2, bg = "white", dpi = 300, width = 600, height = 250, units = "mm")

#####Conceptual figure######
subset.SP.norm<-stand.data %>% dplyr::filter(!row_number() %in% c(71, 41, 173, 200, 193, 198, 202, 205, 24, 264, 103))
site.mod<-lmer(log(n_supply)~ stand.biomass +(1|site:date), data = stand.data, REML = FALSE)
summary(site.mod)
ress=resid(site.mod)
print(shapiro.test(ress))
print(lillie.test(ress))
check_heteroscedasticity(site.mod)
options(na.action=na.fail)
d.site<- as.data.frame(subset(dredge(site.mod,extra = c("AIC","r.squaredGLMM"))))

site.mod.turn<-lmer(log(biomass.turnover)~stand.rich + (1|site:date), data = stand.data, REML = FALSE)
summary(site.mod.turn)
ress2=resid(site.mod.turn)
print(shapiro.test(ress2))
print(lillie.test(ress2))
check_heteroscedasticity(site.mod.turn)
d.site2<- as.data.frame(subset(dredge(site.mod.turn,extra = c("AIC","r.squaredGLMM"))))
ci = confint(site.mod.turn, method = "boot")

(biomass.concept.plot<-ggplot(data = subset.SP.norm) +
   geom_point(data = subset.SP.norm, aes(x= stand.rich, y = log(SP)), color = "#008E98",  alpha = 0.4, size =3) +
   geom_point(data = subset.SP.norm, aes(x= stand.biomass, y = log(SP)), color = "#A4A4A4",  alpha = 0.4, size =3) +
   geom_smooth(data = subset.SP.norm, aes(x= stand.biomass, y = log(SP)), color = "#A4A4A4", linetype = "solid", method = "lm", lwd = 2, se = FALSE) +
   geom_smooth(data = subset.SP.norm, aes(x= stand.rich, y = log(SP)), color = "#008E98", linetype = "dashed", method = "lm", lwd = 2, se = FALSE) +
   xlab("Attribute")+
   ylab("Secondary production")+
   theme_bw()+
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_text(size = 20, margin = margin(t = 0, r = 0, b = 0, l = 0)), axis.text.x = element_text(size = 18),
         panel.border = element_rect(linewidth = 1.2), axis.title.y = element_text(size = 20), axis.text.y = element_text(size = 18), margin = margin(0, 0, 0, 0)))

(SP.biomass.concept.plot<-ggplot(data = stand.data) +
    geom_point(data = stand.data, aes(x= stand.rich, y = log(biomass.turnover)), color = "#008E98", alpha = 0.4, size =3) +
    geom_smooth(data = stand.data, aes(x= stand.rich, y = log(biomass.turnover)), method = "lm", color = "#008E98", lwd = 2, se = FALSE) +
    xlab(NULL)+
    ylab(NULL)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 1.2),
          axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))
)
ggsave("biomass.concept.tiff", plot = biomass.concept.plot, bg = "white", dpi = 300, width = 250, height = 250, units = "mm")
ggsave("turnover.concept.tiff", plot = SP.biomass.concept.plot, bg = "white", dpi = 300, width = 250, height = 250, units = "mm")

concept.data<-read.csv("/Users/lanie/Documents/UM/AR_Prod-BioDiv/Book1.csv")
(concept.plot<-ggplot(concept.data, aes(x=Response, y = Effect, color = Attribute))+
    geom_point(size = 3, position = position_dodge(0.1))+
    geom_errorbar(aes(x = Response, ymin = Effect - SE, 
                     ymax = Effect + SE), 
                width = 0.3, linewidth = 0.5, inherit.aes = T, position=position_dodge(0.1))+
    scale_color_manual(name = "Attribute", values = c("#A4A4A4","#008E98")) +
    scale_y_continuous(limits = c(0,2))+
    ylab("Effect size")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(linewidth = 1.2),
          axis.text.y = element_text(size = 12))
)

#### Supplement figures####
(biomass.boxplot<-ggplot(data = stand.data, aes(x = site, y = stand.biomass, fill = site)) +
   geom_boxplot()+
   xlab("Site")+
   ylab("Biomass") +
   theme_bw(base_size = 20)+
   theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

(sp.rich.boxplot<-ggplot(data = stand.data, aes(x = site, y = stand.rich, fill = site)) +
    geom_boxplot()+
    xlab("Site")+
    ylab("Sp richness") +
    theme_bw(base_size = 20)+
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

(evenness.boxplot<-ggplot(data = stand.data, aes(x = site, y = stand.div, fill = site)) +
    geom_boxplot()+
    xlab("Site")+
    ylab("Sp evenness") +
    theme_bw(base_size = 20)+
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

(lmax.boxplot<-ggplot(data = stand.data, aes(x = site, y = stand.meanLmax, fill = site)) +
    geom_boxplot()+
    xlab("Site")+
    ylab("Lmax") +
    theme_bw(base_size = 20)+
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

(skewness.boxplot<-ggplot(data = stand.data, aes(x = site, y = stand.skew, fill = site)) +
    geom_boxplot()+
    xlab("Site")+
    ylab("Skewness") +
    theme_bw(base_size=20)+
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

(FGevenness.boxplot<-ggplot(data = stand.data, aes(x = site, y = stand.FFG, fill = site)) +
    geom_boxplot()+
    xlab("Site")+
    ylab("FG evenness") +
    theme_bw(base_size =20)+
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

(TL.boxplot<-ggplot(data = stand.data, aes(x = site, y = stand.tp, fill = site)) +
    geom_boxplot()+
    xlab("Site")+
    ylab("FG evenness") +
    theme_bw(base_size=20)+
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

supp.grid<-plot_grid(biomass.boxplot, sp.rich.boxplot, evenness.boxplot, lmax.boxplot, skewness.boxplot, 
                     FGevenness.boxplot, TL.boxplot)

ggsave("Supp.grid.tiff", plot = supp.grid, bg = "white", dpi = 300, width = 500, height = 400, units = "mm")

(biomass.rich.plot<- ggplot(Q1.data, (aes(x = sp_richness, y = log(biomass)))) +
    xlab("Sp richness")+
    ylab("Biomass")+
    geom_point(color = "#9C008B", size = 3, alpha = 0.7)+
    geom_smooth(method = "lm", color = "darkgray", fill = "lightgray", lwd = 1)+
    scale_x_continuous(expand = c(0.01,0.01))+
    scale_y_continuous(expand = c(0.01, 0.01))+
    theme_bw(base_size = 20)+
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text = element_text(size = 18), axis.title.y = element_text(margin = margin(0, 0, 0, 0)))) 

(n.supp.biomass.plot<- ggplot(stand.data, (aes(x = stand.biomass, y = log(SP)))) +
    xlab("Biomass")+
    ylab("SP")+
    geom_point(aes(color = SP.prod2.area), size = 3, alpha = 0.7)+
    geom_smooth(method = "lm", color = "darkgray", fill = "lightgray", lwd = 1)+
    #scale_x_continuous(expand = c(0.01,0.01))+
    #scale_y_continuous(expand = c(0.01, 0.01))+
    theme_bw(base_size = 20)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text = element_text(size = 18), axis.title.y = element_text(margin = margin(0, 0, 0, 0)))) 
