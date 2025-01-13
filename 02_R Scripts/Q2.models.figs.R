library(tidyverse)
library(lme4)
library(MuMIn)
library(partR2)
library(performance)
library(nortest)
devtools::install_github("hohenstein/remef")
library(remef)


####Q2 LMEs######################
Q2.data<-read.csv("/Users/lanie/Documents/UM/Code for Lanie for Global and YoMama/all.data2.csv") %>%
  separate(col = reef.year, into = c("reef", "year"), sep = "/") %>%
  mutate(site = block)

stand <- function(X) { (X-mean(X,na.rm=T))/(2*sd(X,na.rm=T)) }#function
stand.data2<-Q2.data %>% mutate(stand.rich = stand(log(rich_mean))) %>% 
  mutate(stand.div = stand(D_mean^2)) %>%
  mutate(stand.FFG = stand(log(D.ffg_mean + 1))) 

Q2.m1<-lmer(log(plant.production) ~ stand.rich + stand.div + stand.FFG + (1|site:year), data = stand.data2, REML = FALSE)
summary(Q2.m1)
Q2rm1=resid(Q2.m1)
print(shapiro.test(Q2rm1))
print(lillie.test(Q2rm1))
check_heteroscedasticity(Q2.m1)
options(na.action=na.fail)
Q2.model1 <- as.data.frame(subset(dredge(Q2.m1,extra = c("AIC","r.squaredGLMM")),delta < 2))

Q2.model1.1<-data.frame(Q2.model1[,1:8]) %>%
  summarise_at(vars(c(r.squaredGLMM1, r.squaredGLMM2)), list(mean = mean, sd = sd))

bar.order<-c("r.squaredGLMM1_mean", "r.squaredGLMM2_mean", "delta.r")

Q2.m1.graph.data<-Q2.model1.1 %>% 
  pivot_longer(c(r.squaredGLMM1_mean:r.squaredGLMM2_mean), names_to = "r2", values_to = "fit") %>%
  pivot_longer(c(r.squaredGLMM1_sd:r.squaredGLMM2_sd), names_to = "r2.sd", values_to = "fit.sd")
Q2.m1.graph.data2<-data.frame(Q2.m1.graph.data[-c(2:3),])

(Q2.m1.R2.plot2<-ggplot(Q2.m1.graph.data2) +
    geom_col(aes(x=r2, y = fit), fill = c("#D3D3D3", "#808080"), color = "black")+
    geom_errorbar(aes(x = r2, ymin = fit - fit.sd, ymax = fit + fit.sd), width = 0, linewidth = 0.5, inherit.aes = F) +
    scale_x_discrete(labels = c("r.squaredGLMM1_mean" = expression(paste("R"^"2", "m")), "r.squaredGLMM2_mean" = expression(paste("R"^"2", "c"))))+
    ylab("Variance explained")+
    xlab("Prediction")+
    scale_y_continuous(expand = c(0,0), limits = c(0, 0.75), breaks = c(0, 0.25, 0.5, 0.75), position = "right")+
    theme_bw(base_size = 22)+
    theme(axis.text.x=element_text(size = 17), axis.text.y.right =element_text(size = 17), panel.grid.major = element_blank(), plot.margin = unit(c(10,10,5,10), "pt"),
          panel.grid.minor = element_blank(), axis.title.x = element_text(color = "white"), axis.title.y.right = element_text(angle = 90, margin = margin(0, 0, 0, 6)))
)

Q2.m2<-lmer(plant.rich ~ stand.rich + stand.div + stand.FFG + (1|site:year), data = stand.data2, REML = FALSE)
summary(Q2.m2)
Q2rm2=resid(Q2.m2)
print(shapiro.test(Q2rm2))
print(lillie.test(Q2rm2))
check_heteroscedasticity(Q2.m2)
Q2.model2 <- as.data.frame(subset(dredge(Q2.m2,extra = c("AIC","r.squaredGLMM")),delta < 2))

Q2.model2.1<-data.frame(Q2.model2[,1:8]) %>%
  summarise_at(vars(c(r.squaredGLMM1, r.squaredGLMM2)), list(mean = mean, sd = sd))

Q2.m2.graph.data<-Q2.model2.1 %>% 
  pivot_longer(c(r.squaredGLMM1_mean:r.squaredGLMM2_mean), names_to = "r2", values_to = "fit") %>%
  pivot_longer(c(r.squaredGLMM1_sd:r.squaredGLMM2_sd), names_to = "r2.sd", values_to = "fit.sd")
Q2.m2.graph.data2<-data.frame(Q2.m2.graph.data[-c(2:3),])

(Q2.m2.R2.plot2<-ggplot(Q2.m2.graph.data2) +
    geom_col(aes(x=r2, y = fit), fill = c("#D3D3D3", "#808080"), color = "black")+
    geom_errorbar(aes(x = r2, ymin = fit - fit.sd, ymax = fit + fit.sd), width = 0, linewidth = 0.5, inherit.aes = F) +
    scale_x_discrete(labels = c("r.squaredGLMM1_mean" = expression(paste("R"^"2", "m")), "r.squaredGLMM2_mean" = expression(paste("R"^"2", "c"))))+
    ylab("Variance explained")+
    xlab(NULL)+
    scale_y_continuous(expand = c(0,0), limits = c(0, 0.75), breaks = c(0, 0.25, 0.5, 0.75), position = "right")+
    theme_bw(base_size = 22)+
    theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), axis.text.y.right =element_text(size = 17), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.title.y.right = element_text(angle = 90, margin = margin(0, 0, 0, 6)))
)

subset.data2<-stand.data2 %>% filter(!row_number() %in% c(36, 42)) 
Q2.m3<-lmer((plant.div +1)^2 ~ stand.rich + stand.div + stand.FFG + (1|site:year), data = subset.data2, REML = FALSE)
summary(Q2.m3)
Q2rm3=resid(Q2.m3)
print(shapiro.test(Q2rm3))
print(lillie.test(Q2rm3))
check_heteroscedasticity(Q2.m3)
Q2.model3 <- as.data.frame(subset(dredge(Q2.m3,extra = c("AIC","r.squaredGLMM")),delta < 2))

Q2.model3.1<-data.frame(Q2.model3[,1:8]) %>% 
  summarise_at(vars(c(r.squaredGLMM1, r.squaredGLMM2)), list(mean = mean, sd = sd))

Q2.m3.graph.data<-Q2.model3.1 %>% 
  pivot_longer(c(r.squaredGLMM1_mean:r.squaredGLMM2_mean), names_to = "r2", values_to = "fit") %>%
  pivot_longer(c(r.squaredGLMM1_sd:r.squaredGLMM2_sd), names_to = "r2.sd", values_to = "fit.sd")

Q2.m3.graph.data2<-data.frame(Q2.m3.graph.data[-c(2:3),])

(Q2.m3.R2.plot2<-ggplot(Q2.m3.graph.data2) +
    geom_col(aes(x=r2, y = fit), fill = c("#D3D3D3", "#808080"), color = "black")+
    scale_x_discrete(labels = c("r.squaredGLMM1_mean" = expression(paste("R"^"2", "m")), "r.squaredGLMM2_mean" = expression(paste("R"^"2", "c"))))+
    ylab("Variance explained")+
    xlab(NULL)+
    scale_y_continuous(expand = c(0,0), limits = c(0, 0.75), breaks = c(0, 0.25, 0.5, 0.75), position = "right")+
    theme_bw(base_size = 22)+
    theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), axis.text.y.right =element_text(size = 17), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.title.y.right = element_text(angle = 90, margin = margin(0, 0, 0, 6)))
)

#########Partial effects plots##############
Q2.m1.rich<- lmer(log(plant.production) ~ stand.div + stand.FFG + (1|site:year), data = stand.data2, REML = FALSE)
Q2.m1.rich2<- lmer(stand.rich ~ stand.div + stand.FFG + (1|site:year), data = stand.data2, REML = FALSE)
x<-resid(Q2.m1.rich2)
y<-resid(Q2.m1.rich)
plot(x,y, xlab = "richness", ylab = "log(Producer productivity)", type = "n")
part <- lm(y~x)
wx = par("usr")[1:2]
new.x = seq(wx[1],wx[2],len=45)
pred = predict(part, new=data.frame(x=new.x), interval="conf")
lines(new.x,pred[,"fit"],lwd=2)
lines(new.x,pred[,"lwr"],lty=3)
lines(new.x,pred[,"upr"],lty=3)
points(x,y,pch=16,col="steelblue")

(Q2.m1.rich.plot<-ggplot(data = NULL, aes(x=x, y=y))+
    geom_point(size = 4, pch =21, fill = "#008E98", color = "darkgrey", alpha = 0.7)+
    geom_line(aes(x=new.x, y= pred[,"fit"]), lwd = 1, color = "darkgrey", linetype = "dashed")+
    geom_ribbon(aes(x = new.x, ymin = pred[,"lwr"], ymax = pred[,"upr"]), fill = "lightgray", alpha = 0.4) + 
    xlab("Fish richness")+
    ylab("Seagrass production")+
    scale_x_continuous(expand = c(0.0,0.0)) +
    scale_y_continuous(expand = c(0.025,0.025)) +
    theme_bw(base_size = 22)+
    theme(panel.grid.major = element_blank(), axis.text.x = element_text(size = 17), axis.text.y = element_text(size = 17), axis.title.x = element_text(margin = margin(t = 6, r = 0, b = 0, l = 0)),
          panel.grid.minor = element_blank()))

Q2.m1.div<- lmer(log(plant.production) ~ stand.rich + stand.FFG + (1|site:year), data = stand.data2, REML = FALSE)
Q2.m1.div2<- lmer(stand.div ~ stand.rich + stand.FFG + (1|site:year), data = stand.data2, REML = FALSE)
x2<-resid(Q2.m1.div2)
y2<-resid(Q2.m1.div)
plot(x2,y2, xlab = "diversity", ylab = "log(Plant production)", type = "n")
part2 <- lm(y2~x2)
wx2 = par("usr")[1:2]
new.x2 = seq(wx2[1],wx2[2],len=45)
pred2 = predict(part2, new=data.frame(x2=new.x2), interval="conf")
lines(new.x2,pred2[,"fit"],lwd=2)
lines(new.x2,pred2[,"lwr"],lty=3)
lines(new.x2,pred2[,"upr"],lty=3)
points(x2,y2,pch=16,col="steelblue")

(Q2.m1.div.plot<-ggplot(data = NULL, aes(x=x2, y=y2))+
    geom_point(size = 4, pch= 21, fill = "#4B0055", color = "darkgrey", alpha = 0.7)+
    geom_line(aes(x=new.x2, y= pred2[,"fit"]), lwd = 1, color = "darkgrey", linetype = "dashed")+
    geom_ribbon(aes(x = new.x2, ymin = pred2[,"lwr"], ymax = pred2[,"upr"]), fill = "lightgray", alpha = 0.4) + 
    xlab("Fish evenness")+
    ylab(NULL)+
    scale_x_continuous(expand = c(0.0,0.0), breaks = c(-0.6, -0.3, 0.00, 0.3)) +
    scale_y_continuous(expand = c(0.025,0.025)) +
    theme_bw(base_size = 22)+
    theme(panel.grid.major = element_blank(), 
          axis.ticks.y=element_blank(), axis.text.x=element_text(size = 17), axis.text.y = element_blank(),  
          panel.grid.minor = element_blank(),axis.title.x = element_text(margin = margin(6, 0, 0, 0))))

Q2.m1.ffg<- lmer(log(plant.production) ~ stand.rich + stand.div + (1|site:year), data = stand.data2, REML = FALSE)
Q2.m1.ffg2<- lmer(stand.FFG ~ stand.rich + stand.div + (1|site:year), data = stand.data2, REML = FALSE)
x4<-resid(Q2.m1.ffg2)
y4<-resid(Q2.m1.ffg)
plot(x4, y4, xlab = "Max length", ylab = "log(Plant production)", type = "n")
part4<- lm(y4~x4)
wx4 = par("usr")[1:2]
new.x4 = seq(wx4[1],wx4[2],len=45)
pred4 = predict(part4, new=data.frame(x4=new.x4), interval="conf")
lines(new.x4,pred4[,"fit"],lwd=2)
lines(new.x4,pred4[,"lwr"],lty=3)
lines(new.x4,pred4[,"upr"],lty=3)
points(x4,y4,pch=16,col="steelblue")

(Q2.m1.ffg.plot<-ggplot(data = NULL, aes(x=x4, y=y4))+
    geom_point(size = 4, pch =21, fill = "#BBDD38", color = "darkgrey", alpha = 0.7)+
    geom_ribbon(aes(x = new.x4, ymin = pred4[,"lwr"], ymax = pred4[,"upr"]), fill = "lightgray", alpha = 0.4) + 
    geom_line(aes(x=new.x4, y= pred4[,"fit"]), lwd = 1, color = "darkgrey", linetype = "dashed")+
    xlab("Fish FG evenness")+
    ylab(NULL)+
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0.025,0.025)) +
    theme_bw(base_size = 22)+
    theme(panel.grid.major = element_blank(), axis.title.x = element_text(margin = margin(t = 6, r = 0, b = 0, l = 0)),
          axis.ticks.y=element_blank(), axis.text.y=element_blank(),  axis.text.x = element_text(size = 17),
          panel.grid.minor = element_blank()))

Q2.m2.rich<- lmer(plant.rich ~ stand.div + stand.FFG + (1|site:year), data = stand.data2, REML = FALSE)
Q2.m2.rich2<- lmer(stand.rich ~ stand.div + stand.FFG + (1|site:year), data = stand.data2, REML = FALSE)
x5<-resid(Q2.m2.rich2)
y5<-resid(Q2.m2.rich)
plot(x5,y5, xlab = "richness", ylab = "Plant richness", type = "n")
part5 <- lm(y5~x5)
wx5 = par("usr")[1:2]
new.x5 = seq(wx5[1],wx5[2],len=45)
pred5 = predict(part5, new=data.frame(x5=new.x5), interval="conf")
lines(new.x5,pred5[,"fit"],lwd=2)
lines(new.x5,pred5[,"lwr"],lty=3)
lines(new.x5,pred5[,"upr"],lty=3)
points(x5,y5,pch=16,col="steelblue")

(Q2.m2.rich.plot<-ggplot(data = NULL, aes(x=x5, y=y5))+
    geom_point(size = 4, pch = 21, fill = "#008E98", color = "darkgrey", alpha = 0.7)+
    geom_ribbon(aes(x = new.x5, ymin = pred5[,"lwr"], ymax = pred5[,"upr"]), fill = "lightgray", alpha = 0.4) + 
    geom_line(aes(x=new.x5, y= pred5[,"fit"]), lwd = 1, color = "darkgrey", linetype = "dashed")+
    xlab(NULL)+
    ylab("Producer richness")+
    scale_x_continuous(expand = c(0.0,0.0)) +
    scale_y_continuous(expand = c(0.025,0.025)) +
    theme_bw(base_size = 22)+
    theme(panel.grid.major = element_blank(), axis.text.y = element_text(size = 17), axis.ticks.x=element_blank(), axis.text.x=element_blank(), 
          panel.grid.minor = element_blank(), axis.title.y = element_text(margin = margin(0, 17, 0, 0)),))

Q2.m2.div<- lmer(plant.rich ~ stand.rich + stand.FFG + (1|site:year), data = stand.data2, REML = FALSE)
Q2.m2.div2<- lmer(stand.div ~ stand.rich + stand.FFG + (1|site:year), data = stand.data2, REML = FALSE)
x6<-resid(Q2.m2.div2)
y6<-resid(Q2.m2.div)
plot(x6,y6, xlab = "diversity", ylab = "Plant richness", type = "n")
part6 <- lm(y6~x6)
wx6 = par("usr")[1:2]
new.x6 = seq(wx6[1],wx6[2],len=45)
pred6 = predict(part6, new=data.frame(x6=new.x6), interval="conf")
lines(new.x6,pred6[,"fit"],lwd=2)
lines(new.x6,pred6[,"lwr"],lty=3)
lines(new.x6,pred6[,"upr"],lty=3)
points(x6,y6,pch=16,col="steelblue")

(Q2.m2.div.plot<-ggplot(data = NULL, aes(x=x6, y=y6))+
    geom_point(size = 4, pch = 21, fill = "#4B0055", color = "darkgrey", alpha = 0.7)+
    xlab(NULL)+
    ylab(NULL)+
    scale_x_continuous(expand = c(0.025,0.025)) +
    scale_y_continuous(expand = c(0.025,0.025)) +
    theme_bw(base_size = 22)+
    theme(panel.grid.major = element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.y=element_blank(), axis.text.y=element_blank(), panel.grid.minor = element_blank()))

Q2.m2.ffg<- lmer(plant.rich ~ stand.rich + stand.div + (1|site:year), data = stand.data2, REML = FALSE)
Q2.m2.ffg2<- lmer(stand.FFG ~ stand.rich + stand.div + (1|site:year), data = stand.data2, REML = FALSE)
x8<-resid(Q2.m2.ffg2)
y8<-resid(Q2.m2.ffg)
plot(x8, y8, xlab = "Max length", ylab = "Plant richness", type = "n")
part8<- lm(y8~x8)
wx8 = par("usr")[1:2]
new.x8 = seq(wx8[1],wx8[2],len=45)
pred8 = predict(part8, new=data.frame(x8=new.x8), interval="conf")
lines(new.x8,pred8[,"fit"],lwd=2)
lines(new.x8,pred8[,"lwr"],lty=3)
lines(new.x8,pred8[,"upr"],lty=3)
points(x8,y8,pch=16,col="steelblue")

(Q2.m2.ffg.plot<-ggplot(data = NULL, aes(x=x8, y=y8))+
    geom_point(size = 4, pch = 21, fill = "#BBDD38", color = "darkgrey", alpha = 0.7)+
    geom_ribbon(aes(x = new.x8, ymin = pred8[,"lwr"], ymax = pred8[,"upr"]), fill = "lightgray", alpha = 0.4) +
    geom_line(aes(x=new.x8, y= pred8[,"fit"]), lwd = 1, color = "darkgrey", linetype = "dashed")+
    xlab(NULL) +
    ylab(NULL) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0.025,0.025)) +
    theme_bw(base_size = 22)+
    theme(panel.grid.major = element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.y=element_blank(), axis.text.y=element_blank(), panel.grid.minor = element_blank()))

Q2.m3.rich<- lmer((plant.div +1)^2 ~ stand.div  + stand.FFG + (1|site:year), data = subset.data2, REML = FALSE)
Q2.m3.rich2<- lmer(stand.rich ~ stand.div +  stand.FFG + (1|site:year), data = subset.data2, REML = FALSE)
x9<-resid(Q2.m3.rich2)
y9<-resid(Q2.m3.rich)
plot(x9,y9, xlab = "richness", ylab = "Plant diversity", type = "n")
part9 <- lm(y9~x9)
wx9 = par("usr")[1:2]
new.x9 = seq(wx9[1],wx9[2],len=43)
pred9 = predict(part9, new=data.frame(x9=new.x9), interval="conf")
lines(new.x9,pred9[,"fit"],lwd=2)
lines(new.x9,pred9[,"lwr"],lty=3)
lines(new.x9,pred9[,"upr"],lty=3)
points(x9,y9,pch=16,col="steelblue")

(Q2.m3.rich.plot<-ggplot(data = NULL, aes(x=x9, y=y9))+
    geom_point(size = 4, pch = 21, fill = "#008E98", color = "darkgrey", alpha = 0.7)+
    geom_ribbon(aes(x = new.x9, ymin = pred9[,"lwr"], ymax = pred9[,"upr"]), fill = "lightgray", alpha = 0.4) + 
    geom_line(aes(x=new.x9, y= pred9[,"fit"]), lwd = 1, color = "darkgrey")+
    xlab(NULL)+
    ylab("Producer evenness")+
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0.025,0.025)) +
    theme_bw(base_size = 22)+
    theme(panel.grid.major = element_blank(), axis.text.y = element_text(size = 17),
          axis.ticks.x=element_blank(), axis.text.x=element_blank(),
          panel.grid.minor = element_blank()))

Q2.m3.div<- lmer((plant.div+1)^2 ~ stand.rich + stand.FFG + (1|site:year), data = subset.data2, REML = FALSE)
Q2.m3.div2<- lmer(stand.div ~ stand.rich + stand.FFG + (1|site:year), data = subset.data2, REML = FALSE)
x10<-resid(Q2.m3.div2)
y10<-resid(Q2.m3.div)
plot(x10,y10, xlab = "diversity", ylab = "Plant diversity", type = "n")
part10 <- lm(y10~x10)
wx10 = par("usr")[1:2]
new.x10 = seq(wx10[1],wx10[2],len=43)
pred10 = predict(part10, new=data.frame(x10=new.x10), interval="conf")
lines(new.x10,pred10[,"fit"],lwd=2)
lines(new.x10,pred10[,"lwr"],lty=3)
lines(new.x10,pred10[,"upr"],lty=3)
points(x10,y10,pch=16,col="steelblue")

(Q2.m3.div.plot<-ggplot(data = NULL, aes(x=x10, y=y10))+
    geom_point(size = 4, pch = 21, fill = "#4B0055", color = "darkgrey", alpha = 0.7)+
    xlab(NULL)+
    ylab(NULL)+
    scale_x_continuous(expand = c(0.025,0.025)) +
    scale_y_continuous(expand = c(0.025,0.025)) +
    theme_bw(base_size = 22)+
    theme(panel.grid.major = element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.y=element_blank(), axis.text.y=element_blank(), panel.grid.minor = element_blank()))

Q2.m3.ffg<- lmer((plant.div+1)^2 ~ stand.rich + stand.div + (1|site:year), data = subset.data2, REML = FALSE)
Q2.m3.ffg2<- lmer(stand.FFG ~ stand.rich + stand.div + (1|site:year), data = subset.data2, REML = FALSE)
x12<-resid(Q2.m3.ffg2)
y12<-resid(Q2.m3.ffg)
plot(x12, y12, xlab = "Fish Max length", ylab = "Plant diversity", type = "n")
part12<- lm(y12~x12)
wx12 = par("usr")[1:2]
new.x12 = seq(wx12[1],wx12[2],len=43)
pred12 = predict(part12, new=data.frame(x12=new.x12), interval="conf")
lines(new.x12,pred12[,"fit"],lwd=2)
lines(new.x12,pred12[,"lwr"],lty=3)
lines(new.x12,pred12[,"upr"],lty=3)
points(x12,y12,pch=16,col="steelblue")

(Q2.m3.ffg.plot<-ggplot(data = NULL, aes(x=x12, y=y12))+
    geom_point(size = 4, pch = 21, fill = "#BBDD38", color = "darkgrey", alpha = 0.7)+
    xlab(NULL)+
    ylab(NULL)+
    scale_x_continuous(expand = c(0.025,0.025)) +
    scale_y_continuous(expand = c(0.025,0.025)) +
    theme_bw(base_size = 22)+
    theme(panel.grid.major = element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.y=element_blank(), axis.text.y=element_blank(), panel.grid.minor = element_blank()))

library(showtext)

showtext_auto()
showtext_opts(dpi = 300)

Q2cowplot.grid<-plot_grid(Q2.m2.rich.plot, Q2.m2.div.plot, Q2.m2.ffg.plot, Q2.m2.R2.plot2, Q2.m3.rich.plot, Q2.m3.div.plot, Q2.m3.ffg.plot, Q2.m3.R2.plot2,
                        Q2.m1.rich.plot, Q2.m1.div.plot, Q2.m1.ffg.plot, Q2.m1.R2.plot2, ncol = 4,
                        rel_widths = c(0.79, 0.65, 0.65, 0.45), rel_heights = c(0.85, 0.85, 1), scale = 0.99)

ggsave("Q2.plot.grid11.tiff", plot = Q2cowplot.grid, bg = "white", dpi = 300, width = 400, height = 250, units = "mm")


#######Sediment nutrients analysis################
library(lmPerm)
library(rcompanion)

####Site YM2 = fertilized reefs###
sed.nut.data<-Q2.data %>% filter(sed.n > 0) %>%
  mutate(block = case_when(
    block == "YoMama" ~ "YM",
    block == "YoMama2" ~ "YM - fert",
    T ~ block
  ))

div.sub<-sed.nut.data %>% filter(!row_number() %in% c(12, 18))
sea.div<-lm(plant.div ~ sed.n + site, data = sed.nut.data)
summary(sea.div)
resss=resid(sea.div)
print(shapiro.test(resss))
check_heteroscedasticity(sea.div)

sea.prod<-lm(log(plant.production) ~ sed.n + site, data = sed.nut.data)
summary(sea.prod)
resss=resid(sea.prod)
print(shapiro.test(resss))
check_heteroscedasticity(sea.prod)
options(na.action=na.fail)
Q2.model1 <- as.data.frame(subset(dredge(sea.prod,extra = c("AIC","r.squaredGLMM")),delta < 2))

sea.rich<-lm(plant.rich ~ sed.p + site, data = sed.nut.data)
summary(sea.rich)
resss=resid(sea.rich)
print(shapiro.test(resss))
check_heteroscedasticity(sea.rich)
Q2.model2 <- as.data.frame(subset(dredge(sea.rich,extra = c("AIC","r.squaredGLMM")),delta < 2))

perm.anova1<-aovp(sed.n ~ site, sed.nut.data)
summary(perm.anova1)
pairwisePermutationTest(sed.n ~ site, sed.nut.data)

perm.anova2<-aovp(sed.p ~ site, sed.nut.data)
summary(perm.anova2)
pairwisePermutationTest(sed.p ~ site, sed.nut.data)

###########Sediment nutrient plots####################
(n.boxplot<-ggplot(data = sed.nut.data, aes(x = block, y = sed.n, fill = block)) +
    geom_boxplot()+
    scale_fill_manual(values = c("#81C07A", "#8175CB", "#A79FE1", "#CAC5F3", "#E6C6F3"))+
    xlab("Site")+
    ylab("Sediment %N") +
    annotate("text", x = 1, y =0.127, label = "a", size = 3, fontface = "bold") +
    annotate("text", x = 2, y =0.0335, label = "ab", size = 3, fontface = "bold") +
    annotate("text", x = 3, y =0.0663, label = "a", size = 3, fontface = "bold") +
    annotate("text", x = 4, y =0.0575, label = "ac", size = 3, fontface = "bold") +
    annotate("text", x = 5, y =0.0637, label = "ac", size = 3, fontface = "bold") +
    theme_bw()+
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

ggsave("sed.n.boxplot.tiff", plot = n.boxplot, dpi = 150, width = 75, height = 85, units = "mm")

(p.boxplot<-ggplot(data = sed.nut.data, aes(x = block, y = sed.p, fill = block)) +
    geom_boxplot()+
    scale_fill_manual(values = c("#81C07A", "#8175CB", "#A79FE1", "#CAC5F3", "#E6C6F3"))+
    xlab("Site")+
    ylab("Sediment %P") +
    annotate("text", x = 1, y =0.0218, label = "a", size = 3, fontface = "bold") +
    annotate("text", x = 2, y =0.0108, label = "b", size = 3, fontface = "bold") +
    annotate("text", x = 3, y =0.0117, label = "b", size = 3, fontface = "bold") +
    annotate("text", x = 4, y =0.0117, label = "b", size = 3, fontface = "bold") +
    annotate("text", x = 5, y =0.0138, label = "b", size = 3, fontface = "bold") +
    theme_bw()+
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

ggsave("sed.p.boxplot.tiff", plot = p.boxplot, dpi = 150, width = 75, height = 85, units = "mm")
