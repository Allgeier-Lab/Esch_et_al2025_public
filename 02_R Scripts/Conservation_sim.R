rm(list = ls())
Q1.sim.data<-read.csv("/Users/lanie/Documents/UM/AR_Prod-BioDiv/Q1.sim.data.csv")
Q2.sim.data<-read.csv("/Users/lanie/Documents/UM/AR_Prod-BioDiv/Q2.sim.data.csv")

f1 <- Q1.sim.data[,-1]
f2 <- Q2.sim.data[,-1]
names(f2)[c(3,6)] <- c("richness", "biomass")


f1 = data.frame(richness = tapply(f1$richness, f1$reef, mean),
           biomass = tapply(f1$biomass, f1$reef, mean),
           SP = tapply(f1$SP, f1$reef, mean),
           country = tapply(f1$country, f1$reef, unique),
           reef = tapply(f1$reef, f1$reef, unique))



f1.list <- list(f1,
                f1[which(f1$country == "Bahamas"),],
                f1[which(f1$country == "Haiti"),] )

f2.list <- list(f2,
                f2[which(f2$country == "Bahamas"),],
                f2[which(f2$country == "Haiti"),] )

#i=1
sims <- 100 #number of iterations
tophowmany <- 3 # this is the number of reefs being conserved
region.list.RR <- list() #makes list to compile regions
region.list.Perc <- list() #makes list to compile regions

for(r in 1:3){
  #r=2
  #isolate the region
  f <- f1.list[[r]]
  fpp <- f2.list[[r]]

  #isolate the strategy to bootstrap (richness, biomass, SP)
  strategies <- c("biomass", "richness", "SP", "PP")
  strategies.list.RR <- list()
  strategies.list.Perc <- list()
  
  for(strat in 1:3){
    #strat=1
    RR.mat <- matrix(nrow = length(strategies), ncol = sims)
    Perc.mat <- matrix(nrow = length(strategies), ncol = sims)
        
        for(its in 1:sims){
          #for the others
          s <- f[-sample(1:length(f[,1]),1),]
          t3 <- s[order(s[,strategies[strat]],decreasing = T),][1:tophowmany,]
          tlast <- s[order(s[,strategies[strat]],decreasing = T),][(tophowmany+1):length(s[,1]),]
          
          #for the PP
          spp <- fpp[-sample(1:length(fpp[,1]),1),]
          t3pp <- spp[order(spp[,strategies[strat]],decreasing = T),][1:tophowmany,]
          tlastpp <- spp[order(spp[,strategies[strat]],decreasing = T),][(tophowmany+1):length(spp[,1]),]
          
          #load the RRs into the RR matrix
          RR.mat[1,its] <- log(mean(t3$biomass)) - log(mean(tlast$biomass))
          RR.mat[2,its] <- log(mean(t3$richness)) - log(mean(tlast$richness))
          RR.mat[3,its] <- log(mean(t3$SP)) - log(mean(tlast$SP))
          RR.mat[4,its] <- log(mean(t3pp$PP)) - log(mean(tlastpp$PP))
          #load the RRs into the Perc matrix
          Perc.mat[1,its] <- (mean(t3$biomass)/mean(tlast$biomass))
          Perc.mat[2,its] <- (mean(t3$richness)/mean(tlast$richness))
          Perc.mat[3,its] <- (mean(t3$SP)/mean(tlast$SP))
          Perc.mat[4,its] <- (mean(t3pp$PP)/mean(tlastpp$PP))
        }
          #this compiles the matrix outputs for the three strategies
          strategies.list.RR[[strat]] <- RR.mat
          strategies.list.Perc[[strat]] <- Perc.mat
          
      }
  #if(strat == 4) PP
  RR.matpp <- matrix(nrow = length(strategies), ncol = sims)
  Perc.matpp <- matrix(nrow = length(strategies), ncol = sims)
    for(its in 1:sims){
      #its = 1
      
      #for the PP
      spp <- fpp[-sample(1:length(fpp[,1]),1),]
      t3pp <- spp[order(spp[,strategies[4]],decreasing = T),][1:tophowmany,]
      tlastpp <- spp[order(spp[,strategies[4]],decreasing = T),][(tophowmany+1):length(spp[,1]),]
      
      #load the RRs into the RR matrix
      RR.matpp[1,its] <- log(mean(t3pp$biomass)) - log(mean(tlastpp$biomass))
      RR.matpp[2,its] <- log(mean(t3pp$richness)) - log(mean(tlastpp$richness))
      RR.matpp[3,its] <- log(mean(t3pp$SP)) - log(mean(tlastpp$SP))
      RR.matpp[4,its] <- log(mean(t3pp$PP)) - log(mean(tlastpp$PP))
      #load the RRs into the Perc matrix
      Perc.matpp[1,its] <- (mean(t3pp$biomass)/mean(tlastpp$biomass))
      Perc.matpp[2,its] <- (mean(t3pp$richness)/mean(tlastpp$richness))
      Perc.matpp[3,its] <- (mean(t3pp$SP)/mean(tlastpp$SP))
      Perc.matpp[4,its] <- (mean(t3pp$PP)/mean(tlastpp$PP))
    }
    #this compiles the matrix outputs for the three strategies
    strategies.list.RR[[4]] <- RR.matpp
    strategies.list.Perc[[4]] <- Perc.matpp
  
  #this compiles the stragetics within the regions
  region.list.RR[[r]] <- strategies.list.RR
  region.list.Perc[[r]] <- strategies.list.Perc
    
}

names(region.list.RR) <- c("All", "Bahamas", "Haiti")
names(region.list.Perc) <- c("All", "Bahamas", "Haiti")



#plot these thangs up RR

#quick summary

#All: strategy = Biomass - responses c("biomass", "richness", "SP")
apply(region.list.RR[[3]][[1]],1,mean)
#All: strategy = Richenss - responses c("biomass", "richness", "SP")
apply(region.list.RR[[3]][[2]],1,mean)
#All: strategy = SP - responses c("biomass", "richness", "SP")
apply(region.list.RR[[3]][[3]],1,mean)

#Bahamas: strategy = Biomass - responses c("biomass", "richness", "SP")
apply(region.list.RR[[1]][[1]],1,mean)
#Bahamas: strategy = Richenss - responses c("biomass", "richness", "SP")
apply(region.list.RR[[1]][[2]],1,mean)
#Bahamas: strategy = SP - responses c("biomass", "richness", "SP")
apply(region.list.RR[[1]][[3]],1,mean)

#Haiti: strategy = Biomass - responses c("biomass", "richness", "SP")
apply(region.list.RR[[2]][[1]],1,mean)
#Haiti: strategy = Richenss - responses c("biomass", "richness", "SP")
apply(region.list.RR[[2]][[2]],1,mean)
#Haiti: strategy = SP - responses c("biomass", "richness", "SP")
apply(region.list.RR[[2]][[3]],1,mean)

RRRS<-unlist(region.list.RR)
RRRS<-as.data.frame(RRRS)


quartz(width = 10, height = 4)
par(mfrow = c(1,3))
cors <- c( "#ca0020", "#92c5de", "#0571b0", "#f4a582") #biomass, Rich, SP, PP

for(r in 1:3){
  
  #r=1 #Bahamas
  now <- region.list.RR[[r]]
  
  for(s in 1:4){#strategies - order biomass, richness, sp, pp
    #s=1
    s.now <- now[[s]]
    
    if(s == 1) {
      plot(1:10,
         xlim = c(0,8),
         ylim = c(min(s.now)*.9, max(s.now)*1.1),
         col = "white") }
      
     if(s == 1) {x = c(.5, .75, 1, 1.25) } #biomass
     if(s == 2) {x = c(2.5, 2.75, 3, 3.25) } #richness
     if(s == 3) {x = c(4.5, 4.75, 5, 5.25) } #SP
    if(s == 4) {x = c(6.5, 6.75, 7, 7.25) } #PP

    #biomass
    points(
      jitter(rep(x[1],length(s.now[1,])), factor = .5), 
             s.now[1,], col = cors[1],
      pch = 16, cex = .5)
    points(x[1],mean(s.now[1,]),col = cors[1],
           pch = 16, cex = 2.5)
    #richness
    points(
      jitter(rep(x[2],length(s.now[2,])), factor = 1), 
      s.now[2,], col = cors[2] ,
      pch = 16, cex = .5)
    points(x[2],mean(s.now[2,]),col = cors[2],
           pch = 16, cex = 2.5)
    #SP
    points(
      jitter(rep(x[3],length(s.now[3,])), factor = .5), 
      s.now[3,], col = cors[3],
      pch = 16, cex = .5)
    points(x[3],mean(s.now[3,]),col = cors[3],
           pch = 16, cex = 2.5)
    #PP
    points(
      jitter(rep(x[4],length(s.now[4,])), factor = 1), 
      s.now[4,], col = cors[4],
      pch = 16, cex = .5)
    points(x[4],mean(s.now[4,]),col = cors[4],
           pch = 16, cex = 2.5)
    
  }
  
  
  
  
}



#perc
cors <- c( "#ca0020", "#92c5de", "#0571b0", "#f4a582") #biomass, Rich, SP, PP

for(r in 1:3){
  
  #r=1 #Bahamas
  now <- region.list.Perc[[r]]
  
  for(s in 1:4){#strategies - order biomass, richness, sp, pp
    #s=1
    s.now <- now[[s]]
    
    if(s == 1) {
      plot(1:10,
           xlim = c(0,8),
           ylim = c(min(s.now)*.9, max(s.now)*1.1),
           col = "white") }
    
    if(s == 1) {x = c(.5, .75, 1, 1.25) }
    if(s == 2) {x = c(2.5, 2.75, 3, 3.25) }
    if(s == 3) {x = c(4.5, 4.75, 5, 5.25) }
    if(s == 4) {x = c(6.5, 6.75, 7, 7.25) }
    
    #biomass
    points(
      jitter(rep(x[1],length(s.now[1,])), factor = .5), 
      s.now[1,], col = cors[1],
      pch = 16, cex = .5)
    points(x[1],mean(s.now[1,]),col = cors[1],
           pch = 16, cex = 2.5)
    #richness
    points(
      jitter(rep(x[2],length(s.now[2,])), factor = 1), 
      s.now[2,], col = cors[2] ,
      pch = 16, cex = .5)
    points(x[2],mean(s.now[2,]),col = cors[2],
           pch = 16, cex = 2.5)
    #SP
    points(
      jitter(rep(x[3],length(s.now[3,])), factor = .5), 
      s.now[3,], col = cors[3],
      pch = 16, cex = .5)
    points(x[3],mean(s.now[3,]),col = cors[3],
           pch = 16, cex = 2.5)
    #PP
    points(
      jitter(rep(x[4],length(s.now[4,])), factor = 1), 
      s.now[4,], col = cors[4],
      pch = 16, cex = .5)
    points(x[4],mean(s.now[4,]),col = cors[4],
           pch = 16, cex = 2.5)
    
  }
  
}

###transform list into df for plotting###
region.df.rr <- bind_rows(
  bind_rows(
  as.data.frame(region.list.RR[[1]][1]),
  as.data.frame(region.list.RR[[1]][2]),
  as.data.frame(region.list.RR[[1]][3]),
  as.data.frame(region.list.RR[[1]][4])
  ) %>% mutate(Region = "All"),
  bind_rows(
    as.data.frame(region.list.RR[[2]][1]),
    as.data.frame(region.list.RR[[2]][2]),
    as.data.frame(region.list.RR[[2]][3]),
    as.data.frame(region.list.RR[[2]][4])
  ) %>% mutate(Region = "Bahamas"), 
  bind_rows(
    as.data.frame(region.list.RR[[3]][1]),
    as.data.frame(region.list.RR[[3]][2]),
    as.data.frame(region.list.RR[[3]][3]),
    as.data.frame(region.list.RR[[3]][4])
  ) %>% mutate(Region = "Haiti")
) %>% 
  mutate(Metric = rep(c("Biomass", "Sp richness", "SP", "PP"), 12),
         Cons_Strat = rep(c("Biomass", "Sp richness", "SP", "PP"), each = 4, times = 3)) %>% 
  dplyr::select(Region, Cons_Strat, Metric, everything()) %>% 
  pivot_longer(X1:X100, names_to = "rep_num", values_to = "rep_value") %>% 
  group_by(Region, Cons_Strat, Metric) %>% 
  summarise(mean = mean(rep_value), sd = sd(rep_value))

####final RR figure###
region.df.rr$Cons_Strat<-factor(region.df.rr$Cons_Strat, levels = c("Sp richness", "Biomass", "SP", "PP"))
region.df.rr$Metric<-factor(region.df.rr$Metric, levels = c("Sp richness", "Biomass", "SP", "PP"))

(region.plots<-ggplot(data = region.df.rr, 
                    aes(x = Cons_Strat, y = mean, fill = Metric)) +
    geom_col(position=position_dodge(0.9)) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.3, linewidth = 0.5, inherit.aes = T, 
                 position = position_dodge(0.9)) +
    scale_fill_manual(values = c("#008E98","#A4A4A4","#8B0069", "#F8C149"))+
    facet_wrap(~Region, ncol = 3) +
    geom_hline(yintercept = 0, linetype = "dashed")+
    ylab("Response ratio")+
    xlab("Conservation strategy")+
    theme_bw(base_size = 15)+
    theme(legend.key.size = unit(0.75, 'cm'),legend.spacing.y = unit(5, 'cm'), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          strip.background = element_blank(),axis.text.y = element_text(size = 16),
          strip.text.x = element_blank())
)
