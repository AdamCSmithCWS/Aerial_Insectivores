### Friday am convert this to a function

####
socb.jags = function(
indices = dat,
reg = "Contiental",
group = "Continental",
ind = "index",
lci = "lci",
uci = "uci",
year = "year",
base.yr = bsyr,
writedat = F,
ylm = NA,
hght = 12,
wdth = 18){
# require(dplyr)
require(rjags)
require(ggplot2)
  require(ggrepel)
  require(RColorBrewer)
  require(colorscience)
  

  set.seed(2019)
  
 
  
  
#  if(group %in% c("introduced","Other.birds")){ add12 = F }
  indices <- indices[order(indices$English_Name,indices$year),]
  
   splist <- unique(indices[,c("English_Name","Sort_Order","French_Name","Scientific_Name")])
  # 
  # splist$spfactor <- factor(splist$English_Name,
  #                           levels = splist$English_Name,
  #                           ordered = T)
  
 # splist$spfact <- as.integer(splist$spfactor) 
  
  

yrs <- min(indices[,year]):max(indices[,year])

indices[,"se"] <- ((indices[,uci]-indices[,lci])/(1.96*2))





indices[,"index.s"] <- NA
indices[,"cv"] <- NA
indices[,"cvar.s"] <- NA
indices[,"se.s"] <- NA
indices[,"prec.logthetahat"] <- NA
indices[,"logthetahat"] <- NA

# 
# splist = merge(splist,
#                spslist[,c("English_Name","English_Name_ID","Sort_Order")],
#                by.x = "English_Name",
#                by.y = "English_Name")

base.i <- rep(NA, length = length(unique(indices$English_Name)))
names(base.i) <- unique(indices$English_Name)
base.se.sp <- rep(NA, length = length(unique(indices$English_Name)))
names(base.se.sp) <- unique(indices$English_Name)

# mn.i <- vector(length = length(unique(indices$English_Name)))
# names(mn.i) <- unique(indices$English_Name)


for (sn in 1:nrow(splist)) {
  
  s = splist[sn,"English_Name"]
  
 
  se = "se"
  r <- which(indices[,"English_Name"] == s)
  

  
  
  base.s <- indices[which(indices$English_Name == s & indices$year == base.yr),ind] #stores the base index value 
  base.se <- indices[which(indices$English_Name == s & indices$year == base.yr),se] #stores the base se value 
  if(is.na(base.s)){
    byr <- min(indices[which(indices$English_Name == s & !is.na(indices[,ind])),"year"],na.rm = T)
    base.s <- indices[which(indices$English_Name == s & indices$year == byr),ind] #stores the base index value 
    base.se <- indices[which(indices$English_Name == s & indices$year == byr),se] #stores the base se value 
    
  }
  #mn.i[s]    <- mean(indices[which(indices$English_Name == s),ind],na.rm = T)
  for (y in r) {
    
    indices[y,"index.s"] <- (indices[y,ind])/base.s # standardized index WRT base year
    indices[y,"cv"] <- indices[y,se]/indices[y,ind] #coefficient of variation of the original index
    indices[y,"cvar.s"] <- (((indices[y,se]^2)/(indices[y,ind]^2))+((base.se^2)/(base.s^2))) 
    ### this line sums the ratios of the variance over the index for year y and 
    ### for the base year -
    ### the cvar.s calculation assumes that index estimates are independent among years, which is clearly false, however it is reasonable given that:
    ## there are no estimates of the covariance among years for indices from any of the indices sources, and
    ## the scaling of variances in the final analysis are all relative not absolute and in the absence of more information, assuming that the covariance of annual indices is the same among different English_Name seems reasonable.
    indices[y,"se.s"] <- (sqrt(indices[y,"cvar.s"])*indices[y,"index.s"])
    indices[y,"logthetahat"] <- log(indices[y,"index.s"])      
    indices[y,"prec.logthetahat"] <- 1/log(1+((indices[y,"se.s"]^2)/(indices[y,"index.s"]^2)))

      
  }


    
  print(s)
}




#### insert code to fill all missing English_Name data with 1-values
### adding the assumption that English_Name with missing data have not changed
### rather arbitrary assumption, but avoids having jagged changes in the
### indices as new English_Name are added in.
for(y in base.yr:2017){
  rs <- which(indices[,year] == y)
  # nspecies <- nrow(splist[which(splist$firstyear <= y &
  #                                 splist$English_Name %in% indices[rs,"English_Name"]),])
  # 
  # 
  nspecies = length(rs)
  
  # logthetahat <- rep(NA,nspecies)
  # prec.logthetahat <- logthetahat 
  # 
  logthetahat <- indices[rs,"logthetahat"]
  
  prec.logthetahat <- indices[rs,"prec.logthetahat"]
  
  nmspecies = which(is.na(logthetahat))
  #nspecies = which(!is.na(logthetahat))
  
  if(length(nmspecies) > 0){
    
    print(paste("ERROR species",nmspecies,"missing data"))
    break
  }
  
  data.jags = list(nspecies = nspecies,
                   logthetahat = logthetahat,
                   prec.logthetahat = prec.logthetahat)
  
  # inits = function() {
  #   prec.logtheta = runif(1,1,2)
  #   mu.logtheta = runif(1,0.4,0.5)
  #   logtheta = rep(0,nspecies)
  #   logthetahat = vector(length = nspecies)
  #   logthetahat[mspecies] <- NA
  #   logthetahat[nmspecies] <- 0
  #  
  #   list(prec.logtheta = prec.logtheta, 
  #        mu.logtheta = mu.logtheta,
  #        logtheta = logtheta,
  #        logthetahat = logthetahat,
  #        prec.logthetahat = logthetahat+1)
  # }
  # 
  params = c("logtheta",
             "mu.logtheta",
             "prec.logtheta",
             "expmu")
  
  mod = "c:/SOCB/model_logtheta.txt"
  
  
  adaptSteps = 500              # Number of steps to "tune" the samplers.
  burnInSteps = 20000            # Number of steps to "burn-in" the samplers.
  nChains = 3                   # Number of chains to run.
  numSavedSteps=10000           # Total number of steps to save.
  thinSteps=10                   # Number of steps to "thin" (1=keep every step).
  nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
  
  t1 = Sys.time()
  
  
  jagsMod = jags.model( mod, 
                        data= data.jags ,  
                        n.chains= nChains , 
                        n.adapt= adaptSteps)
  
  # Burn-in:
  cat( paste(y,"Burning in the MCMC chain...\n" ))
  update( jagsMod , n.iter=burnInSteps )
  # The saved MCMC chain:
  cd13 = coda.samples( jagsMod , variable.names=params ,
                       n.iter=nIter , thin=thinSteps )
  
  sums <- summary(cd13)
  expmut <- as.data.frame(t(sums$quantiles["expmu",]))
  expmut[,"year"] <- y
  expmut[,"nspecies"] <- nspecies
  expmut[,"species.list"] <- paste(splist[which(splist$English_Name %in% indices[rs,"English_Name"]),"English_Name"],
                                   collapse = " - ")
  
  logtht <- as.data.frame(sums$quantiles[paste0("logtheta[",
                                                1:nspecies,
                                                "]"),])
  logtht[,"node"] <- row.names(logtht)
  
  logtht[,"year"] <- y
  
  logtht[,"spfact"] <- 1:nspecies
  logtht[,"spfactout"] <- 1:nspecies
  splist[,"spfactsplist"] <- as.integer(splist$English_Name)
  logtht <- merge(logtht,splist,
                  by.x = "spfact",
                  by.y = "spfactsplist")
  
  if(y == base.yr){
    expmu <- expmut
    logth <- logtht
  }else{
    expmu <- rbind(expmu,expmut)
    logth <- rbind(logth,logtht)
  }
  


}### y


expmuout = expmu
expmuout$Indicator = ((expmuout[,"50%"])-1)*100
expmuout$Indicator.lci = ((expmuout[,"2.5%"])-1)*100
expmuout$Indicator.uci = ((expmuout[,"97.5%"])-1)*100

expmuout = expmuout[,c("year","Indicator","Indicator.lci","Indicator.uci","nspecies","species.list")]
 
 sp.col = brewer.pal(8,"Set1")[-6]

    graphname = paste0("output/",reg,"/",group," SOCB index")

  ib = indices[which(indices$year == base.yr),]
  ie = indices[which(indices$year == 2017),"index.s"]
  
  # nlargedec1 = length(which(ie < 0.5))
  # nmoddec1 = length(which(ie < 0.75 & ie > 0.5))
  # ndec = length(which(ie < 1 & ie > 0.75))
  # sumdec = sum(nlargedec1,nmoddec1,ndec)
  # nlargedec = sumdec
  # nmoddec = sum(ndec,nmoddec1)
  # 
  # 
  # nlargeinc1 = length(which(ie > 2))
  # nmodinc1 = length(which(ie > 1.33 & ie < 2))
  # ninc = length(which(ie > 1 & ie < 1.33))
  # suminc = sum(nlargeinc1,nmodinc1,ninc)
  # nlargeinc = suminc
  # nmodinc = sum(ninc,nmodinc1)
  # 
  # 
  # cmykmat <- read.csv("ColourScale.csv")
  # cmykmat <- cmykmat/100
  # 
  # cmymat <- CMYK2CMY(cmykmat)
  # rgbmat <- CMY2RGB(cmymat)
  # 
  # colset <- rgb(red = pmin(255,rgbmat[,"R"]),
  #               green = pmin(255,rgbmat[,"G"]),
  #               blue = pmin(255,rgbmat[,"B"]),
  #               maxColorValue = 255)
  # 
  # 
  # 
  # nplots = list(nlargedec = nlargedec*-1,
  #               nlargeinc = nlargeinc,
  #               nmoddec = nmoddec*-1,
  #               nmodinc = nmodinc,
  #               ndec = ndec*-1,
  #               ninc = ninc)
  # 
  # nupdown = c(nlargedec1,nlargeinc1,nmoddec1,nmodinc1,ndec,ninc)
  # colsnplots = colset[c(16,1,9,2,6,5)]
  # pdf(paste(graphname,"stacked bar.pdf"),
  #     height = 8,
  #     width = 4)
  # yup = max(c(suminc,sumdec))
  # 
  # plot(-99,
  #      xlim = c(0,1),
  #      ylim = c(yup*-1,yup),
  #      ylab = "",
  #      xlab = "",
  #      yaxt = "n",
  #      bty = "n",
  #      xaxt = "n")
  # polyx = c(0.2,0.4,0.4,0.2)
  # i = 0
  # for(j in names(nplots)){
  #   i = i+1
  #   y1 = c(rep(nplots[[j]],2),0,0)
  #   
  #   polygon(x = polyx,
  #           y = y1,
  #           col = colsnplots[i],
  #           border = NA)
  #   
  # }
  # abline(h = 0)
  # axis(side = 2,
  #      at = seq(yup*-1,yup,by = max(c(1,floor(yup/6)))),
  #      labels = abs(seq(yup*-1,yup,by = max(c(1,floor(yup/6))))),
  #      las = 2)
  # mtext(paste(c("# English_Name increased since","# English_Name decreased since"),base.yr),
  #       at = c(yup/2,(yup/2)*-1),
  #       side = 2,
  #       line = 3)
  # lablsy = c(nplots[["nlargedec"]]+(nlargedec1/2),
  #            nplots[["nlargeinc"]]-(nlargeinc1/2),
  #            nplots[["nmoddec"]]+(nmoddec1/2),
  #            nplots[["nmodinc"]]-(nmodinc1/2),
  #            nplots[["ndec"]]+(ndec/2),
  #            nplots[["ninc"]]-(ninc/2)
  # )
  # text(x = 0.41,
  #      y = lablsy[which(nupdown != 0)],
  #      c(       "large decrease",
  #               "large increase",
  #               "moderate decrease",
  #               "moderate increase",
  #               "small decrease",
  #               "small increase"
  #      )[which(nupdown != 0)],
  #      pos = 4)
  # 
  # dev.off()
  # 
  # 
  # ############### end stacked bar
  # 

  
  if(base.yr == 1970){xlm = c(1970,2017)}else{xlm = c(base.yr,2017)}
  
  for(addsp in c(1,2)){
    addsplines = T
    if(addsp == 1) {
      pdf(paste0(graphname,"present w sp.pdf"),
          height = hght,
          width = wdth)
      par(mar = c(4,5,1,9),
          xpd = T)
      
    }else{
      pdf(paste0(graphname,"present.pdf"),
          height = hght,
          width = wdth,
          pointsize = 30)
      par(mar = c(4,5,1,9),
          xpd = T)
      
      
    }
    # par(mar = c(4,5,1,1))
  if(any(is.na(ylm))){
    if(base.yr == 2006){
      
      ylm <- c(0.5,2)
    }else{
      ylm <- c(min(indices[,"index.s"],na.rm = T),max(indices[,"index.s"],na.rm = T))
      
    }
  }
    xplot = c(base.yr:2017)
    
    plot(1,1,
         xlim = xlm,
         ylim = ylm,
         ylab = "",
         xlab = "",
         log = "y",
         bty = "l",
         xaxs = "i",
         yaxt = "n",
         xaxt = "n",
         yaxs = "i")
    
    ylabs = c(-95,-90,-85,-80,-75,-67,-50,-33,-20,0,25,50,100,200,400,600,800)
    ylabsy = (ylabs/100)+1
    ylabs = ylabs[which(ylabsy > min(ylm) & ylabsy < max(ylm))]
    ylabsy = ylabsy[which(ylabsy > min(ylm) & ylabsy < max(ylm))]
    
if(length(ylabs) > 9){
  ylabs = c(-95,-90,-80,-67,-50,0,100,400,800)
  ylabsy = (ylabs/100)+1
  ylabs = ylabs[which(ylabsy > min(ylm) & ylabsy < max(ylm))]
  ylabsy = ylabsy[which(ylabsy > min(ylm) & ylabsy < max(ylm))]
  
  
}
    xlabs = c(seq(1970,2010,by = 10),2017)
    xlabs = xlabs[which(xlabs > base.yr)]
    xlabs = c(base.yr,xlabs)
    axis(side = 2,
         at = ylabsy,
         labels = paste0(ylabs,"%"),
         las = 2)
    axis(side = 1,
         at = xlabs)
    mtext(side = 2,
          at = median(ylabsy),
          paste("Average % change since",base.yr),
          line = 4)
    if(addsp == 1){
    polygon(x = c(xplot,rev(xplot)),
            y = c(expmu[,"2.5%"],rev(expmu[,"97.5%"])),
            col = transp.func(grey(0.8),0.2),
            border = NA)
    }
    lines(x = xplot,
          y = expmu[,"50%"],
          lwd = 3)
    lines(x = xlm,
          y = c(1,1),
          lty = 3,
          lwd = 2)

    if(addsp == 1){
      j = 0
      for(s in splist$English_Name){
        j = j+1
        tmp <- indices[which(indices$English_Name == s &
                               indices$year >= min(xlm)),]
        tmpo = tmp[,c("year","logthetahat","prec.logthetahat")]
        tmpo$prec.logthetahat = 1/(sqrt(tmpo$prec.logthetahat))
        names(tmpo) = c("year",s,paste0(s,".se"))
        
        tmpo2 = tmpo
        tmpo2[,paste0(s,".index")] = (exp(tmpo2[,s])-1)*100
        tmpo2[,paste0(s,".lci")] = (exp(tmpo2[,s]-(1.96*tmpo2[,paste0(s,".se")]))-1)*100
        tmpo2[,paste0(s,".uci")] = (exp(tmpo2[,s]+(1.96*tmpo2[,paste0(s,".se")]))-1)*100
        tmpo2 = tmpo2[,c("year",paste0(s,c(".index",".lci",".uci")))]
        if(j == 1){
          alldatout = merge(expmu,
                            tmpo,
                            by = "year",
                            all.x = T)
          alldatout2 = merge(expmuout,
                            tmpo2,
                            by = "year",
                            all.x = T)
          alldatout3 = tmp
    
        }else{
          alldatout = merge(alldatout,
                            tmpo,
                            by = "year",
                            all.x = T)
          alldatout2 = merge(alldatout2,
                            tmpo2,
                            by = "year",
                            all.x = T)

          alldatout3 = rbind(alldatout3,
                             tmp)
        }
        lines(x = tmp$year,
              y = tmp[,"index.s"],
              lwd = 1,#2,#1,
              col = transp.func(rep(sp.col,length = nrow(splist))[j],0.4))
        
        if(tmp[nrow(tmp),"index.s"] > max(ylm) | 
           tmp[nrow(tmp),"index.s"] < min(ylm)){
          
          if(tmp[nrow(tmp),"index.s"] > max(ylm)){
            w = max(which(tmp$index.s < max(ylm)))
            text(s,
                 x = tmp[w,"year"],
                 y = tmp[w,"index.s"],
                 pos = 4,
                 col = transp.func(rep(sp.col,length = nrow(splist))[j],0.6),
                 cex = 1)#1)#
          }else{
            w = max(which(tmp$index.s > min(ylm)))
            text(s,
                 x = tmp[w,"year"],
                 y = tmp[w,"index.s"],
                 pos = 4,
                 col = transp.func(rep(sp.col,length = nrow(splist))[j],0.4),
                 cex = 1)#1)#
          }
        }else{
          text(s,
               x = 2017,
               y = tmp[nrow(tmp),"index.s"],
               pos = 4,
               col = transp.func(rep(sp.col,length = nrow(splist))[j],0.4),
               cex = 1)#1)#
        }
        
      }
    }else{
      
      j = 0
      for(s in splist$English_Name){
        j = j+1
        tmp <- indices[which(indices$English_Name == s &
                               indices$year >= min(xlm)),]
        # tmpo = tmp[,c("year","logthetahat","prec.logthetahat")]
        # tmpo$prec.logthetahat = 1/(sqrt(tmpo$prec.logthetahat))
        # names(tmpo) = c("year",s,paste0(s,".se"))
        # 
        # tmpo2 = tmpo
        # tmpo2[,paste0(s,".index")] = (exp(tmpo2[,s])-1)*100
        # tmpo2[,paste0(s,".lci")] = (exp(tmpo2[,s]-(1.96*tmpo2[,paste0(s,".se")]))-1)*100
        # tmpo2[,paste0(s,".uci")] = (exp(tmpo2[,s]+(1.96*tmpo2[,paste0(s,".se")]))-1)*100
        # tmpo2 = tmpo2[,c("year",paste0(s,c(".index",".lci",".uci")))]
        # # if(j == 1){
        #   alldatout = merge(expmu,
        #                     tmpo,
        #                     by = "year",
        #                     all.x = T)
        #   alldatout2 = merge(expmuout,
        #                      tmpo2,
        #                      by = "year",
        #                      all.x = T)
        #   alldatout3 = tmp
        #   
        # }else{
        #   alldatout = merge(alldatout,
        #                     tmpo,
        #                     by = "year",
        #                     all.x = T)
        #   alldatout2 = merge(alldatout2,
        #                      tmpo2,
        #                      by = "year",
        #                      all.x = T)
        #   
        #   alldatout3 = rbind(alldatout3,
        #                      tmp)
        # }
        lines(x = tmp$year,
              y = tmp[,"index.s"],
              lwd = 4,#1,
              col = transp.func(rep(grey(0.5),length = nrow(splist))[j],0.4))
        
        # if(tmp[nrow(tmp),"index.s"] > max(ylm) | 
        #    tmp[nrow(tmp),"index.s"] < min(ylm)){
        #   
        #   if(tmp[nrow(tmp),"index.s"] > max(ylm)){
        #     w = max(which(tmp$index.s < max(ylm)))
        #     text(s,
        #          x = tmp[w,"year"],
        #          y = tmp[w,"index.s"],
        #          pos = 4,
        #          col = transp.func(rep(sp.col,length = nrow(splist))[j],0.6),
        #          cex = 1)#1)#
        #   }else{
        #     w = max(which(tmp$index.s > min(ylm)))
        #     text(s,
        #          x = tmp[w,"year"],
        #          y = tmp[w,"index.s"],
        #          pos = 4,
        #          col = transp.func(rep(sp.col,length = nrow(splist))[j],0.4),
        #          cex = 1)#1)#
        #   }
        # }else{
        #   text(s,
        #        x = 2017,
        #        y = tmp[nrow(tmp),"index.s"],
        #        pos = 4,
        #        col = transp.func(rep(sp.col,length = nrow(splist))[j],0.4),
        #        cex = 1)#1)#
        # }
        
      }
      # if(add12){
      #   lines(x = oldind[which(!is.na(oldind[,group])),"Year"],
      #         y = oldind[which(!is.na(oldind[,group])),group],
      #         lwd = 2,
      #         col = grey(0.3),
      #         lty = 2)
      # }
      # polygon(x = c(xplot,rev(xplot)),
      #         y = c(expmu[,"2.5%"],rev(expmu[,"97.5%"])),
      #         col = transp.func(grey(0.8),0.2),
      #         border = NA)
      lines(x = xplot,
            y = expmu[,"50%"],
            lwd = 3)
    }
    dev.off()
  }#end addlines loop
  
  expmup = expmu
  names(expmup)[1:5] = c("lci","lqrt","med","uqrt","uci")
  alldatout3$English_Name = as.character(alldatout3$English_Name)
  #### spaghetti plots
  pdf(paste0(graphname,"spaghetti w sp.pdf"),
      height = hght,
      width = wdth)
  par(mar = c(4,5,1,1),
      xpd = T)
  labs = alldatout3[which(alldatout3$year == max(alldatout3$year)),]
  ylabla = c(-90,-80,-67,-50,-33)
  ylabbrk = c(((ylabla/100)+1),1,rev(1/((ylabla/100)+1)))
  ylabl = (ylabbrk-1)*100
  ylabl[which(ylabl > 100)] = round(ylabl[which(ylabl > 100)]/100)*100
  ylabl[which(ylabl>0)] = ceiling(ylabl[which(ylabl>0)])
  ylabl = paste0(ylabl,"%")
  
  hl = data.frame(year = range(alldatout3$year),index.s = c(1,1))
  p2 = ggplot(data = alldatout3,aes(x = year,y = index.s))+
    theme_classic()+
    theme(legend.position = "none",
          axis.text = element_text(size = 22))+
    scale_y_continuous(trans = "log",breaks = ylabbrk,labels = ylabl,name = NULL)+
    scale_x_continuous(limits = c(base.yr,2030),breaks = c(seq(1970,2010,by = 5),2017),expand = expand_scale(mult = 0, add = 0),name = NULL)+
    geom_ribbon(data = expmup,aes(x = year,ymin = lci,ymax = uci),inherit.aes = F,fill = grey(0.7),alpha = 0.5)+
    #geom_line(data = alldatout3,aes(group = English_Name),colour = grey(0.8))+
    geom_line(data = alldatout3,aes(group = English_Name,colour = English_Name))+
    geom_line(data = expmup,aes(x = year,y = med))+
    geom_line(data = hl,aes(x = year,y = index.s),linetype = 3,colour = grey(0.2),inherit.aes = F)+
    geom_text_repel(data = labs,aes(x = year,y = index.s,label = English_Name,colour = English_Name),size = 6,xlim = c(2015,2030))
    
  print(p2)
  dev.off()#spaghetti end
  
 
  
  
  
  
  
   
  for(cc in c("2.5%","25%","50%","75%","97.5%")){
    expmu[,paste("percent change",cc)] = (expmu[,cc]-1)*100
  }
  
  
  write.csv(expmu[which(expmu$year >= base.yr),],paste(graphname,"values.csv"))
  write.csv(alldatout,paste(graphname,"values incl English_Name.csv"),row.names = F)
  write.csv(alldatout2,paste(graphname,"values incl English_Name for distribution.csv"),row.names = F)
  
  return(list(alldatout = alldatout,
              graphname = graphname,
              indices = indices))
  
}




  
  
  
  
  
  
  


