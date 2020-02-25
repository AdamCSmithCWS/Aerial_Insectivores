#### jags function for SOCB model
setwd("c:/SOCB/grassland/")

pkgs = c("sp","maptools","rgdal","maptools","RColorBrewer","GISTools","ggplot2",
         "dplyr","tidyr","scales","geofacet","stringr","tools",
         "ggrepel","e1071","rjags")
inst.pck = installed.packages()
inst.pck2 = names(inst.pck[,2])

if(length(pkgs[-which(pkgs %in% inst.pck2)])>0){
  install.packages(pkgs[-which(pkgs %in% inst.pck2)])
}

for(i in 1:length(pkgs)){
  library(pkgs[i],character.only = T)
}

####add the last few years of the 1970 index to the 2006 plots


spslist = read.csv(paste0("final survey and species lists SOCB 2018.csv"),
                   stringsAsFactors = F)


source("c:/Functions/transparency function.r")

spgrps = read.csv("spslists.csv",stringsAsFactors = F)
grpnms = names(spgrps)[-1]

basedir = getwd()
for(gr in grpnms){
setwd(basedir)
alldata = read.csv(paste0("BBS 2017 GAM all indices.csv"))

alldata = alldata[which(alldata$year > 1969),]

sp.inc = spgrps[which(spgrps[,gr]),"species"] #species selection
alldata = alldata[which(alldata$species %in% sp.inc),] #

outdir = paste0(basedir,"/",gr,"/")
dir.create(outdir)
setwd(outdir)

dir.create(paste0(outdir,"output"))

dir.create(paste0(outdir,"output/State"))
dir.create(paste0(outdir,"output/Strat.name"))
dir.create(paste0(outdir,"output/countrynum"))
dir.create(paste0(outdir,"output/Continental"))

nybysp = table(alldata$species)
nspeciestotal = length(nybysp)
fybysp = tapply(alldata$year,alldata$species,min)


# alldata$year = alldata[,yearcol]
 alldata$index = alldata[,"med"]
 # alldata$se = (alldata[,"lci"]-alldata[,"uci"])/(1.96*4)

# alldata$lci = alldata[,lcicol]
# alldata$uci = alldata[,ucicol]
# alldata[which(alldata$index == 0),"index"] = 0.001 
# alldata[which(alldata$lci == 0),"lci"] = 0.0005 
# alldata[which(alldata$uci == 0),"uci"] = 0.0015 

today = Sys.Date()
alldataout = alldata[,c("species",
                        "year",
                        "index",
                        "lci",
                        "uci",
                        "region",
                        "type",
                        "o.mean")]
names(alldataout)[c(1)] <- c("English_Name")

alldataout = merge(alldataout,spslist[,c("English_Name","French_Name","Sort_Order","Scientific_Name")], by = "English_Name")
alldataout = alldataout[order(alldataout$Sort_Order,alldataout$type,alldataout$region,alldataout$year),]

alldataout$English_Name = factor(alldataout$English_Name,ordered = T)

#regs = unique(alldataout[which(alldataout$type %in% c("Continental","countrynum")),"region"])
regs = unique(alldataout[,"region"])

source(paste0(basedir,"/socb.jags.func.grass.r"))
allResultExport = list()
length(allResultExport) = length(regs)


baseline = 1970
## loop through both the long-term and short-term baselines
  q = 0
for(wind in regs){

    q = q+1

    pg = unique(alldataout[which(alldataout$region == wind),"type"])
    
  names(allResultExport)[q] = wind
  
    spsel = unique(alldataout$English_Name)

    dat = alldataout[which(alldataout$region == wind),]
    
  spmiss = spsel[-which(spsel %in% dat$species)]
  
  nspy = table(dat$year)
  psp = nspy/max(nspy)
  if(max(nspy) < 5 | mean(psp) < 0.5){next}
 
  fyr = max(c(min(as.integer(names(nspy)[which(nspy > 1)])),min(as.integer(names(psp)[which(psp > 0.5)])))) 
  # if(length(spmiss) > 0){
  #   print(paste("DATA missing for",spmiss,"in",ind))
  # }
  # 
  bsyr = baseline

   dat = dat[which(dat$year >= bsyr),] 
   
   
   
   
   
   
   
   
   
   
  
  allResultExport[[wind]] = socb.jags(indices = dat,
            reg = pg,
            group = wind,
            ind = "index",
            lci = "lci",
            uci = "uci",
            year = "year",
            base.yr = bsyr)
  
  
  
  
  
}
  
 
  save(allResultExport,
       file = paste("SCB results GAM",baseline,".RData"))
  
   
  
  
}#grpnms end modeling loop







### begin plotting loops

for(gr in grpnms){
  setwd(basedir)
  alldata = read.csv(paste0("BBS 2017 GAM all indices.csv"))
  
  alldata = alldata[which(alldata$year > 1969),]
  
  sp.inc = spgrps[which(spgrps[,gr]),"species"] #species selection
  alldata = alldata[which(alldata$species %in% sp.inc),] #
  
  outdir = paste0(basedir,"/",gr,"/")
  setwd(outdir)
  
  
  nybysp = table(alldata$species)
  nspeciestotal = length(nybysp)
  fybysp = tapply(alldata$year,alldata$species,min)
  
  
  # alldata$year = alldata[,yearcol]
  alldata$index = alldata[,"med"]
  # alldata$se = (alldata[,"lci"]-alldata[,"uci"])/(1.96*4)
  
  # alldata$lci = alldata[,lcicol]
  # alldata$uci = alldata[,ucicol]
  # alldata[which(alldata$index == 0),"index"] = 0.001 
  # alldata[which(alldata$lci == 0),"lci"] = 0.0005 
  # alldata[which(alldata$uci == 0),"uci"] = 0.0015 
  
  today = Sys.Date()
  alldataout = alldata[,c("species",
                          "year",
                          "index",
                          "lci",
                          "uci",
                          "region",
                          "type",
                          "o.mean")]
  names(alldataout)[c(1)] <- c("English_Name")
  
  alldataout = merge(alldataout,spslist[,c("English_Name","French_Name","Sort_Order","Scientific_Name")], by = "English_Name")
  alldataout = alldataout[order(alldataout$Sort_Order,alldataout$type,alldataout$region,alldataout$year),]
  
  alldataout$English_Name = factor(alldataout$English_Name,ordered = T)
  
  #regs = unique(alldataout[which(alldataout$type %in% c("Continental","countrynum")),"region"])
  regs = unique(alldataout[,"region"])
  baseline = 1970
  
  
  load(paste("SCB results GAM",baseline,".RData"))
 


  
  ############# facet plot for strata
  ############# facet plot for strata
  ############# facet plot for strata
  ############# facet plot for strata
  ############# facet plot for strata
  ############# facet plot for strata
  rc = 0
for(rr in unique(alldataout[which(alldataout$type == "strat.name"),"region"])){

  tmp = allResultExport[[rr]]$alldatout
  if(length(tmp) < 1){next}
  rc = rc+1
  tmp$region = rr
  tmp$state = strsplit(rr,fixed = T,split = "-")[[1]][1]
  tmp$BCR = strsplit(rr,fixed = T,split = "-")[[1]][2]
  names(tmp)[c(2,4,6,7,8)] = c("lci","index","uci","nspecies","splist")
  tmp = tmp[,c("year","lci","index","uci","nspecies","region","state","BCR","splist")]
  if(rc == 1){
    tpout = tmp
  }else{
    tpout = rbind(tpout,tmp)
  }
  #### compile all regions predictions into a single data frame with state names as well as strata names
}
  
  
  stprovfacet = read.csv(paste0(basedir,"/BBB_StateProvCWS_facet_grid2.csv"),stringsAsFactors = F)
  #stprovfacet = stprovfacet[-which(stprovfacet$code == "HI"),]
  ### trajectory plots by prov state
  stprovfacet[nrow(stprovfacet)+1,"row"] = 1
  stprovfacet[nrow(stprovfacet),"col"] = 5
  stprovfacet[nrow(stprovfacet),"code"] = "BCR7"
  stprovfacet[nrow(stprovfacet),"name"] = "Bcr7" 
  
  country = rep("USA",nrow(stprovfacet))
  country[which(stprovfacet$code %in% c("BC",
                                        "YT",
                                        "NT",
                                        "NU",
                                        "BCR7",
                                        "NSPE",
                                        "PE",
                                        "NS",
                                        "ON",
                                        "AB",
                                        "SK",
                                        "MB",
                                        "QC",
                                        "NB",
                                        "NL"))] = "CAN"
  strcol = rep("darkslategray1",length(country))
  strcol[which(stprovfacet$country == "CAN")] = "lightpink"
  

    # indt$prst = gsub(str_extract(indt$stratcode,"-(.*)-"),pattern = "-",replacement = "",fixed = T)
    # indt$prst = as.character(factor(indt$prst,levels = unique(stprovfacet$code)))
  indt = tpout
    indt$state = toTitleCase(tolower(indt$state))
    indt$prst = as.character(factor(indt$state,levels = unique(stprovfacet$name)))
    labsbcr = indt[which(indt$year == max(indt$year,na.rm = T)),]
    
    cols = function(x){
      cl = "stable"
      if(x > 1){cl = "inc"}
      if(x < 0.5 & x > 0.25){cl = "dec"}
      if(x < 0.25){cl = "Ldec"}
      
      return(cl)
    }
    
    for(pp in labsbcr$region){
      indt[which(indt$region == pp),"col"] = cols(labsbcr[which(labsbcr$region == pp),"index"])
    }
    
    
    
     ptraj <- ggplot(data = indt) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(colour = grey(0.2)),
            strip.background = element_rect(fill = grey(0.97)),#strcol #, colour = grey(0.9), size = NULL, linetype = NULL, color = NULL, inherit.blank = FALSE
            #axis.line = element_line(colour = "black"),
            legend.position = "none") +
      labs(title = paste("Grassland birds indicator by strata (separate BCRS) within Provinces and States"), x = "", y = "mean % population change since 1970") +

      #geom_pointrange(data = indt, aes(x = year, y = mean,ymin = lci,ymax = uci),alpha = 0.3) +
       geom_line(data = indt, aes(x = year, y = index,group = region,colour = col)) +#, colour = stratcode)) +
       geom_ribbon(data = indt, aes(x = year, ymin = lci, ymax = uci,group = region,fill = col),alpha = 0.12)+ #, fill = stratcode , fill = grey(0.8)
       geom_text_repel(data = labsbcr, aes(x = year, y = index,label = BCR),
                       size = 2)+
       #geom_hline(yintercept = 0,colour = gray(0.8))+
       facet_geo(facets = ~ prst,grid = stprovfacet,scales = "free_y", label = "code")+
       coord_cartesian(ylim = c(0.15,1.75))+
       scale_y_continuous(trans = "log", breaks = ((c(-75,-50,0)/100)+1))+ #,labels = c("-75","-50","0")
    scale_x_continuous(breaks = c(1970,2017),minor_breaks = c(1980,1990,2000,2010)) #limits = c(1966, 2017), oob = squish, 
  pdf(file = paste0(outdir,"Grassland indicator geofacet plot by stratum.pdf"),
      height = 8.5,width = 11)
    print(ptraj)
    dev.off()
      # for(p in stprovfacet$code)){
    #
    # }#p


    
    
    
    
    ############# by state-prov
    rc = 0
    for(rr in unique(alldataout[which(alldataout$type == "State"),"region"])){
      
      tmp = allResultExport[[rr]]$alldatout
      if(length(tmp) < 1){next}
      rc = rc+1
      tmp$region = rr
      tmp$state = strsplit(rr,fixed = T,split = "-")[[1]][1]
      tmp$BCR = strsplit(rr,fixed = T,split = "-")[[1]][2]
      names(tmp)[c(2,4,6,7,8)] = c("lci","index","uci","nspecies","splist")
      tmp = tmp[,c("year","lci","index","uci","nspecies","region","state","BCR","splist")]
      if(rc == 1){
        tpout = tmp
      }else{
        tpout = rbind(tpout,tmp)
      }
      #### compile all regions predictions into a single data frame with state names as well as strata names
    }
    
    
    indt = tpout
    indt$state = toTitleCase(tolower(indt$state))
    indt$prst = as.character(factor(indt$state,levels = unique(stprovfacet$name)))
    labsbcr = indt[which(indt$year == max(indt$year,na.rm = T)),]

    for(pp in labsbcr$region){
      indt[which(indt$region == pp),"col"] = cols(labsbcr[which(labsbcr$region == pp),"index"])
    }
    
    
    indt$col = factor(indt$col,levels = c("Ldec","dec","stable","inc"),ordered = T)
    
    
    myColors <- brewer.pal(10,"RdYlBu")[c(1,3,4,10)]
    names(myColors) <- (levels(indt$col))
    
    
    colScale <- scale_colour_manual(values = myColors, aesthetics = c("colour","fill"))
    
    
    
    
    ptraj <- ggplot(data = indt) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(colour = grey(0.2)),
            strip.background = element_rect(fill = grey(0.97)),#strcol #, colour = grey(0.9), size = NULL, linetype = NULL, color = NULL, inherit.blank = FALSE
            #axis.line = element_line(colour = "black"),
            legend.position = "none") +
      labs(title = paste("Grassland birds indicator by Provinces and States"), x = "", y = "mean % population change since 1970") +
      
      #geom_pointrange(data = indt, aes(x = year, y = mean,ymin = lci,ymax = uci),alpha = 0.3) +
      geom_line(data = indt, aes(x = year, y = index,group = state,colour = col)) +
      geom_ribbon(data = indt, aes(x = year, ymin = lci, ymax = uci,group = region,fill = col),alpha = 0.12)+ #, fill = stratcode , fill = grey(0.8)
      colScale + #, colour = stratcode)) +
      
      
      #geom_pointrange(data = indt, aes(x = year, y = mean,ymin = lci,ymax = uci),alpha = 0.3) +
      geom_line(data = indt, aes(x = year, y = index,group = region,colour = col)) +#, colour = stratcode)) +
     #geom_text_repel(data = labsbcr, aes(x = year, y = index,label = BCR),
      #                size = 2)+
      #geom_hline(yintercept = 0,colour = gray(0.8))+
      coord_cartesian(ylim = c(0.15,1.75))+
      scale_y_continuous(trans = "log", breaks = ((c(-75,-50,0)/100)+1))+ #,labels = c("-75","-50","0"))+
      facet_geo(facets = ~ prst,grid = stprovfacet,scales = "free_y", label = "code")+
      scale_x_continuous(breaks = c(1970,2017),minor_breaks = c(1980,1990,2000,2010)) #limits = c(1966, 2017), oob = squish, 
    pdf(file = paste0(outdir,"Grassland indicator geofacet plot by stateprov.pdf"),
        height = 8.5,width = 11)
    print(ptraj)
    dev.off()



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ############# commmunity composition plots
   # if(gr == "all"){
    
    
    
    
    
    
    
    
    
    ############# species-level trajectory plots.
    
    
    rc = 0
    for(rr in unique(alldataout[which(alldataout$type == "State"),"region"])){
      
      tmp1 = allResultExport[[rr]]$alldatout
      if(length(tmp1) < 1){next}
      clx = grep(names(tmp1),pattern = ".se",fixed = T)
      spsl = gsub(names(tmp1)[clx],pattern = ".se",replacement = "",fixed = T)
      for(ss in spsl){
        tmp2 = tmp1[,c("year",ss,paste0(ss,".se"))]
        names(tmp2) = c("year","index","se")
        tmp2$species = ss
        if(ss == spsl[1]){
          
          tmp = tmp2
        }else{
          tmp = rbind(tmp,tmp2)
        }
      }
      
      tmp$lci = tmp$index-(tmp$se)*1.96
      tmp$uci = tmp$index+(tmp$se)*1.96
      
      tmp$index = exp(tmp$index)
      tmp$lci = exp(tmp$lci)
      tmp$uci = exp(tmp$uci)
      
      rc = rc+1
      tmp$region = rr
      tmp$state = strsplit(rr,fixed = T,split = "-")[[1]][1]
      tmp$BCR = strsplit(rr,fixed = T,split = "-")[[1]][2]
      if(rc == 1){
        tpout = tmp
      }else{
        tpout = rbind(tpout,tmp)
      }
      #### compile all regions predictions into a single data frame with state names as well as strata names
    }
    
    
    indt1 = tpout
    indt1$state = toTitleCase(tolower(indt1$state))
    indt1$prst = as.character(factor(indt1$state,levels = unique(stprovfacet$name)))
  indt1 = indt1[order(indt1$prst,indt1$year,indt1$species),]
 

    pdf(file = paste0(outdir,"Grassland species lines geofacet plot by stateprov.pdf"),
        height = 8.5,width = 11)

    for(ss in unique(indt1$species)){
      indt = indt1[which(indt1$species == ss),]
    
      
      labsbcr = indt[which(indt$year == max(indt$year,na.rm = T)),]
      labsbcr = merge(labsbcr,spslist[,c("English_Name","Species_ID")],by.x = "species",by.y = "English_Name")
      

      for(pp in labsbcr$region){
        indt[which(indt$region == pp),"col"] = cols(labsbcr[which(labsbcr$region == pp),"index"])
      }
      
      
      indt$col = factor(indt$col,levels = c("Ldec","dec","stable","inc"),ordered = T)
      
      
      myColors <- brewer.pal(10,"RdYlBu")[c(1,3,4,10)]
      names(myColors) <- (levels(indt$col))
      
      
      colScale <- scale_colour_manual(values = myColors, aesthetics = c("colour","fill"))
      
      
      
      
    ptraj <- ggplot(data = indt) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(colour = grey(0.2)),
            strip.background = element_rect(fill = grey(0.97)),#strcol #, colour = grey(0.9), size = NULL, linetype = NULL, color = NULL, inherit.blank = FALSE
            #axis.line = element_line(colour = "black"),
            legend.position = "none") +
      labs(title = paste(ss,"Population trajectories by Provinces and States"), x = "", y = "mean % population change since 1970") +
      
      #geom_pointrange(data = indt, aes(x = year, y = mean,ymin = lci,ymax = uci),alpha = 0.3) +
      geom_line(data = indt, aes(x = year, y = index,group = state,colour = col)) +
      #geom_point(data = indt, aes(x = year, y = index,group = state,colour = species)) +#, colour = stratcode)) +
      
      
      geom_ribbon(data = indt, aes(x = year, ymin = lci, ymax = uci,group = state,fill = col),alpha = 0.12)+ #, fill = stratcode , fill = grey(0.8)
      colScale + #, colour = stratcode)) +
      # geom_text_repel(data = labsbcr, aes(x = year, y = index,label = Species_ID,colour = Species_ID),
      #                 size = 1.5)+
      #geom_hline(yintercept = 0,colour = gray(0.8))+
      coord_cartesian(ylim = c(0.05,4))+
      scale_y_continuous(trans = "log", breaks = ((c(-75,-50,0)/100)+1))+ #,labels = c("-75","-50","0"))+
      facet_geo(facets = ~ prst,grid = stprovfacet,scales = "free_y", label = "code")+
      scale_x_continuous(breaks = c(1970,2017),minor_breaks = c(1980,1990,2000,2010)) #limits = c(1966, 2017), oob = squish, 
    # pdf(file = paste0(ss,"Grassland species lines geofacet plot by stateprov.pdf"),
    #     height = 8.5,width = 11)
     print(ptraj)
    # dev.off()
    # 
    }
  dev.off()
  
  #  }##3 end if gr == all
    
    
    
    
}#end gr loop


#### combine indicators into plots that look like the 2012 plots
#### combine indicators into plots that look like the 2012 plots
#### combine indicators into plots that look like the 2012 plots
#### combine indicators into plots that look like the 2012 plots
#### combine indicators into plots that look like the 2012 plots
#### combine indicators into plots that look like the 2012 plots
#### combine indicators into plots that look like the 2012 plots
#### combine indicators into plots that look like the 2012 plots
#### combine indicators into plots that look like the 2012 plots
#### combine indicators into plots that look like the 2012 plots
#### combine indicators into plots that look like the 2012 plots
#### combine indicators into plots that look like the 2012 plots
#### combine indicators into plots that look like the 2012 plots
#### combine indicators into plots that look like the 2012 plots




################ added afterwards:
# 
# tmp = tmp[-which(abs(tmp$Trend) > 12),]
# 
# pp = ggplot(data = tmp,aes(x = PopUsCa,y = Trend))+
#   geom_point()+
#   scale_x_continuous(trans = "log",breaks = c(1000,10000,100000,1000000,10000000,100000000))+
#   geom_smooth(method = "lm")+
#   geom_hline(yintercept = 0)+
#   labs(title = "Long-term Trend decreases with population size", x = "Log(population size)", y = "Trend") 
# 
# 
# print(pp)

