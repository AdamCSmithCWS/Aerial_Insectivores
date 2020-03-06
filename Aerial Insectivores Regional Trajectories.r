#### jags function for SOCB model

pkgs = c("sp","maptools","rgdal","maptools","RColorBrewer","GISTools","ggplot2",
         "dplyr","tidyr","scales","geofacet","stringr","tools",
         "ggrepel","e1071","rjags")


# inst.pck = installed.packages()
# inst.pck2 = names(inst.pck[,2])
# 
# if(length(pkgs[-which(pkgs %in% inst.pck2)])>0){
#   install.packages(pkgs[-which(pkgs %in% inst.pck2)])
# }

for(i in 1:length(pkgs)){
  library(pkgs[i],character.only = T)
}



spslist_Rosenberg = read.csv(paste0("data/Rosenberg et al species list.csv"),
                   stringsAsFactors = F)
SOCBlist = read.csv("data/final survey and species lists SOCB 2018.csv",stringsAsFactors = F)

load("data/combined ai species list.RData")
allai = spai

allai = allai[-which(allai == "Cave Swallow")]




# Species groups n = 4 ----------------------------------------------------


Swallows_Swifts_Nightjars = c(allai[grep(pattern = " Sw",x = allai)],
                              allai[grep(pattern = " Night",x = allai)],
                              allai[grep(pattern = "oor",x = allai)],
                              allai[grep(pattern = "widow",x = allai)],
                              allai[grep(pattern = "Martin",x = allai)])
Flycatchers = c(allai[-which(allai %in% Swallows_Swifts_Nightjars)])
aiSOCB = SOCBlist[which(SOCBlist$national.Aerial.insectivores == "Included in group"),"English_Name"]
Raptors = c("Merlin","American Kestrel","Cooper's Hawk","Sharp-shinned Hawk","Peregrine Falcon")

splists = list(Swallows_Swifts_Nightjars = Swallows_Swifts_Nightjars,
               Flycatchers =  Flycatchers,
               aiSOCB = aiSOCB,
               Raptors = Raptors)



grpnms = c("Swallows_Swifts_Nightjars",
           "Flycatchers",
           "aiSOCB",
           "Raptors")
regtypes = c("continental","national","prov_state","stratum")


V = "GAM Smooth Only"
V = "with YE"
basedir = getwd()


source(paste0(basedir,"/socb.jags.func.r"))
source(paste0(basedir,"/transparency function.r"))



for(gr in grpnms){
#setwd(basedir)

 if(gr == "Raptors"){
   if(V == "GAM Smooth Only"){
     alldata_N = read.csv(paste0("data/raptor annual indices SmoothOnly BBS 2018.csv"),stringsAsFactors = F)
   v = ""
     }else{
     alldata_N = read.csv(paste0("data/raptor annual indices BBS 2018.csv"),stringsAsFactors = F)
     v = "YE"
     }
 }else{
   if(V == "GAM Smooth Only"){
     alldata_N = read.csv(paste0("data/AI annual indices SmoothOnly BBS 2018.csv"),stringsAsFactors = F)
     v = ""
     }else{
     alldata_N = read.csv(paste0("data/AI annual indices BBS 2018.csv"),stringsAsFactors = F)
     v = "YE"
     } 
 } 

  
  

sp.inc = splists[[gr]]


alldata_N = alldata_N[which(alldata_N$species %in% sp.inc &
                              alldata_N$Region_type %in% regtypes &
                              alldata_N$Trend_Time == "Long-term"),] #

outdir = paste0(basedir,"/",gr,v,"/")
dir.create(outdir)
#setwd(outdir)

dir.create(paste0(outdir,"output"))

for(regt in regtypes){
  
dir.create(paste0(outdir,"output/",regt))
  dirsave = paste0(outdir,"output/",regt,"/")
  alldata = alldata_N[which(alldata_N$Region_type == regt),] 
  
  
nybysp = table(alldata$species, exclude = c(NA))
nspeciestotal = length(nybysp)
fybysp = tapply(alldata$Year,alldata$species,min)



 alldata$lci = alldata[,"Index_q_0.025"]
 alldata$uci = alldata[,"Index_q_0.975"]
 
# alldata[which(alldata$Index == 0),"Index"] = 0.001 
# alldata[which(alldata$lci == 0),"lci"] = 0.0005 
# alldata[which(alldata$uci == 0),"uci"] = 0.0015 

today = Sys.Date()
alldataout = alldata[,c("species",
                        "espece",
                        "Year",
                        "Index",
                        "lci",
                        "uci",
                        "Region",
                        "Region_type",
                        "obs_mean")]

alldataout = alldataout[order(alldataout$Region_type,alldataout$Region,alldataout$Year),]

alldataout$species = factor(alldataout$species,ordered = T)

#regs = unique(alldataout[which(alldataout$Region_type %in% c("Continental","countrynum")),"Region"])
regs = unique(alldataout[,"Region"])


allResultExport = list()
length(allResultExport) = length(regs)


baseline = 1970
  q = 0
for(wind in regs){

    q = q+1

#    pg = unique(alldataout[which(alldataout$Region == wind),"Region_type"])
    
  names(allResultExport)[q] = wind
  
    spsel = unique(alldataout$species)

    dat = alldataout[which(alldataout$Region == wind),]
    
  spmiss = spsel[-which(spsel %in% dat$species)]
  
  nspy = table(dat$Year)
  psp = nspy/max(nspy)
  if(max(nspy) < 3 | mean(psp) < 0.5){next}
 
  fyr = max(c(min(as.integer(names(nspy)[which(nspy > 1)])),min(as.integer(names(psp)[which(psp > 0.5)])))) 
  # if(length(spmiss) > 0){
  #   print(paste("DATA missing for",spmiss,"in",ind))
  # }
  # 
  bsyr = baseline

   dat = dat[which(dat$Year >= bsyr),] 
   
  
  allResultExport[[wind]] = socb.jags(indices = dat,
            reg = regt,
            group = wind,
            ind = "Index",
            lci = "lci",
            uci = "uci",
            Year = "Year",
            base.yr = bsyr)
  
  print(paste(wind))
  
  
  
}
  
 
  save(allResultExport,
       file = paste0(dirsave,gr," ",regt," SCB results GAM ",baseline,".RData"))
  
 
} 
  
}#grpnms end modeling loop







# Plotting of geofacet plots ----------------------------------------------


for(gr in grpnms){
  #setwd(basedir)
  
  if(gr == "Raptors"){
    if(V == "GAM Smooth Only"){
      alldata_N = read.csv(paste0("data/raptor annual indices SmoothOnly BBS 2018.csv"),stringsAsFactors = F)
      v = ""
    }else{
      alldata_N = read.csv(paste0("data/raptor annual indices BBS 2018.csv"),stringsAsFactors = F)
      v = "YE"
    }
  }else{
    if(V == "GAM Smooth Only"){
      alldata_N = read.csv(paste0("data/AI annual indices SmoothOnly BBS 2018.csv"),stringsAsFactors = F)
      v = ""
    }else{
      alldata_N = read.csv(paste0("data/AI annual indices BBS 2018.csv"),stringsAsFactors = F)
      v = "YE"
    } 
  } 
  
  

  sp.inc = splists[[gr]]
  
  regtypes = c("continental","national","prov_state","stratum")
  
  alldata_N = alldata_N[which(alldata_N$species %in% sp.inc &
                                alldata_N$Region_type %in% regtypes &
                                alldata_N$Trend_Time == "Long-term"),] #
  
  outdir = paste0(basedir,"/",gr,v,"/")
  #setwd(outdir)
  
 
  
  baseline = 1970
  
  for(regt in regtypes[3:4]){
    dirsave = paste0(outdir,"output/",regt,"/")
    
    
    alldata = alldata_N[which(alldata_N$Region_type == regt),] 
    
    
    nybysp = table(alldata$species, exclude = c(NA))
    nspeciestotal = length(nybysp)
    fybysp = tapply(alldata$Year,alldata$species,min)
    
    
    
    alldata$lci = alldata[,"Index_q_0.025"]
    alldata$uci = alldata[,"Index_q_0.975"]
    
    # alldata[which(alldata$Index == 0),"Index"] = 0.001 
    # alldata[which(alldata$lci == 0),"lci"] = 0.0005 
    # alldata[which(alldata$uci == 0),"uci"] = 0.0015 
    
    alldataout = alldata[,c("species",
                            "espece",
                            "Year",
                            "Index",
                            "lci",
                            "uci",
                            "Region",
                            "Region_type",
                            "obs_mean")]
    
    alldataout = alldataout[order(alldataout$Region_type,alldataout$Region,alldataout$Year),]
    
    alldataout$species = factor(alldataout$species,ordered = T)
    
    #regs = unique(alldataout[which(alldataout$Region_type %in% c("Continental","countrynum")),"Region"])
    regs = unique(alldataout[,"Region"])
    
  
  
  
  load(paste0(dirsave,gr," ",regt," SCB results GAM ",baseline,".RData"))
  
  

# geofacet grid arrangement -----------------------------------------------

  
  stprovfacet = read.csv(paste0(basedir,"/data/BBB_StateProvCWS_facet_grid2.csv"),stringsAsFactors = F)
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
  
  

# Change colour scale ----------------------------------------------------

  
  cols = function(x){
    cl = "stable"
    if(x > 1){cl = "inc"}
    if(x < 0.5 & x > 0.25){cl = "dec"}
    if(x < 0.25){cl = "Ldec"}
    
    return(cl)
  }
  
  
  if(regt == "stratum"){
    ############# facet plot for strata
  ############# facet plot for strata
  ############# facet plot for strata
  ############# facet plot for strata
  ############# facet plot for strata
  ############# facet plot for strata
  rc = 0
for(rr in unique(alldataout[which(alldataout$Region_type == "stratum"),"Region"])){

  tmp = allResultExport[[rr]]$alldatout
  if(length(tmp) < 1){next}
  rc = rc+1
  tmp$Region = rr
  tmp$state = strsplit(rr,fixed = T,split = "-")[[1]][2]
  tmp$BCR = strsplit(rr,fixed = T,split = "-")[[1]][3]
  names(tmp)[c(2,4,6,7,8)] = c("lci","Index","uci","nspecies","splist")
  tmp = tmp[,c("Year","lci","Index","uci","nspecies","Region","state","BCR","splist")]
  if(rc == 1){
    tpout = tmp
  }else{
    tpout = rbind(tpout,tmp)
  }
  #### compile all Regions predictions into a single data frame with state names as well as strata names
}
  
  
 

    # indt$prst = gsub(str_extract(indt$stratcode,"-(.*)-"),pattern = "-",replacement = "",fixed = T)
    # indt$prst = as.character(factor(indt$prst,levels = unique(stprovfacet$code)))
  indt = tpout
    indt$prst = as.character(factor(indt$state,levels = unique(stprovfacet$code)))
    labsbcr = indt[which(indt$Year == max(indt$Year,na.rm = T)),]
       
    for(pp in labsbcr$Region){
      indt[which(indt$Region == pp),"col"] = cols(labsbcr[which(labsbcr$Region == pp),"Index"])
    }
    
    indt$col = factor(indt$col,levels = c("Ldec","dec","stable","inc"),ordered = T)
    
    
    myColors <- brewer.pal(10,"RdYlBu")[c(1,3,4,10)]
    names(myColors) <- (levels(indt$col))
    
    
    colScale <- scale_colour_manual(values = myColors, aesthetics = c("colour","fill"))
    ylm = range(indt$Index)
    
    
     suppressMessages(
       ptraj <- ggplot(data = indt) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(colour = grey(0.2)),
            strip.background = element_rect(fill = grey(0.97)),#strcol #, colour = grey(0.9), size = NULL, lineRegion_type = NULL, color = NULL, inherit.blank = FALSE
            #axis.line = element_line(colour = "black"),
            legend.position = "none") +
      labs(title = paste(gr,"indicator by strata (separate BCRS) within Provinces and States"), x = "", y = "mean change since 1970") +

      #geom_pointrange(data = indt, aes(x = Year, y = mean,ymin = lci,ymax = uci),alpha = 0.3) +
       geom_line(data = indt, aes(x = Year, y = Index,group = Region,colour = col)) +#, colour = stratcode)) +
       geom_ribbon(data = indt, aes(x = Year, ymin = lci, ymax = uci,group = Region,fill = col),alpha = 0.12)+ #, fill = stratcode , fill = grey(0.8)
       geom_text_repel(data = labsbcr, aes(x = Year, y = Index,label = BCR),
                       size = 2)+
       colScale+
       #geom_hline(yintercept = 0,colour = gray(0.8))+
       facet_geo(facets = ~ prst,grid = stprovfacet,scales = "free_y", label = "code")+
       coord_cartesian(ylim = ylm)+
       scale_y_continuous(trans = "log", breaks = ((c(-98,-90,-75,-50,0,200,400)/100)+1))+ #,labels = c("-75","-50","0"))+
       scale_x_continuous(breaks = c(1970,2018),minor_breaks = c(1980,1990,2000,2010)) #limits = c(1966, 2018), oob = squish, 
 )
    pdf(file = paste0(outdir,gr," indicator geofacet plot by ",regt,".pdf"),
         height = 8.5,width = 11)
    print(ptraj)
    dev.off()
      # for(p in stprovfacet$code)){
    #
    # }#p

  }else{
  

    
    
    
    
    ############# by state-prov
    rc = 0
    for(rr in unique(alldataout[which(alldataout$Region_type == regt),"Region"])){
      
      tmp = allResultExport[[rr]]$alldatout
      if(length(tmp) < 1){next}
      rc = rc+1
      tmp$Region = rr
      tmp$state = rr
      names(tmp)[c(2,4,6,7,8)] = c("lci","Index","uci","nspecies","splist")
      tmp = tmp[,c("Year","lci","Index","uci","nspecies","Region","state","splist")]
      if(rc == 1){
        tpout = tmp
      }else{
        tpout = rbind(tpout,tmp)
      }
      #### compile all Regions predictions into a single data frame with state names as well as strata names
    }
    
    
    
    indt = tpout
    indt$prst = as.character(factor(indt$state,levels = unique(stprovfacet$code)))
    labsbcr = indt[which(indt$Year == max(indt$Year,na.rm = T)),]
    
    
    for(pp in labsbcr$Region){
      indt[which(indt$Region == pp),"col"] = cols(labsbcr[which(labsbcr$Region == pp),"Index"])
    }
    
    
    indt$col = factor(indt$col,levels = c("Ldec","dec","stable","inc"),ordered = T)
    
    
    myColors <- brewer.pal(10,"RdYlBu")[c(1,3,4,10)]
    names(myColors) <- (levels(indt$col))
    
    
    colScale <- scale_colour_manual(values = myColors, aesthetics = c("colour","fill"))
    
    
    ylm = range(indt$Index)
    
    suppressMessages(
      ptraj <- ggplot(data = indt) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(colour = grey(0.2)),
            strip.background = element_rect(fill = grey(0.97)),#strcol #, colour = grey(0.9), size = NULL, linetype = NULL, color = NULL, inherit.blank = FALSE
            #axis.line = element_line(colour = "black"),
            legend.position = "none") +
      labs(title = paste(gr,"indicator by strata (separate BCRS) within Provinces and States"), x = "", y = "mean change since 1970") +
      
      #geom_pointrange(data = indt, aes(x = Year, y = mean,ymin = lci,ymax = uci),alpha = 0.3) +
      geom_line(data = indt, aes(x = Year, y = Index,group = state,colour = col)) +
      geom_ribbon(data = indt, aes(x = Year, ymin = lci, ymax = uci,group = Region,fill = col),alpha = 0.12)+ #, fill = stratcode , fill = grey(0.8)
      colScale + #, colour = stratcode)) +
      
      
      #geom_pointrange(data = indt, aes(x = Year, y = mean,ymin = lci,ymax = uci),alpha = 0.3) +
      geom_line(data = indt, aes(x = Year, y = Index,group = Region,colour = col)) +#, colour = stratcode)) +
     #geom_text_repel(data = labsbcr, aes(x = Year, y = Index,label = BCR),
      #                size = 2)+
      #geom_hline(yintercept = 0,colour = gray(0.8))+
      coord_cartesian(ylim = ylm)+
      scale_y_continuous(trans = "log", breaks = ((c(-98,-90,-75,-50,0,200,400)/100)+1))+ #,labels = c("-75","-50","0"))+
      facet_geo(facets = ~ prst,grid = stprovfacet,scales = "free_y", label = "code")+
      scale_x_continuous(breaks = c(1970,2018),minor_breaks = c(1980,1990,2000,2010)) #limits = c(1966, 2018), oob = squish, 
)
    pdf(file = paste0(outdir,gr," indicator geofacet plot by ",regt,".pdf"),
        height = 8.5,width = 11)
    print(ptraj)
    dev.off()


  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ############# commmunity composition plots
   # if(gr == "all"){
    
    
    
    
    
    
    
    
    
    ############# species-level trajectory plots.
    
    
    rc = 0
    for(rr in unique(alldataout[which(alldataout$Region_type == regt),"Region"])){
      
      tmp1 = allResultExport[[rr]]$alldatout
      if(length(tmp1) < 1){next}
      clx = grep(names(tmp1),pattern = ".se",fixed = T)
      spsl = gsub(names(tmp1)[clx],pattern = ".se",replacement = "",fixed = T)
      for(ss in spsl){
        tmp2 = tmp1[,c("Year",ss,paste0(ss,".se"))]
        names(tmp2) = c("Year","Index","se")
        tmp2$species = ss
        if(ss == spsl[1]){
          
          tmp = tmp2
        }else{
          tmp = rbind(tmp,tmp2)
        }
      }
      
      tmp$lci = tmp$Index-(tmp$se)*1.96
      tmp$uci = tmp$Index+(tmp$se)*1.96
      
      tmp$Index = exp(tmp$Index)
      tmp$lci = exp(tmp$lci)
      tmp$uci = exp(tmp$uci)
      
      rc = rc+1
      tmp$Region = rr
      tmp$state = rr
       if(rc == 1){
        tpout = tmp
      }else{
        tpout = rbind(tpout,tmp)
      }
      #### compile all Regions predictions into a single data frame with state names as well as strata names
    }
    
    
    indt1 = tpout
    indt1$prst = as.character(factor(indt1$state,levels = unique(stprovfacet$code)))
    labsbcr = indt1[which(indt1$Year == max(indt1$Year,na.rm = T)),]
    indt1 = indt1[order(indt1$prst,indt1$Year,indt1$species),]
 

    pdf(file = paste0(outdir,gr,"species lines geofacet plot by stateprov.pdf"),
        height = 8.5,width = 11)

    for(ss in unique(indt1$species)){
      indt = indt1[which(indt1$species == ss),]
    
      
      labsbcr = indt[which(indt$Year == max(indt$Year,na.rm = T)),]
      #labsbcr = merge(labsbcr,spslist[,c("English_Name","Species_ID")],by.x = "species",by.y = "English_Name")
      

      for(pp in labsbcr$Region){
        indt[which(indt$Region == pp),"col"] = cols(labsbcr[which(labsbcr$Region == pp),"Index"])
      }
      
      
      indt$col = factor(indt$col,levels = c("Ldec","dec","stable","inc"),ordered = T)
      
      
      # myColors <- brewer.pal(10,"RdYlBu")[c(1,3,4,10)]
      # names(myColors) <- (levels(indt$col))
      # 
      
      colScale <- scale_colour_manual(values = myColors, aesthetics = c("colour","fill"))
      
      ylm = range(indt$Index)
      
      
      suppressMessages(
    ptraj <- ggplot(data = indt) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(colour = grey(0.2)),
            strip.background = element_rect(fill = grey(0.97)),#strcol #, colour = grey(0.9), size = NULL, linetype = NULL, color = NULL, inherit.blank = FALSE
            #axis.line = element_line(colour = "black"),
            legend.position = "none") +
      labs(title = paste(ss,"Population trajectories by Provinces and States"), x = "", y = "mean change since 1970") +
      
      #geom_pointrange(data = indt, aes(x = Year, y = mean,ymin = lci,ymax = uci),alpha = 0.3) +
      geom_line(data = indt, aes(x = Year, y = Index,group = state,colour = col)) +
      #geom_point(data = indt, aes(x = Year, y = Index,group = state,colour = species)) +#, colour = stratcode)) +
      
      
      geom_ribbon(data = indt, aes(x = Year, ymin = lci, ymax = uci,group = state,fill = col),alpha = 0.12)+ #, fill = stratcode , fill = grey(0.8)
      colScale + #, colour = stratcode)) +
      # geom_text_repel(data = labsbcr, aes(x = Year, y = Index,label = Species_ID,colour = Species_ID),
      #                 size = 1.5)+
      #geom_hline(yintercept = 0,colour = gray(0.8))+
      coord_cartesian(ylim = ylm)+
      scale_y_continuous(trans = "log", breaks = ((c(-98,-90,-75,-50,0,200,400)/100)+1))+ #,labels = c("-75","-50","0"))+
      facet_geo(facets = ~ prst,grid = stprovfacet,scales = "free_y", label = "code")+
      scale_x_continuous(breaks = c(1970,2018),minor_breaks = c(1980,1990,2000,2010)) #limits = c(1966, 2018), oob = squish, 
)
      # pdf(file = paste0(ss,"Grassland species lines geofacet plot by stateprov.pdf"),
    #     height = 8.5,width = 11)
     print(ptraj)
    # dev.off()
    # 
    }
  dev.off()
  
  #  }##3 end if gr == all
    
    
  }##end if regt == state_prov
  
  
  }#end regt loop 
    
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

