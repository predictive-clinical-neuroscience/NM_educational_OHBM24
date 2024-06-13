# If using these tools or functions please cite Tervo-Clemmens et al., 2023, Nature Communications https://www.nature.com/articles/s41467-023-42540-8


### load packages
if(require(pacman)==FALSE){
  install.packages("pacman")
}
pacman::p_load(scales, mgcv, lme4, gamm4, dplyr, ggplot2, gratia, reshape, cowplot, tidyr, install = TRUE)

####convenience functions######
ttor<-function(ts,nobs){
  rfromt<-sqrt((t^2)/(t^2+nobs))
}
find_covars_gam <- function(fml, ...) {
  ind <- as.character(fml)[3]
  # s(x1) + x2 +s(x3,"re") -> x1, x2
  vars <- unlist(strsplit(ind, "\\+"))  # formula split by +
  vars <- gsub(" ", "", vars) # remove spaces
  vars <- gsub("\\w+\\((.*?)\\)", "\\1", vars) # remove surrounding s()
  # remove random effect
  no_re <- grep("re[\"']", vars, value=T, invert=T)
  no_re <- gsub(",.*", "", no_re) # nothing after comma
  # remove anything else (likely passed in agevar)
  if (length(list(...)) != 0L) {
    no_re <- no_re[ ! no_re  %in% c(...)]
  }
  return(no_re)
}
too_small <- function(x) abs(x) < 10^-15
clip_on_sig <- function(ci){
  # if confidence interval includes zero
  # signs of x and y will be different, -x * +y  < 0
  # or if both high and low are extremly close to zero
  not_sig <- ci$ci_low * ci$ci_high < 0 |
    (too_small(ci$ci_low) & too_small(ci$ci_high))
  ci$mean_dff_clip <- ci$mean_dff
  ci$mean_dff_clip[not_sig] <- 0
  return(ci)
}
clip_on_siggratia <- function(ci){
  # if confidence interval includes zero
  # signs of x and y will be different, -x * +y  < 0
  # or if both high and low are extremly close to zero
  not_sig <- ci$lower * ci$upper < 0 |
    (too_small(ci$lower) & too_small(ci$upper))
  ci$sig <- 1
  ci$sig[not_sig] <- 0
  return(ci)
}
##############Generate Fits#######################
##################################################
###convenience function to predict from mgcv models
mgcvgampreddata<-function(model,predvar,idvar=NULL,interval_inc=.1,varycovarsname=NULL,varycovarlevels=NULL){
  if (class(model)[1]=="gamm"){
    model<-model$gam
  }
  modeldata<-data.frame(ydata=model$y, predvar=model$model[, predvar])
  preddata<-data.frame(var=seq(min(modeldata$predvar),max(modeldata$predvar),by=interval_inc))
  names(preddata)[1]<-predvar
  if (!identical(find_covars_gam(model$formula, predvar),character(0))){
    for (cv in find_covars_gam(model$formula,predvar)){
      x <- model$model[, cv]
      if (is.character(x) || is.factor(x)){
        warning("gam w/character or factor covar, setting all sim to the first obs for character/first level if factor")
        y <- x[1]
        if (class(x)=="factor"){y<-levels(x)[1]}
        print(sprintf("covar % set to level %s",cv,y))
      } else {
        y <- mean(x, na.rm=T)
      }
      preddata[, cv] <- y
    }
    #}
  }else if(is.null(varycovarsname)){
    preddata$cov<-"no cov"
  }
  names(preddata)<-c(predvar,find_covars_gam(model$formula, predvar))
  if (identical(find_covars_gam(model$formula, predvar),character(0)) && is.null(varycovarsname)){
    names(preddata)<-c(predvar,"nullcovar") 
  }
  
  if (!is.null(varycovarsname)){
    require(reshape)
    orignameorder<-names(preddata)
    preddata[,varycovarsname]<-NULL
    preddata<-reshape::expand.grid.df(preddata,data.frame(varycovar=varycovarlevels))
    names(preddata)[names(preddata)=="varycovar"]<-varycovarsname
  }
  
  yhats <- predict(model,preddata,se.fit=TRUE)
  preddata<-cbind(preddata,yhats$fit,yhats$se.fit)
  names(preddata)<-c(predvar,find_covars_gam(model$formula, predvar),"fit","se")
  if (identical(find_covars_gam(model$formula, predvar),character(0)) && is.null(varycovarsname)){
    names(preddata)<-c(predvar,"nullcovar","fit","se")
  }else if(identical(find_covars_gam(model$formula, predvar),character(0)) && !is.null(varycovarsname)){
    names(preddata)<-c(predvar,varycovarsname,"fit","se")
  }
  preddata$CI<-2*preddata$se
  return(preddata)
}

###convenience function to iteratively run mgcv gam/gamm and save outputs
mgcvscalefits<-function(df,outcomevars,predvars,idvar=NULL,mformula,interval_inc=.1,scale=TRUE){
  pairs<-as.data.frame(expand.grid(outcomevars,predvars))
  names(pairs)<-c("outcome","pred")
  allscalefits<-lapply(1:nrow(pairs),function(p){
    if (scale){
      df$outcome<-scale(unlist(df[,as.character(pairs$outcome[p])])) 
    }else{
      df$outcome<-unlist(df[,as.character(pairs$outcome[p])]) 
    }
    df$pred<-unlist(df[,as.character(pairs$pred[p])])
    #print(pairs[p,])
    model<-mgcv::gam(mformula,data=df)

    modelpred<-mgcvgampreddata(model,predvar="pred",idvar=idvar,interval_inc = interval_inc)
    modelpred$fitscale<-scales::rescale(modelpred$fit,c(0,1))
    modelpred$fitz<-base::scale(modelpred$fit)[,1]
    modelpred$outcome<-as.character(pairs$outcome[p])
    
    modelpred$samplesize<-nrow(df)
    nps<-length(summary(model)$s.pv)
      if (nps==1){
        modelpred$modelp<-summary(model)$s.pv
        modelpred$anovamodelp<-anova.gam(model)$s.table[1,4]
        modelpred$Fstat<-summary(model)$s.table[1,3]
        modelpred$edf<-summary(model)$edf
      } else {
        modelpred$modelp<-summary(model)$s.pv[which(grepl("pred",row.names(summary(model)$s.table)))]
        modelpred$anovamodelp<-anova.gam(model)$s.table[which(grepl("pred",row.names(anova.gam(model)$s.table))),4]
        modelpred$Fstat<-anova.gam(model)$s.table[which(grepl("pred",row.names(anova.gam(model)$s.table))),3]
        modelpred$edf<-summary(model)$edf[which(grepl("pred",row.names(summary(model)$s.table)))]
      }
    return(modelpred)
  })
  allscalefits<-dplyr::bind_rows(allscalefits)
  return(allscalefits) 
}

##############Plotting Functions##################
##################################################
ggfinalize<-function (p, ajust = 0.05) {
  p$labels$y <- paste0(p$labels$y, "\n")
  p$labels$x <- paste0("\n", p$labels$x)
  p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                         text = element_text(size = 20), axis.text.y = element_text(hjust = ajust, 
                                                                                    color = "black"), axis.text.x = element_text(vjust = ajust, 
                                                                                                                                 color = "black"))
}
plotfitsfunc<-function(fitsall,outcomegroupdata,datasetdata,points=TRUE,rawdat=FALSE){
  
fitsall<-merge(fitsall,outcomegroupdata,by=c("outcome"))
fitsall<-fitsall %>% group_by(outcome) %>% mutate(adjusted.p=p.adjust(unique(modelp)))
fitsallsig<-fitsall
fitsallsig$sig<-as.character(dplyr::if_else(fitsallsig$adjusted.p < .05,1,0))

###scaled points#######
accvars<-outcomegroupdata$outcome[outcomegroupdata$type=="acc" & !is.na(outcomegroupdata$outcome)]
latvars<-outcomegroupdata$outcome[outcomegroupdata$type=="lat" & !is.na(outcomegroupdata$outcome)]

pointsacc<-mgcvscalefits(datasetdata,outcomevars = c(accvars),predvars = "age",mformula =  as.formula('outcome~s(pred)'),interval_inc = 3)
pointsacc$type<-"acc"

pointslat<-mgcvscalefits(datasetdata,outcomevars = c(latvars),predvars = "age",mformula =  as.formula('outcome~s(pred)'),interval_inc = 3)
pointslat$type<-"lat"

pointssall<-rbind(pointsacc,pointslat)
pointssallsig<-pointssall[pointssall$outcome %in% unique(pointssall$outcome),]
pointssallsig<-merge(pointssall,outcomegroupdata,by=c("outcome","type"))

pointssallsig<-pointssallsig[pointssallsig$outcome %in% unique(fitsall$outcome),]

if (points==TRUE){
ggfitall<-ggplot()+
  geom_line(data=fitsallsig[fitsallsig$group=="Composite",],aes(x=pred,y=fit),colour="black",size=1.5,alpha=.5)+
  geom_line(data=fitsallsig[fitsallsig$group!="Composite",],aes(x=pred,y=fit,colour=group),size=.5,alpha=1)+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+
  geom_point(data=pointssallsig[pointssallsig$group!="Composite",],aes(x=pred,y=fit,shape=group),size=2)+
  facet_grid(rows=vars(type),scales="free")+theme(strip.text.y= element_text(angle = 180))+
  scale_shape_manual(values=c(0,1,2,3,4,5,6,7))+
  scale_colour_manual(values=rep("black",length(unique(fitsallsig$group))))
ggfitall<-ggfinalize(ggfitall)+xlab("Age (years)")+ylab("Executive Function\n (z-score)")+theme(legend.position = "top",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))
}else{
  ggfitall<-ggplot()+
    geom_line(data=fitsallsig[fitsallsig$group=="Composite",],aes(x=pred,y=fit),colour="black",size=1.5,alpha=.5)+
    geom_line(data=fitsallsig[fitsallsig$group!="Composite",],aes(x=pred,y=fit,colour=group),size=.5,alpha=1)+
    scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+
    #geom_point(data=pointssallsig[pointssallsig$group!="Composite",],aes(x=pred,y=fit,shape=group),size=2)+
    facet_grid(rows=vars(type),scales="free")+theme(strip.text.y= element_text(angle = 180))+
    #scale_shape_manual(values=c(0,1,2,3,4,5,6,7))+
    scale_colour_manual(values=rep("black",length(unique(fitsallsig$group))))
  ggfitall<-ggfinalize(ggfitall)+xlab("Age (years)")+ylab("Executive Function\n (z-score)")+theme(legend.position = "top",legend.title = element_blank())+
    theme(strip.text.y= element_text(angle = 360))
}

if(rawdat==TRUE){
  datasetdatalong<-datasetdata %>% tidyr::pivot_longer ( cols = -age,names_to = "outcome",values_to = "value")
  datasetdatalong<-datasetdatalong[datasetdatalong$outcome %in% unique(fitsallsig$outcome),]
  datasetdatalong<-merge(datasetdatalong,outcomegroupdata,by="outcome")
  ggfitall<-ggplot()+
    geom_line(data=fitsallsig[fitsallsig$group=="Composite",],aes(x=pred,y=fit),colour="black",size=2.0,alpha=.5)+
    geom_line(data=fitsallsig[fitsallsig$group!="Composite",],aes(x=pred,y=fit,colour=group),size=1.5,alpha=1)+
    scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+
    geom_point(data=datasetdatalong,aes(x=age,y=value,shape=group),size=2,alpha=.33)+
    facet_grid(rows=vars(type),scales="free")+theme(strip.text.y= element_text(angle = 180))+
    scale_colour_manual(values=rep("black",length(unique(fitsallsig$group))))
  ggfitall<-ggfinalize(ggfitall)+xlab("Age (years)")+ylab("Executive Function\n (z-score)")+theme(legend.position = "top",legend.title = element_blank())+
    theme(strip.text.y= element_text(angle = 360))
}

return(ggfitall)
}



#############Derivative Functions#################
##################################################
###convenience plot funcs#######
too_small <- function(x) abs(x) < 10^-15
clip_on_sig <- function(ci){
  # if confidence interval includes zero
  # signs of x and y will be different, -x * +y  < 0
  # or if both high and low are extremly close to zero
  not_sig <- ci$ci_low * ci$ci_high < 0 |
    (too_small(ci$ci_low) & too_small(ci$ci_high))
  ci$mean_dff_clip <- ci$mean_dff
  ci$mean_dff_clip[not_sig] <- 0
  return(ci)
}
lunaize_geomraster<-function(x){
  x+
    theme_bw()+
    theme(
      axis.line = element_line(colour = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.title.y     = element_blank(),
      axis.ticks.y     = element_blank(),
      axis.text.y      = element_blank(),
      axis.title.x     = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.text.x      = element_blank(),
      legend.position  = "none")+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
}
lunaize_geomrasterxkeep<-function(x){
  x+
    theme_bw()+
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.title.y     = element_blank(),
      axis.ticks.y     = element_blank(),
      axis.text.y      = element_blank(),
      legend.position  = "none")+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
}
##################################
mgcvgam_growthrate<-function(m, agevar, idvar = NULL, n.iterations = 10000, qnt = c(0.025,0.975),interval_inc=.1,scaletomaximum=FALSE) 
{
    require(gratia)
    #deriv_gratia<-fderiv(m)
    v <- m$model[, agevar]
    agevarintervaldata <- mgcvgampreddata(m,agevar,idvar,interval_inc=.1)
    CIlevel<-qnt[2]-qnt[1]
    f1deriv<-gratia::derivatives(m,term=agevar,partial_match=TRUE,interval="simultaneous",level=CIlevel,type="backward",newdata=agevarintervaldata)
    names(f1deriv)<-c("smooth","term","ages","mean_dff","se","crit","ci_low","ci_high")
    f1deriv$ages<-agevarintervaldata[,agevar]
    ciclip<-clip_on_sig(f1deriv)
  return(ciclip)
}
mgcvgam_growthrate_multiplot<-function(df,outcomevars,predvars,plotxlim=c(8,36),breaks=(c(10,15,20,25,30,35)),
                                                        idvar=NULL,mformula,zscale=TRUE,totaltiles=NULL,xscaleplot=TRUE,derivrangemanual=c(-.2,.2),
                                                        derivcolourmanual=c("#2f42bde6","#c9270e"),
                                                        xrangeticks=TRUE,xrangetickvals=NULL,datarangegrey=FALSE,devageguides=TRUE
){ ###derivative raster 
  if(is.null(xrangetickvals)){xrangetickvals<-plotxlim}
  pairs<-as.data.frame(expand.grid(outcomevars,predvars))
  names(pairs)<-c("outcome","pred")
  ggrout<-NULL
  for (p in 1:nrow(pairs)){
    if (zscale){
      df$outcome<-base::scale(unlist(df[,as.character(pairs$outcome[p])]))
    }else{
      df$outcome<-unlist(df[,as.character(pairs$outcome[p])])
    }
    df$pred<-unlist(df[,as.character(pairs$pred[p])])
    model<-mgcv::gam(mformula,data=df)
    gammodel<-model
    ggr<-mgcvgam_growthrate(model,agevar="pred",idvar=idvar)
    ggr$var<-as.character(pairs$outcome[p])
    ggrout<-rbind(ggrout,ggr)
  }
  ggrouttoplot<-ggrout[!is.na(ggrout$mean_dff),]
  ggrouttoplotclip<-clip_on_sig(ggrouttoplot)
  ggrouttoplotclip$se<-(ggrouttoplotclip$ci_high-ggrouttoplotclip$mean_dff)/2
  ggrouttoplotclip$sd<-ggrouttoplotclip$se*sqrt(10000) #n interations from model ###change to sample size
  ggrouttoplotclip$tstat<-(ggrouttoplotclip$mean_dff-0)/(ggrouttoplotclip$sd/sqrt(10000))
  
  ggrouttoplotclip$tstat[ggrouttoplotclip$mean_dff_clip==0]<-0
  deriv_range<-range(ggrouttoplotclip$tstat)
  fillvar="tstat"
  if (zscale){
    deriv_range<-range(ggrouttoplotclip$mean_dff_clip)
    fillvar="mean_dff_clip"
  }
  
  if (sign(deriv_range[1])==sign(deriv_range[2])){
    if (sign(max(deriv_range))==1){
      deriv_range[deriv_range==min(deriv_range)]<-0
    } else {
      deriv_range[deriv_range==max(deriv_range)]<-0
    }
  }
  
  if(!is.null(derivrangemanual)){
    deriv_range<-derivrangemanual
    #print(sprintf("using manual deriv range %s %s",deriv_range[1],deriv_range[2]))
    #print("values outside of range will be set to max via oob=squish")
  }
  
  allplots<-list()
  for (vi in 1:length(unique(ggrouttoplotclip$var))){
    #print(vi)
    v<-unique(ggrouttoplotclip$var)[vi]
    ci<-ggrouttoplotclip[ggrouttoplotclip$var==v,]
    tile <- ggplot(ci) + aes_string(x = "ages", y = 1, fill = fillvar) + geom_raster(interpolate = TRUE) +
      scale_fill_gradient2(low = derivcolourmanual[1], mid = "white", high = derivcolourmanual[2], midpoint = 0, space = "Lab", breaks = sort(c(0, deriv_range)), limits = deriv_range,oob=squish)+xlim(plotxlim)
    #+ggtitle(v)
    if (vi !=length(unique(ggrouttoplotclip$var))){
      tile_luna <- lunaize_geomraster(tile) + theme(text = element_text(size = 20))+theme(axis.title.x=element_blank())#,axis.text.x=element_text(color="white"),axis.ticks.x=element_line(color="white"))
    }else{
      #print("final")
      tile_luna <- lunaize_geomraster(tile) + theme(text = element_text(size = 20))+theme(axis.title.x=element_blank())
    }
    if(datarangegrey){
      require(ggpattern)
      tile_luna<-tile_luna+geom_rect(xmin = xrangetickvals[1], xmax = min(df$pred,na.rm=TRUE),   ymin =0.5, ymax = 1.5,   fill = "grey88") 
      tile_luna<-tile_luna+geom_rect(xmin = max(df$pred,na.rm=TRUE), xmax = xrangetickvals[2],   ymin =0.5, ymax = 1.5,   fill = "grey88") 
    }
    if(xrangeticks){
      tile_luna<-tile_luna+geom_segment(linetype = 1, size=.5,colour = "black", aes(x = min(df$pred,na.rm=TRUE), xend =  min(df$pred,na.rm=TRUE), y = 0.5, yend = 1.5))
      tile_luna<-tile_luna+geom_segment(linetype = 1, size=.5,colour = "black", aes(x = max(df$pred,na.rm=TRUE), xend = max(df$pred,na.rm=TRUE), y = 0.5, yend = 1.5))
    }
    
    allplots[[vi]]<-tile_luna
    legend <- cowplot::get_legend(tile)
    tilegrid<-cowplot::plot_grid(plotlist=allplots,ncol = 1)
    returnplot<-tilegrid
  }
  fillplots<-NA
  pairswithfill<-pairs
  if (!is.null(totaltiles)){
    fillplots<-list()
    print(totaltiles)
    addtiles<-totaltiles-length(allplots)
    cifill<-ci
    cifill$tstat<-0
    cifill$mean_dff_clip<-0
    tilefill <- ggplot(cifill) + aes_string(x = "ages", y = 1, fill = fillvar) + geom_raster(interpolate = TRUE) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, space = "Lab", breaks = sort(c(0, deriv_range)), limits = deriv_range)+xlim(plotxlim)
    tilefillluna <- lunaize_geomraster(tilefill) + theme(text = element_text(size = 20))+theme(axis.title.x=element_blank())
    for (ai in seq(1,addtiles)){
      print(ai)
      #ati<-length(allplots)+ai
      fillplots[[ai]]<-tilefillluna
    }
    tilegrid<-cowplot::plot_grid(plotlist=c(allplots,fillplots),ncol = 1)
    returnplot<-tilegrid
    fillforpais<-data.frame(matrix(ncol = ncol(pairs), nrow = addtiles))
    names(fillforpais)<-names(pairs)
    pairswithfill<-rbind(pairs,fillforpais)
  }
  tilescaleplot<-NA
  if(xscaleplot){
    ciscaleplot<-ci
    ciscaleplot$tstat<-0
    ciscaleplot$mean_dff_clip<-0
    tilescaleplot <- ggplot(ciscaleplot) + aes_string(x = "ages", y = 1, fill = fillvar) + geom_raster(interpolate = TRUE) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, space = "Lab", breaks = sort(c(0, deriv_range)), limits = deriv_range)+
      scale_x_continuous(limits = plotxlim,breaks=breaks)
    tilescaleplot <- lunaize_geomrasterxkeep(tilescaleplot)+theme(text = element_text(size = 20))+theme(axis.title.x=element_blank())
    if(xrangeticks){
      tilescaleplot<-tilescaleplot+geom_segment(linetype = 1, colour = "black", aes(x = min(df$pred,na.rm=TRUE), xend =  min(df$pred,na.rm=TRUE), y = 0.5, yend = 1.5))
      tilescaleplot<-tilescaleplot+geom_segment(linetype = 1, colour = "black", aes(x = max(df$pred,na.rm=TRUE), xend = max(df$pred,na.rm=TRUE), y = 0.5, yend = 1.5))
    }
    if (devageguides){
      tilescaleplot<-tilescaleplot+geom_segment(linetype = 2, colour = "black", aes(x = 12, xend =  12, y = 0.5, yend = 1.5))
      tilescaleplot<-tilescaleplot+geom_segment(linetype = 2, colour = "black", aes(x = 15, xend=15, y = 0.5, yend = 1.5))
      tilescaleplot<-tilescaleplot+geom_segment(linetype = 2, colour = "black", aes(x = 18, xend=18, y = 0.5, yend = 1.5))
    }
    if(datarangegrey){
      tilescaleplot<-tilescaleplot+geom_rect(xmin = xrangetickvals[1], xmax = min(df$pred,na.rm=TRUE),   ymin =0.5, ymax = 1.5,   fill = "grey88") 
      tilescaleplot<-tilescaleplot+geom_rect(xmin = max(df$pred,na.rm=TRUE), xmax = xrangetickvals[2],   ymin =0.5, ymax = 1.5,   fill = "grey88") 
    }
    tilescaleplotlist<-list(tilescaleplot)
    if(!is.null(totaltiles)){
      returnplot<-cowplot::plot_grid(plotlist=c(allplots,fillplots,tilescaleplotlist),ncol = 1)
    }else{
      returnplot<-cowplot::plot_grid(plotlist=c(allplots,tilescaleplotlist),ncol = 1)
    }
  }
  
  
  return(list(returnplot=returnplot,tilelegend=legend,derivdata=ggrouttoplotclip,pairs=pairs,pairswithfill=pairswithfill,allplots=allplots,fillplots=fillplots,tilescaleplot=tilescaleplot))
  
}
#####deriv unthresholded plot####

derivplot_unthresholded<-function(derivdat){

derivdat$mean_dff_clip_bin<-dplyr::if_else(derivdat$mean_dff_clip!=0,1,0)
derivdat_diffcliponly<-derivdat[derivdat$mean_dff_clip_bin==1,]
derivdat_diffcliponly$mean_dff_clip_sign<-sign(derivdat_diffcliponly$mean_dff_clip)

if(length(unique(derivdat_diffcliponly$mean_dff_clip_sign))==1){
  if(unique(derivdat_diffcliponly$mean_dff_clip_sign)==1){
    colorone="#c9270e"
    colortwo="grey80"
  }else{
    colorone="#2f42bde6"
    colortwo="grey80"
  }
}else{
  colorone="#2f42bde6"
  colortwo="#c9270e"
}
   
  

ggderivCI<-ggplot(derivdat,aes(x=ages,y=mean_dff))+geom_line(size=.5)+
  geom_ribbon(aes(x=ages,ymin=ci_low,ymax=ci_high),fill="grey55",colour="grey55",alpha=.2)+
  geom_ribbon(data=derivdat_diffcliponly,aes(x=ages,ymin=ci_low,ymax=ci_high,fill=as.factor(mean_dff_clip_sign),colour=as.factor(mean_dff_clip_sign)),alpha=.2)+
  geom_line(data=derivdat_diffcliponly,aes(x=ages,y=mean_dff,colour=as.factor(mean_dff_clip_sign)))+
  scale_colour_manual(values=c(colorone,colortwo))+
  scale_fill_manual(values=c(colorone,colortwo))+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+
  geom_hline(yintercept = 0,linetype="dashed")+
  theme(legend.position = "none")

ggderivCI<-ggfinalize(ggderivCI)+xlab("Age (years)")+ylab("First Derivative of Age Trajectory\n (z-score)")+theme(legend.position = "none",legend.title = element_blank())+
  theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),strip.text.x=element_text(size=10))+
  theme(plot.title = element_text(hjust = 0.5,size=24))

return(ggderivCI)
}
