## MUVR functions ##

#'Get VIP ranks
#'@description Get VIP ranks from MUVR object.
#'@usage MUVR_getvip(MVObj, model="min")
#'@param MVObj MVObj object contains list of data.
#'@param model name indicating Which model to use ("min" [default], "mid", or "max").
#'@details Overriding \code{\link[MUVR:getVIP]{MUVR::getVIP()}}.
#'@return data frame with order, name and average rank of features/variables.
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@references Lin Shi, et al. MUVR (2019) \url{https://doi.org/10.1093/bioinformatics/bty710}.
#'@seealso \code{\link[MUVR:getVIP]{MUVR::getVIP()}}
#'@examples
#'#res_muvr = run_muvr(METBObj, method="PLS")
#'#MUVR_getvip(res_muvr)
#'@export
MUVR_getvip <- function(MVObj, model="min") {
  nVar=round(MVObj$nVar[model])
  VIPs=sort(MVObj$VIP[,model])[1:nVar]
  VIPs=data.frame(order=1:nVar,name=names(VIPs),rank=VIPs)
  VIPs$name=as.character(VIPs$name)
  return(VIPs)
}

#'VIP rank plot
#'@description Plot VIP ranking in MUVR object.
#'@usage MUVR_plotvip(MVObj,model="min",plot_title="")
#'@param MVObj MVObj object contains list of data.
#'@param model name indicating Which model to use ("min" [default], "mid", or "max").
#'@param plot_title text indicating plot title.
#'@details Overriding \code{\link[MUVR:plotVIP]{MUVR::plotVIP()}}.
#'@return ggplot object
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@references Lin Shi, et al. MUVR (2019) \url{https://doi.org/10.1093/bioinformatics/bty710}.
#'@seealso \code{\link[MUVR:plotVIP]{MUVR::plotVIP()}}
#'@examples
#'#res_muvr = run_muvr(METBObj, method="PLS")
#'#MUVR_plotvip(res_muvr)
#'@import ggplot2
#'@export
MUVR_plotvip <- function(MVObj, model="min", plot_title="") {
  nModel=ifelse(model=='min',1,ifelse(model=='mid',2,3))
  nFeat=round(MVObj$nVar[nModel])
  VIP=MVObj$VIP[,nModel]
  VIPRep=MVObj$VIPPerRep[[nModel]]
  VIPRep=VIPRep[order(VIP),][1:nFeat,]
  dat=cbind(data.frame(var=row.names(VIPRep)),VIPRep)
  dat$var = factor(dat$var,levels = rev(dat$var))
  dat=reshape2::melt(dat, id.vars=c("var")) #for ggplot
  ggplot(dat, aes(x = var, y = value)) + geom_boxplot(outlier.shape = NA) + ggtitle(plot_title) +
    geom_jitter(width=0.1, alpha=0.2) + theme_light() + coord_flip() + labs(x = "", y = "VIP rank") +
    theme(plot.title = element_text(size = 12), axis.title=element_text(size=10), axis.text=element_text(size=8))
}

#'Validation plot
#'@description Plot validation metric from MUVR object.
#'@usage MUVR_plotval(MVObj,plot_title="")
#'@param MVObj MVObj object contains list of data.
#'@param plot_title text indicating plot title.
#'@details Overriding \code{\link[MUVR:plotVAL]{MUVR::plotVAL()}}.
#'@return ggplot object
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@references Lin Shi, et al. MUVR (2019) \url{https://doi.org/10.1093/bioinformatics/bty710}.
#'@seealso \code{\link[MUVR:plotVAL]{MUVR::plotVAL()}}
#'@examples
#'#res_muvr = run_muvr(METBObj, method="PLS")
#'#MUVR_plotval(res_muvr)
#'@import ggplot2
#'@export
MUVR_plotval <- function(MVObj, plot_title="") {
  VAL=MVObj$VAL$VAL
  metric=MVObj$VAL$metric
  count=as.numeric(colnames(VAL))
  nRep=dim(VAL)[3]
  outplot = ggplot() + theme_classic()
  colors <- c("validations" = "#fec44f", "repetitions" = "#737373", "overall" = "red",
              "min"="#4393c3","mid"="#7fbc41","max"="#c51b7d")
  for (r in 1:nRep) {
    for(o in 1:nrow(VAL[,,r])){
      all_validate = data.frame(count=factor(count),value=VAL[o,,r])
      outplot = outplot+geom_line(data=all_validate,aes(x=count,y=value,group=1,color="validations"),linetype="dashed")
    }
  }
  for (r in 1:nRep) {
    avg_rep = data.frame(count=factor(count),value=colMeans(VAL[,,r]))
    outplot = outplot+geom_line(data=avg_rep,aes(x=count,y=value,group=1,color="repetitions"),size=0.7)
  }
  avg_all = data.frame(count=factor(count),value=apply(VAL,2,mean))
  outplot = outplot+geom_line(data=avg_all,aes(x=count,y=value,group=1,color="overall"),size=0.9)
  xindex = lapply(MVObj$nVar, function(x){#get x-axis index
    tmp = which((rev(count)) == x)
    if(length(tmp) == 0){ tmp=min(which(rev(count)>x))-0.5}
    tmp
  })
  vlines = data.frame(xint=unlist(xindex), value=MVObj$nVar, model=names(MVObj$nVar)) #get v-lines
  for (i in 1:3) {
    outplot = outplot+geom_vline(data=vlines[i,],aes(xintercept = xint, color=model), linetype="dotdash",size=0.8)+
      geom_text(data=vlines[i,], aes(x = xint+0.1, y=round(max(VAL),1), color=model, label=value), angle=90)
  }
  outplot = outplot+labs(x="Number of variables", y=metric) +
    scale_colour_manual(name="",values=colors,label=c("Validations", "Repetitions","Overall","Min-relevant","Mid-relevant","All-relevant")) + ggtitle(plot_title) +
    theme(plot.title = element_text(size = 12), axis.title=element_text(size=10), axis.text=element_text(size=8))
  return(outplot)
}
