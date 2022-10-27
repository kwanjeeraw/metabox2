## packMBPLSDA functions ##

#'Scree plot
#'@description Plot eigenvalues against the corresponding principle component number.
#'@usage mbplsda_screeplot(x_input,plot_title="",ptsize=2)
#'@param x_input numeric vectors of eigenvalues.
#'@param plot_title text indicating plot title.
#'@param ptsize a number of geom_point size.
#'@return ggplot object.
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@examples
#'#sugar_dt = set_input_obj(sugar, 1,2,5)
#'#out = mbplsda_analyze(sugar_dt$inputdata[,2:3], list(omics = sugar_dt$X),
#'#nrepet=5, npermut=5, nboot=5)
#'#mbplsda_screeplot(out$result$base_model$eig)
#'@import ggplot2
#'@export
mbplsda_screeplot <- function(x_input, plot_title="", ptsize=2) {
  if(!is.null(x_input)){
    dat = data.frame(pc=seq(length(x_input)), eigen=x_input) #working data
    ggplot(dat, aes(x = pc, y = eigen)) + geom_point(alpha=0.5,size=ptsize) + theme_light() +
      labs(title = plot_title, x = "PC", y = "Eigenvalues", color = "") +
      scale_x_continuous(breaks=dat$pc) +
      theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title=element_text(size=10),
            axis.text=element_text(size=8))
  }else{
    ggplot()+labs(title = "Could not compute model,no plot returned")
  }
}

#'Prediction error rate plot
#'@description Plot prediction error rate against the corresponding principle component number.
#'@usage mbplsda_plottestdim(x_input,plotvalida=TRUE,ptsize=2)
#'@param x_input \code{res_optimal} object result from \code{mbplsda_analyze()}.
#'@param plotvalida a logical variable to plot validation (TRUE) or calibration subset (FALSE).
#'@param ptsize a number of geom_point size.
#'@details Modify from \code{\link[packMBPLSDA:plot_testdim_mbplsda]{packMBPLSDA::plot_testdim_mbplsda()}}.
#'@return ggplot object.
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@references Brandolini-Bunlon M., et al. (2019) Multi-block PLS discriminant analysis. Metabolomics, 15(10):134.
#'@seealso \code{\link[packMBPLSDA:testdim_mbplsda]{packMBPLSDA::testdim_mbplsda()}},\code{\link[packMBPLSDA:plot_testdim_mbplsda]{packMBPLSDA::plot_testdim_mbplsda()}}
#'@examples
#'#sugar_dt = set_input_obj(sugar, 1,2,5)
#'#out = mbplsda_analyze(sugar_dt$inputdata[,2:3], list(omics = sugar_dt$X),
#'#nrepet=5, npermut=5, nboot=5)
#'#mbplsda_plottestdim(out$result$res_optimal)
#'@import ggplot2
#'@export
mbplsda_plottestdim <- function (x_input, plotvalida=TRUE, ptsize=2){
  if(length(x_input)>0){
    if(plotvalida){#plot validation graph
      ERvM.mean = data.frame(globalERvM=x_input$ErrorRateVglobal.max[x_input$ErrorRateVglobal.max$variable=="global","mean"], stringsAsFactors = TRUE)
      ncmp = nrow(ERvM.mean)
      ERvM.ciinf = data.frame(globalERvM=x_input$ErrorRateVglobal.max[x_input$ErrorRateVglobal.max$variable=="global","95CIinf"], stringsAsFactors = TRUE)
      ERvM.cisup = data.frame(globalERvM=x_input$ErrorRateVglobal.max[x_input$ErrorRateVglobal.max$variable=="global","95CIsup"], stringsAsFactors = TRUE)
      ciinf = abs(ERvM.mean[,"globalERvM"] - ERvM.ciinf[,"globalERvM"])
      cisup = abs(ERvM.mean[,"globalERvM"] - ERvM.cisup[,"globalERvM"])
      dat = data.frame(pc=1:ncmp, mean_val=ERvM.mean$globalERvM, ciinf=ciinf, cisup=cisup) #working data
      plot_title = "Mean and 95% CI of the global error rate \n from validation subset"
    }else{#plot calibration graph
      ERcM.mean = data.frame(globalERcM=x_input$ErrorRateCglobal.max[x_input$ErrorRateCglobal.max$variable=="global","mean"], stringsAsFactors = TRUE)
      ncmp = nrow(ERcM.mean)
      ERcM.ciinf = data.frame(globalERcM=x_input$ErrorRateCglobal.max[x_input$ErrorRateCglobal.max$variable=="global","95CIinf"], stringsAsFactors = TRUE)
      ERcM.cisup = data.frame(globalERcM=x_input$ErrorRateCglobal.max[x_input$ErrorRateCglobal.max$variable=="global","95CIsup"], stringsAsFactors = TRUE)
      ciinf = abs(ERcM.mean[,"globalERcM"] - ERcM.ciinf[,"globalERcM"])
      cisup = abs(ERcM.mean[,"globalERcM"] - ERcM.cisup[,"globalERcM"])
      dat = data.frame(pc=1:ncmp, mean_val=ERcM.mean$globalERcM, ciinf=ciinf, cisup=cisup) #working data
      plot_title = "Mean and 95% CI of the global error rate \n from calibration subset"
    }
    ggplot(dat, aes(x=pc, y=mean_val)) +
      #geom_errorbar(aes(ymin=mean_val-ciinf, ymax=mean_val+cisup), alpha=0.5, width=0.05, position=position_dodge(0.1)) +
      #geom_line(alpha=0.5) +  geom_point(alpha=0.5,size=ptsize) +
      geom_pointrange(aes(ymin=mean_val-ciinf, ymax=mean_val+cisup), alpha=0.5) +
      scale_x_continuous(breaks=dat$pc) +
      labs(title = plot_title, x = "PC", y = "Error rate", color = "") + theme_light() +
      theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title=element_text(size=10),
            axis.text=element_text(size=8))
  }else{
    ggplot()+labs(title = "Could not compute model,no plot returned")
  }
}

#'Permutation testing plot
#'@description Plot prediction error rate against the percentage of modified Y-block values in permutation testing.
#'@usage mbplsda_plotpermut(x_input,plot_title="",ptsize=2)
#'@param x_input \code{res_permut} object result from \code{mbplsda_analyze()}.
#'@param plot_title text indicating plot title.
#'@param ptsize a number of geom_point size.
#'@details Modify from \code{\link[packMBPLSDA:plot_permut_mbplsda]{packMBPLSDA::plot_permut_mbplsda()}}.
#'@return ggplot object.
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@references Brandolini-Bunlon M., et al. (2019) Multi-block PLS discriminant analysis. Metabolomics, 15(10):134.
#'@seealso \code{\link[packMBPLSDA:permut_mbplsda]{packMBPLSDA::permut_mbplsda()}},\code{\link[packMBPLSDA:plot_permut_mbplsda]{packMBPLSDA::plot_permut_mbplsda()}}
#'@examples
#'#sugar_dt = set_input_obj(sugar, 1,2,5)
#'#out = mbplsda_analyze(sugar_dt$inputdata[,2:3], list(omics = sugar_dt$X),
#'#nrepet=5, npermut=5, nboot=5)
#'#mbplsda_plotpermut(out$result$res_permut)
#'@import ggplot2
#'@import RColorBrewer
#'@export
mbplsda_plotpermut <- function(x_input, plot_title="", ptsize=2){
  if(length(x_input)>0){
    npermut = nrow(x_input$RV.YYpermut.values)-1
    perm_type = factor(c("No permutation",rep("Permutation", times=npermut)))
    mean_val = x_input$ErrorRateVglobal.max[x_input$ErrorRateVglobal.max$variable=="global","mean"]
    dat = cbind(x_input$prctGlob.Ychange.values,mean_val,perm_type) #working data
    colnames(dat) = c("index","mod_val","mean_val","mod_type")
    grcolors = brewer.pal(name="Set1", n = 8)
    ggplot(dat, aes(x=mod_val, y=mean_val, color=mod_type)) + geom_point(alpha=0.5,size=ptsize) +
      geom_smooth(method=lm, se=FALSE, formula = y ~ x, fullrange=TRUE, size=0.5) + lims(x= c(0,1), y = c(0,1)) +
      labs(title = plot_title, x = "% Modified values", y = "Mean error rate", color = "") + theme_light() + scale_color_manual(values = grcolors) +
      theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title=element_text(size=10),
            axis.text=element_text(size=8))
  }else{
    ggplot()+labs(title = "Could not test model,no plot returned")
  }
}

#'Block importance plot
#'@description Plot block importance and confidence intervals from bootstrapping simulations.
#'@usage mbplsda_plotboot_bipc(x_input,plot_title="",ptsize=2)
#'@param x_input \code{res_boot} object result from \code{mbplsda_analyze()}.
#'@param plot_title text indicating plot title.
#'@param ptsize a number of geom_point size.
#'@details Modify from \code{\link[packMBPLSDA:plot_boot_mbplsda]{packMBPLSDA::plot_boot_mbplsda()}}.
#'@return ggplot object.
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@references Brandolini-Bunlon M., et al. (2019) Multi-block PLS discriminant analysis. Metabolomics, 15(10):134.
#'@seealso \code{\link[packMBPLSDA:boot_mbplsda]{packMBPLSDA::boot_mbplsda()}},\code{\link[packMBPLSDA:plot_boot_mbplsda]{packMBPLSDA::plot_boot_mbplsda()}}
#'@examples
#'#sugar_dt = set_input_obj(sugar, 1,2,5)
#'#out = mbplsda_analyze(sugar_dt$inputdata[,2:3], list(omics = sugar_dt$X),
#'#nrepet=5, npermut=5, nboot=5)
#'#mbplsda_plotboot_bipc(out$result$res_boot)
#'@import ggplot2
#'@import RColorBrewer
#'@export
mbplsda_plotboot_bipc <- function(x_input, plot_title="", ptsize=2){
  if(length(x_input)>0){
    dat = x_input$bipc #working data
    if(nrow(dat)>1){#two or more data sets
      ciinf = abs(dat$mean - dat$`95CIinf`)
      cisup = abs(dat$mean - dat$`95CIsup`)
      dat = cbind(dat, ciinf=ciinf, cisup=cisup) #working data
      if(nrow(dat) <= 8){
        grcolors = brewer.pal(name="Set1", n = 8)
      }else{
        grcolors = colorRampPalette(brewer.pal(name="Set1", n = 8))(nrow(dat))
      }
      ggplot(dat, aes(x=blocks, y=mean, color=blocks)) +
        #geom_errorbar(aes(ymin=mean-ciinf, ymax=mean+cisup, color=blocks), alpha=0.5, width=0.05, position=position_dodge(0.1)) +
        #geom_line(alpha=0.5) + geom_point(aes(color=blocks),alpha=0.5,size=ptsize) +
        geom_pointrange(aes(ymin=mean-ciinf, ymax=mean+cisup, color=blocks), alpha=0.5) +
        scale_color_manual(name="Legend",values = grcolors) +
        labs(title = plot_title, x = "Data sets", y = "BIPc") + theme_light() +
        theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title=element_text(size=10),
              axis.text=element_text(size=8), axis.text.x = element_text(angle=90,hjust=1,vjust=0.2))
    }else{#one data set
      cat("\nThere is only one data set. The plot can't be generated.\n")
      NULL
    }
  }else{
    ggplot()+labs(title = "Could not compute model,no plot returned")
  }
}

#'Variable importance plot
#'@description Plot variable importance (VIP) and confidence intervals from bootstrapping simulations.
#'@usage mbplsda_plotboot_vipc(x_input,plot_title="",ptsize=2,propbestvar=0.5)
#'@param x_input \code{res_boot} object result from \code{mbplsda_analyze()}.
#'@param plot_title text indicating plot title.
#'@param ptsize a number of geom_point size.
#'@param propbestvar a number between 0 and 1 (Default = 0.5) indicating the top n-percent of variables with the best VIP values to plot.
#'@details Modify from \code{\link[packMBPLSDA:plot_boot_mbplsda]{packMBPLSDA::plot_boot_mbplsda()}}.
#'@return ggplot object.
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@references Brandolini-Bunlon M., et al. (2019) Multi-block PLS discriminant analysis. Metabolomics, 15(10):134.
#'@seealso \code{\link[packMBPLSDA:boot_mbplsda]{packMBPLSDA::boot_mbplsda()}},\code{\link[packMBPLSDA:plot_boot_mbplsda]{packMBPLSDA::plot_boot_mbplsda()}}
#'@examples
#'#sugar_dt = set_input_obj(sugar, 1,2,5)
#'#out = mbplsda_analyze(sugar_dt$inputdata[,2:3], list(omics = sugar_dt$X),
#'#nrepet=5, npermut=5, nboot=5)
#'#mbplsda_plotboot_vipc(out$result$res_boot)
#'@import ggplot2
#'@import RColorBrewer
#'@export
mbplsda_plotboot_vipc <- function(x_input, plot_title="", ptsize=2, propbestvar=0.5){
  if(length(x_input)>0){
    ## higher than the 1-propbestvar quantile
    dat = x_input$vipc[order(x_input$vipc$mean, decreasing=FALSE),] #working data
    vipcBest = dat[dat$mean>quantile(dat$mean, probs=(1-propbestvar), na.rm=TRUE),] #working data
    ciinf = abs(vipcBest$mean - vipcBest$`95CIinf`)
    cisup = abs(vipcBest$mean - vipcBest$`95CIsup`)
    vipcBest = cbind(vipcBest, ciinf=ciinf, cisup=cisup)
    vipcBest$variables = factor(vipcBest$variables,levels = vipcBest$variables)
    vipcBest$block = factor(vipcBest$block)
    if(nlevels(vipcBest$block) <= 8){
      grcolors = brewer.pal(name="Set1", n = 8)
    }else{
      grcolors = colorRampPalette(brewer.pal(name="Set1", n = 8))(nlevels(vipcBest$block))
    }
    ggplot(vipcBest, aes(x=variables, y=mean, color=block)) +
      #geom_errorbar(aes(ymin=mean-ciinf, ymax=mean+cisup), alpha=0.5, width=0.1, position=position_dodge(0.1)) +
      #geom_line(alpha=0.5) +  geom_point(alpha=0.5,size=ptsize) +
      geom_pointrange(aes(ymin=mean-ciinf, ymax=mean+cisup, color=block), alpha=0.5) +
      scale_color_manual(name="Legend",values = grcolors) +
      labs(title = plot_title, x = "", y = "VIPc", color = "") + theme_light() +
      theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title=element_text(size=10),
            axis.text=element_text(size=8), axis.text.x = element_text(angle=90,hjust=1,vjust=0.2))
  }else{
    ggplot()+labs(title = "Could not compute model,no plot returned")
  }
}

#'Loading plot by PC
#'@description Plot loadings and confidence intervals from bootstrapping simulations by a principle component (pc).
#'@usage mbplsda_plotboot_loading(x_input,pc=1,plot_title="",ptsize=2,propbestvar=0.5)
#'@param x_input \code{res_boot} object result from \code{mbplsda_analyze()}.
#'@param pc a number indicating a principal component to plot. 1st PC is plotted by default.
#'@param plot_title text indicating plot title.
#'@param ptsize a number of geom_point size.
#'@param propbestvar a number between 0 and 1 (Default = 0.5) indicating the top n-percent of variables with the best loading values to plot.
#'@details Modify from \code{\link[packMBPLSDA:plot_boot_mbplsda]{packMBPLSDA::plot_boot_mbplsda()}}.
#'@return ggplot object.
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@references Brandolini-Bunlon M., et al. (2019) Multi-block PLS discriminant analysis. Metabolomics, 15(10):134.
#'@seealso \code{\link[packMBPLSDA:boot_mbplsda]{packMBPLSDA::boot_mbplsda()}},\code{\link[packMBPLSDA:plot_boot_mbplsda]{packMBPLSDA::plot_boot_mbplsda()}}
#'@examples
#'#sugar_dt = set_input_obj(sugar, 1,2,5)
#'#out = mbplsda_analyze(sugar_dt$inputdata[,2:3], list(omics = sugar_dt$X),
#'#nrepet=5, npermut=5, nboot=5)
#'#mbplsda_plotboot_loading(out$result$res_boot)
#'@import ggplot2
#'@import RColorBrewer
#'@export
mbplsda_plotboot_loading <- function(x_input, pc=1, plot_title="", ptsize=2, propbestvar=0.5){
  if(length(x_input)>0){
    #Check argument
    if (length(x_input$faX) < pc){
      stop("argument 'pc' is not valid, the number of principle components is less than 'pc'.")
    }
    matfaX = x_input$faX[[pc]] #working data
    matfaX = matfaX[order(matfaX$mean, decreasing=FALSE),]
    faXBest = matfaX[abs(matfaX$mean)>quantile(abs(matfaX$mean), probs=(1-propbestvar), na.rm=TRUE),]
    ciinf = abs(faXBest$mean - faXBest$`95CIinf`)
    cisup = abs(faXBest$mean - faXBest$`95CIsup`)
    faXBest = cbind(faXBest, ciinf=ciinf, cisup=cisup)
    faXBest$variables = factor(faXBest$variables,levels = faXBest$variables)
    faXBest$block = factor(faXBest$block)
    if(nlevels(faXBest$block) <= 8){
      grcolors = brewer.pal(name="Set1", n = 8)
    }else{
      grcolors = colorRampPalette(brewer.pal(name="Set1", n = 8))(nlevels(faXBest$block))
    }
    ggplot(faXBest, aes(x=variables, y=mean, color=block)) +
      #geom_errorbar(aes(ymin=mean-ciinf, ymax=mean+cisup), alpha=0.5, width=0.1, position=position_dodge(0.1)) +
      #geom_line(alpha=0.5) +  geom_point(alpha=0.5,size=ptsize) +
      geom_pointrange(aes(ymin=mean-ciinf, ymax=mean+cisup, color=block), alpha=0.5) +
      scale_color_manual(name="Legend",values = grcolors) +
      labs(title = plot_title, x = "", y = paste("Loading",pc), color = "") + theme_light() +
      theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title=element_text(size=10),
            axis.text=element_text(size=8), axis.text.x = element_text(angle=90,hjust=1,vjust=0.2))
  }else{
    ggplot()+labs(title = "Could not compute model,no plot returned")
  }
}

#'Regression coefficient plot by factor level
#'@description Plot regression coefficients and confidence intervals from bootstrapping simulations by a factor level (lv).
#'@usage mbplsda_plotboot_coeff(x_input,lv=1,plot_title="",ptsize=2,propbestvar=0.5)
#'@param x_input \code{res_boot} object result from \code{mbplsda_analyze()}.
#'@param lv a number indicating a factor level to plot. 1st level is plotted by default.
#'@param plot_title text indicating plot title.
#'@param ptsize a number of geom_point size.
#'@param propbestvar a number between 0 and 1 (Default = 0.5) indicating the top n-percent of variables with the best coefficient values to plot.
#'@details Modify from \code{\link[packMBPLSDA:plot_boot_mbplsda]{packMBPLSDA::plot_boot_mbplsda()}}.
#'@return ggplot object.
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@references Brandolini-Bunlon M., et al. (2019) Multi-block PLS discriminant analysis. Metabolomics, 15(10):134.
#'@seealso \code{\link[packMBPLSDA:boot_mbplsda]{packMBPLSDA::boot_mbplsda()}},\code{\link[packMBPLSDA:plot_boot_mbplsda]{packMBPLSDA::plot_boot_mbplsda()}}
#'@examples
#'#sugar_dt = set_input_obj(sugar, 1,2,5)
#'#out = mbplsda_analyze(sugar_dt$inputdata[,2:3], list(omics = sugar_dt$X),
#'#nrepet=5, npermut=5, nboot=5)
#'#mbplsda_plotboot_coeff(out$result$res_boot)
#'@import ggplot2
#'@import RColorBrewer
#'@export
mbplsda_plotboot_coeff <- function(x_input, lv=1, plot_title="", ptsize=2, propbestvar=0.5){
  if(length(x_input)>0){
    #Check argument
    if (length(x_input$XYcoef) < lv){
      stop("argument 'lv' is not valid, the number of factor levels is less than 'lv'.")
    }
    matXYcoef = x_input$XYcoef[[lv]] #working data
    matXYcoef = matXYcoef[order(matXYcoef$mean, decreasing=FALSE),]
    XYcoefBest = matXYcoef[abs(matXYcoef$mean)>quantile(abs(matXYcoef$mean), probs=(1-propbestvar), na.rm=TRUE),]
    ciinf = abs(XYcoefBest$mean - XYcoefBest$`95CIinf`)
    cisup = abs(XYcoefBest$mean - XYcoefBest$`95CIsup`)
    XYcoefBest = cbind(XYcoefBest, ciinf=ciinf, cisup=cisup)
    XYcoefBest$variables = factor(XYcoefBest$variables,levels = XYcoefBest$variables)
    XYcoefBest$block = factor(XYcoefBest$block)
    if(nlevels(XYcoefBest$block) <= 8){
      grcolors = brewer.pal(name="Set1", n = 8)
    }else{
      grcolors = colorRampPalette(brewer.pal(name="Set1", n = 8))(nlevels(XYcoefBest$block))
    }
    ggplot(XYcoefBest, aes(x=variables, y=mean, color=block)) +
      #geom_errorbar(aes(ymin=mean-ciinf, ymax=mean+cisup), alpha=0.5, width=0.1, position=position_dodge(0.1)) +
      #geom_line(alpha=0.5) +  geom_point(alpha=0.5,size=ptsize) +
      geom_pointrange(aes(ymin=mean-ciinf, ymax=mean+cisup, color=block), alpha=0.5) +
      scale_color_manual(name="Legend",values = grcolors) +
      labs(title = plot_title, x = "", y = paste0("Regression coefficients\n(",names(x_input$XYcoef)[[lv]],")"), color = "") + theme_light() +
      theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title=element_text(size=10),
            axis.text=element_text(size=8), axis.text.x = element_text(angle=90,hjust=1,vjust=0.2))
  }else{
    ggplot()+labs(title = "Could not compute model,no plot returned")
  }
}
