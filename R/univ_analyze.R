#'Univariate analysis
#'@description Perform univariate analysis.
#'@usage univ_analyze(METBObj, var.equal = FALSE, ispara = FALSE,
#'factor2Col = FALSE, doposthoc = FALSE)
#'@param METBObj METBObj object contains list of data.
#'@param var.equal a logical variable for pairwise comparisons indicating
#'whether to treat the two variances as being equal (see \code{\link[stats:t.test]{stats::t.test()}}).
#'@param ispara a logical variable wheter to perform parametric (TRUE) or non-parametric (FALSE) test. Default is FALSE.
#'@param factor2Col a column number/index of the second factor in the \code{METBObj$inputdata} data frame.
#'This parameter is used in Two-way ANOVA.
#'@param doposthoc a logical variable wheter to perform posthoc test for One-way ANOVA and Two-way ANOVA. Default is FALSE.
#'@details The column number or the column index always begins with 1.
#'Statistical analysis for independent samples will be performed for unequal sample sizes.
#'The friedman.test is the non-parametric alternative to the repeated measures one-way ANOVA and
#'it requires an unreplicated complete block design.
#'@return a list of the following components:
#'
#'stat_summ = a data frame of statistical summary including median, mean, SD (sd) and fold change (fc).
#'If 2 levels, fc=level1/level2, if more than 2 levels, fc=each_level/all_mean. Levels are ordered by characters.
#'
#'p_value = a numeric vector or matrix of p-values.
#'
#'p_adj = a numeric vector or matrix of adjusted p-values.
#'
#'posthoc_data = a list of post-hoc test results. Return empty list for Pairwise comparisons or \code{doposthoc = FALSE}.
#'
#'posthoc_table = a data frame of post-hoc test results. Return empty data frame for Pairwise comparisons or \code{doposthoc = FALSE}.
#'
#'details = a list of analysis details: isParametric, testMethod, pAdjusted, posthocTest.
#'
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@references Dunn, O.J. 1964. Multiple comparisons using rank sums. Technometrics 6:241-252.
#'@references Sokal, R.R. and F.J. Rohlf. 1995. scheirerRayHare. Biometry. 3rd ed. W.H. Freeman, New York.
#'@references http://rcompanion.org/handbook/F_14.html.
#'@seealso \code{\link[stats:t.test]{stats::t.test()}}, \code{\link[stats:wilcox.test]{stats::wilcox.test()}}
#'@examples
#'#sugar_dt = set_input_obj(sugar, 1,2,5)
#'#out=univ_analyze(sugar_dt, factor2Col=3)
#'@export
univ_analyze <- function(METBObj, var.equal = FALSE, ispara = FALSE, factor2Col = FALSE, doposthoc = FALSE){
  #Check argument
  if(sum(is.na(METBObj$X)) > 0){#Data contains missing values
    cat("\nThe data contains missing values. Data was not analyzed.\n")
    return(FALSE)
  }
  #Initialize parameters
  univ_result = list(); methodls = list(); posthoc_data = list(); posthoc_table = data.frame(); posthocTest = "No test"; #Working data
  dat = METBObj$X; isRepeated = METBObj$isRepeated; #working data
  cat("\nChecking essential parameters ...")
  #Check category/factor column
  if (is.factor(METBObj$Y)) {
    F1 = METBObj$Y #working data
    cat('\nFind a factor F1 of ',nlevels(F1),' levels.',sep='')
  }else{
    F1 = as.factor(METBObj$Y) #working data
    cat('\nThe category/factor column is', class(METBObj$Y), '.\nConverting to a factor F1 of ',nlevels(F1),' levels.',sep='')
  }
  #Check no. of levels
  levelrnum = nlevels(F1)
  #Check unequal sample size for F1 repeated samples
  if(isRepeated){#1w-anova
    ID = factor(METBObj$orgID) #working data, don't sort
    nID = nlevels(ID) #no. of unique samples
    nF1 = nlevels(F1) #no. of F1 levels
    # if(nrow(dat) == (nID*nF1)){
    #   isRepeated = TRUE #an unreplicated complete block design
    # }else{
    #   isRepeated = FALSE #assuming independent samples
    #   cat('\nWarning! For F1, data is: \n-not a repeated-measures design or','\n-with replication or',
    #       '\n-with unpaired samples.\nStatistical analysis for independent samples will be performed for One-way ANOVA.',sep='')
    # }
  }
  #Check unequal sample size for F2 repeated samples
  if(factor2Col){#2w-anova
    inputdata = METBObj$inputdata #working data
    factor1Col = METBObj$classCol #working data
    F2 = inputdata[,factor2Col] #2nd factor column
    #Check 2nd factor column
    if (is.factor(F2)) {
      cat('\n\nFind 2nd factor F2 of ',nlevels(F2),' levels.',sep='')
    }else{
      F2 = as.factor(F2)
      cat('\n\nThe 2nd factor column is', class(inputdata[,factor2Col]), '.\nConverting to a factor F2 of ',nlevels(F2),' levels.',sep='')
    }
    ID = factor(METBObj$orgID) #working data, don't sort
    nID = nlevels(ID) #no. of unique samples
    F12 = factor(paste(F1,F2,sep = "-")) #combine factors; no of lv = F1*F2
    nF12 = nlevels(F12) #no. of F2 levels
    if(nrow(dat) == (nID*nF12)){
      isRepeated = TRUE #an unreplicated complete block design
    }else{
      isRepeated = FALSE
      cat('\nWarning! For F2, data is: \n-not a repeated-measures design or','\n-with unequal sample sizes or',
          '\n-with unpaired samples.\nStatistical analysis for independent samples will be performed for Two-way ANOVA.\n',sep='')
    }
  }

  #Univariate analysis
  cat("\n\nPerforming univariate analysis ...")
  stat_summ = get_stat_summary(METBObj) #get stat table
  if(levelrnum == 2 && !(factor2Col)){#t-test
    if(ispara){#parametric test
      t_data = apply(dat, 2, function(x) {t.test(x ~ F1, paired=isRepeated, var.equal=var.equal)})
      p_value = sapply(t_data, function(x) {x$p.value})
      p_adj = p.adjust(p_value, method="fdr") #Benjamini & Hochberg
    }else{#nonparametric test
      t_data = apply(dat, 2, function(x) {wilcox.test(x ~ F1, paired=isRepeated, var.equal=var.equal)})
      p_value = sapply(t_data, function(x) {x$p.value})
      p_adj = p.adjust(p_value, method="fdr") #Benjamini & Hochberg
    }
    cat("\nPairwise comparisons with",t_data[[1]]$method,".\n")
    methodls$isParametric = ispara; methodls$varequal = var.equal; methodls$testMethod = t_data[[1]]$method; methodls$pAdjusted = "Benjamini-Hochberg or FDR"; methodls$posthocTest = posthocTest;
  }
  if((levelrnum>2) && !(factor2Col) && !(isRepeated)){#1w-anova; independent samples
    if(ispara){#parametric test
      aov_data = apply(dat, 2, function(x) {aov(x ~ F1)})
      aov_summ = lapply(aov_data, summary)
      p_value = sapply(aov_summ, function(x) {t(data.frame(unlist(x)))[,"Pr(>F)1"]})
      p_adj = p.adjust(p_value, method="fdr")
      testmethod = "aov(x ~ F1)"
      if(doposthoc){#posthoc test
        posthoc_data = lapply(aov_data, function(x) {data.frame(TukeyHSD(x)$F1)}) #TukeyHSD test
        #Format posthoc table
        posthoc_mx = lapply(posthoc_data, function(x) {t(data.frame(x$p.adj))})
        posthoc_table = data.frame(do.call(rbind, posthoc_mx))
        row.names(posthoc_table) = names(posthoc_mx)
        colnames(posthoc_table) = row.names(posthoc_data[[1]])
        posthocTest = "Tukey multiple comparisons of means"
      }
    }else{#nonparametric test
      aov_data = apply(dat, 2, function(x) {kruskal.test(x ~ F1)})
      p_value = sapply(aov_data, function(x) {x$p.value})
      p_adj = p.adjust(p_value, method="fdr")
      if(doposthoc){#posthoc test
        inst_pkg = NULL
        if(!requireNamespace("FSA", quietly = TRUE)){#check and install required package
          cat("\nMissing the required package 'FSA', trying to install the package ...\n")
          inst_pkg = install_pkgs('FSA')
        }
        if(length(unlist(inst_pkg))){
          posthoc_data = apply(dat, 2, function(x) {pairwise.wilcox.test(x, F1, p.adjust.method = "fdr", paired = FALSE)}) #posthoc; independent samples
          #Format posthoc table
          posthoc_mx = lapply(posthoc_data, function(x) {reshape2::melt(x$p.value, varnames = c("lv1", "lv2"), na.rm = TRUE)})
          posthoc_df = dplyr::bind_rows(posthoc_mx, .id = "var_label")
          posthoc_df$condition = paste0(posthoc_df$lv1,"-",posthoc_df$lv2)
          posthoc_table = reshape2::dcast(posthoc_df, factor(var_label,levels = unique(var_label)) ~ condition, value.var="value")
          row.names(posthoc_table) = posthoc_table[,1]
          posthoc_table = posthoc_table[,-1]
          posthocTest = "Pairwise comparisons using Wilcoxon rank sum exact test; P-adjustment method: fdr"
          cat("\nERROR! Could not install the required package 'FSA'. Post-hoc test for independent samples was calculated.\n")
        }else{
          posthoc_data = apply(dat, 2, function(x) {FSA::dunnTest(x ~ F1, method="bh")$res}) #Dunn test
          #Format posthoc table
          posthoc_mx = lapply(posthoc_data, function(x) {t(data.frame(x$P.adj))})
          posthoc_table = data.frame(do.call(rbind, posthoc_mx))
          row.names(posthoc_table) = names(posthoc_mx)
          colnames(posthoc_table) = posthoc_data[[1]]$Comparison
          posthocTest = "Dunn Kruskal-Wallis multiple comparison; P-adjustment method: fdr"
        }
      }
      testmethod = aov_data[[1]]$method
    }
    cat("\nOne-way ANOVA with",testmethod,".\n")
    methodls$isParametric = ispara; methodls$varequal = FALSE; methodls$testMethod = testmethod; methodls$pAdjusted = "Benjamini & Hochberg or FDR"; methodls$posthocTest = posthocTest;
  }
  if((levelrnum>2) && !(factor2Col) && (isRepeated)){#1w-anova; repeated samples
    if(ispara){#parametric test
      aov_data = apply(dat, 2, function(x) {aov(x ~ F1 + Error(ID))})
      aov_summ = lapply(aov_data, summary)
      p_value = tryCatch({
        sapply(aov_summ, function(x) {t(data.frame(unlist(x)))[,"Error: Within.Pr(>F)1"]})
      },
      error=function(e){
        cat(e$message)
        #message(e)
        cat("\nWarning! Data contains no repeated samples. One-way ANOVA for independent samples was calculated.\n")
        sapply(aov_summ, function(x) {t(data.frame(unlist(x)))[,"Error: ID.Pr(>F)1"]})
      })
      p_adj = p.adjust(p_value, method="fdr")
      if(doposthoc){#posthoc test
        posthoc_data = tryCatch({
          apply(dat, 2, function(x) {pairwise.t.test(x, F1, p.adjust.method = "fdr", paired = TRUE)}) #posthoc; repeated samples
        },
        error=function(e){
          cat(e$message)
          #message(e)
          cat("\nWarning! Data contains no repeated samples or unequal sample sizes. Post-hoc test for independent samples was calculated.\n")
          apply(dat, 2, function(x) {pairwise.t.test(x, F1, p.adjust.method = "fdr", paired = FALSE)}) #posthoc; independent samples
        })
        #Format posthoc table
        posthoc_mx = lapply(posthoc_data, function(x) {reshape2::melt(x$p.value, varnames = c("lv1", "lv2"), na.rm = TRUE)})
        posthoc_df = dplyr::bind_rows(posthoc_mx, .id = "var_label")
        posthoc_df$condition = paste0(posthoc_df$lv1,"-",posthoc_df$lv2)
        posthoc_table = reshape2::dcast(posthoc_df, factor(var_label,levels = unique(var_label)) ~ condition, value.var="value")
        row.names(posthoc_table) = posthoc_table[,1]
        posthoc_table = posthoc_table[,-1]
        posthocTest = paste0(posthoc_data[[1]]$method,"; P-adjustment method:", posthoc_data[[1]]$p.adjust.method)
      }
      testmethod = "aov(x ~ F1 + Error(ID))"
    }else{#nonparametric test
      aov_data = tryCatch({
        apply(dat, 2, function(x) {friedman.test(y=x, groups=F1, blocks=ID)}) #an unreplicated complete block design needed
      },
      error=function(e){
        cat(e$message)
        #message(e)
        cat("\nWarning! Kruskal-Wallis H test for independent samples was calculated.\n")
        isRepeated = FALSE
        apply(dat, 2, function(x) {kruskal.test(x ~ F1)})
      })
      p_value = sapply(aov_data, function(x) {x$p.value})
      p_adj = p.adjust(p_value, method="fdr")
      if(doposthoc){#posthoc test
        posthoc_data = tryCatch({
          apply(dat, 2, function(x) {pairwise.wilcox.test(x, F1, p.adjust.method = "fdr", paired = TRUE)}) #posthoc; repeated samples
        },
        error=function(e){
          cat(e$message)
          #message(e)
          cat("\nWarning! Data contains no repeated samples or unequal sample sizes. Post-hoc test for independent samples was calculated.\n")
          apply(dat, 2, function(x) {pairwise.wilcox.test(x, F1, p.adjust.method = "fdr", paired = FALSE)}) #posthoc; independent samples
        })
        #Format posthoc table
        posthoc_mx = lapply(posthoc_data, function(x) {reshape2::melt(x$p.value, varnames = c("lv1", "lv2"), na.rm = TRUE)})
        posthoc_df = dplyr::bind_rows(posthoc_mx, .id = "var_label")
        posthoc_df$condition = paste0(posthoc_df$lv1,"-",posthoc_df$lv2)
        posthoc_table = reshape2::dcast(posthoc_df, factor(var_label,levels = unique(var_label)) ~ condition, value.var="value")
        row.names(posthoc_table) = posthoc_table[,1]
        posthoc_table = posthoc_table[,-1]
        posthocTest = paste0(posthoc_data[[1]]$method,"; P-adjustment method:", posthoc_data[[1]]$p.adjust.method)
      }
      testmethod = aov_data[[1]]$method
    }
    cat("\nOne-way ANOVA with",testmethod,".\n")
    methodls$isParametric = ispara; methodls$varequal = FALSE; methodls$testMethod = testmethod; methodls$pAdjusted = "Benjamini & Hochberg or FDR"; methodls$posthocTest = posthocTest;
  }
  if((levelrnum>=2) && (factor2Col) && !(isRepeated)){#2w-anova; independent samples
    if(ispara){#parametric test
      aov_data = apply(dat, 2, function(x) {aov(x ~ F1 * F2)})
      aov_summ = lapply(aov_data, summary)
      p_value = sapply(aov_summ, function(x){
        data.frame(x[[1]])[1:3,5]
      })
      row.names(p_value) = c(colnames(inputdata)[c(factor1Col,factor2Col)],paste(colnames(inputdata)[factor1Col],colnames(inputdata)[factor2Col], sep=":"))
      p_adj = t(apply(p_value, 1, function(x){
        p.adjust(x, method="fdr")
      }))
      if(doposthoc){#posthoc test
        posthoc_data = lapply(aov_data, function(x) {TukeyHSD(x)}) #TukeyHSD test
        posthoc_table = data.frame(t(data.frame(sapply(posthoc_data, function(x){c(x$F1[,4],x$F2[,4],x$`F1:F2`[,4])}), check.names = FALSE)), check.names = FALSE)
        colnames(posthoc_table) = unlist(lapply(posthoc_data[[1]],function(x) row.names(x))) #set colname
        posthocTest = "Tukey multiple comparisons of means"
      }
      testmethod = "aov(formula = x ~ F1 * F2)"
    }else{#nonparametric test
      inst_pkg = NULL
      if(!requireNamespace("rcompanion", quietly = TRUE)){#check and install required package
        cat("\nMissing the required package 'rcompanion', trying to install the package ...\n")
        inst_pkg = install_pkgs('rcompanion')
      }
      aov_data = apply(dat, 2, function(x) {rcompanion::scheirerRayHare(y=x,x1=F1,x2=F2, verbose = FALSE)})
      p_value = sapply(aov_data, function(x){
        x[1:3,4]
      })
      row.names(p_value) = c(colnames(inputdata)[c(factor1Col,factor2Col)],paste(colnames(inputdata)[factor1Col],colnames(inputdata)[factor2Col], sep=":"))
      p_adj = t(apply(p_value, 1, function(x){
        p.adjust(x, method="fdr")
      }))
      ph_F1 = apply(dat, 2, function(x){pairwise.wilcox.test(x, F1, p.adjust.method = "fdr", paired = FALSE)}) #posthoc F1; independent samples
      m_F1 = lapply(ph_F1, function(x){ #format matirx
        reshap = reshape2::melt(x$p.value, varnames = c('row', 'col'), na.rm = TRUE)
        matrix(reshap[,3], dimnames=list(paste(reshap$row,reshap$col,sep="-"), "p adj"))
      })
      ph_F2 = apply(dat, 2, function(x) {pairwise.wilcox.test(x, F2, p.adjust.method = "fdr", paired = FALSE)}) #posthoc F2; independent samples
      m_F2 = lapply(ph_F2, function(x){ #format matirx
        reshap = reshape2::melt(x$p.value, varnames = c('row', 'col'), na.rm = TRUE)
        matrix(reshap[,3], dimnames=list(paste(reshap$row,reshap$col,sep="-"), "p adj"))
      })
      ph_F1F2 = apply(dat, 2, function(x) {pairwise.wilcox.test(x, F1:F2, p.adjust.method = "fdr", paired = FALSE)}) #posthoc F1:F2; independent samples
      m_F1F2 = lapply(ph_F1F2, function(x){ #format matirx
        reshap = reshape2::melt(x$p.value, varnames = c('row', 'col'), na.rm = TRUE)
        matrix(reshap[,3], dimnames=list(paste(reshap$row,reshap$col,sep="-"), "p adj"))
      })
      if(doposthoc){#posthoc test
        for(i in 1:ncol(dat)){ #format list
          posthoc_data[[names(ph_F1[i])]] = list(m_F1[[i]],m_F2[[i]],m_F1F2[[i]])
          names(posthoc_data[[i]]) = c('F1','F2','F1:F2')
        }
        posthoc_table = as.data.frame(do.call(rbind,lapply(posthoc_data, function(x){c(x$F1[,1],x$F2[,1],x$`F1:F2`[,1])})))
        colnames(posthoc_table) = unlist(lapply(posthoc_data[[1]],function(x) row.names(x))) #set colname
        posthocTest = paste0(ph_F1[[1]]$method,"; P-adjustment method:", ph_F1[[1]]$p.adjust.method);
      }
      testmethod = "Scheirer Ray Hare test"
    }
    cat("\nTwo-way ANOVA with",testmethod,".\n")
    methodls$isParametric = ispara; methodls$varequal = FALSE; methodls$testMethod = testmethod; methodls$pAdjusted = "Benjamini & Hochberg or FDR"; methodls$posthocTest = posthocTest;
  }
  if((levelrnum>=2) && (factor2Col) && (isRepeated)){#2w-anova; complete-repeated samples; only parametric test
    aov_data = apply(dat, 2, function(x) {aov(x ~ F1 * F2 + Error(ID/(F1*F2)))})
    aov_summ = lapply(aov_data, summary)
    p_value = tryCatch({
      sapply(aov_summ, function(x) {t(data.frame(unlist(x, recursive = F)))[c(10,15,20),1]})
    },
    error=function(e){
      cat(e$message)
      #message(e)
      cat("\nERROR! Data contains no repeated samples or unequal sample sizes. Two-way ANOVA for independent samples was calculated.\n")
      sapply(aov_summ, function(x){data.frame(x[[1]])[1:3,5]})
    })
    row.names(p_value) = c(colnames(inputdata)[c(factor1Col,factor2Col)],paste(colnames(inputdata)[factor1Col],colnames(inputdata)[factor2Col], sep=":"))
    p_adj = t(apply(p_value, 1, function(x){
      p.adjust(x, method="fdr")
    }))
    ph_F1 = apply(dat, 2, function(x){pairwise.t.test(x, F1, p.adjust.method = "fdr", paired = TRUE)}) #posthoc F1; repeated samples
    m_F1 = lapply(ph_F1, function(x){ #format matirx
      reshap = reshape2::melt(x$p.value, varnames = c('row', 'col'), na.rm = TRUE)
      matrix(reshap[,3], dimnames=list(paste(reshap$row,reshap$col,sep="-"), "p adj"))
    })
    ph_F2 = apply(dat, 2, function(x) {pairwise.t.test(x, F2, p.adjust.method = "fdr", paired = TRUE)}) #posthoc F2; repeated samples
    m_F2 = lapply(ph_F2, function(x){ #format matirx
      reshap = reshape2::melt(x$p.value, varnames = c('row', 'col'), na.rm = TRUE)
      matrix(reshap[,3], dimnames=list(paste(reshap$row,reshap$col,sep="-"), "p adj"))
    })
    ph_F1F2 = apply(dat, 2, function(x) {pairwise.t.test(x, F1:F2, p.adjust.method = "fdr", paired = TRUE)}) #posthoc F1:F2; repeated samples
    m_F1F2 = lapply(ph_F1F2, function(x){ #format matirx
      reshap = reshape2::melt(x$p.value, varnames = c('row', 'col'), na.rm = TRUE)
      matrix(reshap[,3], dimnames=list(paste(reshap$row,reshap$col,sep="-"), "p adj"))
    })
    if(doposthoc){#posthoc test
      for(i in 1:ncol(dat)){ #format list
        posthoc_data[[names(ph_F1[i])]] = list(m_F1[[i]],m_F2[[i]],m_F1F2[[i]])
        names(posthoc_data[[i]]) = c('F1','F2','F1:F2')
      }
      posthoc_table = as.data.frame(do.call(rbind,lapply(posthoc_data, function(x){c(x$F1[,1],x$F2[,1],x$`F1:F2`[,1])})))
      colnames(posthoc_table) = unlist(lapply(posthoc_data[[1]],function(x) row.names(x))) #set colname
      posthocTest = paste0(ph_F1[[1]]$method,"; P-adjustment method:", ph_F1[[1]]$p.adjust.method)
    }
    testmethod = "aov(formula = x ~ F1 * F2 + Error(ID/(F1 * F2)))"
    cat("\nTwo-way ANOVA with repeated measures.\n")
    methodls$isParametric = TRUE; methodls$varequal = FALSE; methodls$testMethod = testmethod; methodls$pAdjusted = "Benjamini & Hochberg or FDR"; methodls$posthocTest = posthocTest;
  }
  univ_result$stat_summ = stat_summ; univ_result$p_value = p_value; univ_result$p_adj = p_adj; univ_result$posthoc_data = posthoc_data; univ_result$posthoc_table = posthoc_table; univ_result$details = methodls;
  return(univ_result)
}
