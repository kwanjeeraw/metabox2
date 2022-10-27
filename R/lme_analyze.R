#'Linear mixed-effects modeling
#'@description Fit linear mixed-effects models.
#'@usage lme_analyze(METBObj, fix, random)
#'@param METBObj METBObj object contains list of data.
#'@param fix a number or numeric vector indicating the fixed-effects column(s) in the \code{METBObj$inputdata} data frame.
#'@param random text of random-effect formular. See \code{\link[lmm2met:fitLmm]{lmm2met::fitLmm()}}.
#'@details The column number or the column index always begins with 1.
#'Random-effect terms can be in several forms as given in Table 2 of Bates et al. (2015).
#'@return a list of the following components:
#'
#'completeMod = a list of linear mixed models for all features/variables.
#'
#'testRes = a list of chi-square test results for significant fixed effects.
#'
#'testTable = a data frame of Intercept, coefficients and chi-square test results for significant fixed effects.
#'
#'fittedDat = a data frame of fitted data.
#'
#'details = a list of analysis details: testMethod, fix, random.
#'
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@references Bates D, et al.,lme4 (2015). \url{https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf}.
#'@references Wanichthanarak K, et al. lmm2met (2019) . \url{https://doi.org/10.1016/j.csbj.2019.04.009}.
#'@seealso \code{\link[lme4:lmer]{lme4::lmer()}}, \code{\link[stats:drop1]{stats::drop1()}}
#'@examples
#'#sugar_dt = set_input_obj(sugar, 1,2,5)
#'#out=lme_analyze(sugar_dt, fix=c(2:3), random="(1 | subjectID)")
#'@export
lme_analyze <- function(METBObj, fix, random) {
  cat("\nChecking essential parameters ...")
  if(missing(fix) || missing(random)){
    cat("\nError! Argument is missing, with no default.\nData was not analyzed.\n")
    return(FALSE)
  }
  if(!is.numeric(fix) || !is.vector(fix)){
    cat("\nError! Incorrect type of argument, 'fix' is not a number or numeric vector.\nData was not analyzed.\n")
    return(FALSE)
  }
  if(!is.character(random)){
    cat("\nError! Incorrect type of argument, 'random' is not a character.\nData was not analyzed.\n")
    return(FALSE)
  }
  #Check argument
  if(sum(is.na(METBObj$X)) > 0){#Data contains missing values
    cat("\nThe data contains missing values. Data was not analyzed.\n")
    return(FALSE)
  }

  #set default values
  #Initialize parameters
  lme_result = list(); methodls = list();
  dat = METBObj$X; #working data
  cnames = colnames(dat) #keep original colnames
  colnames(dat) = paste0('X',1:ncol(dat)) #working names
  metadat = METBObj$inputdata[,1:METBObj$xCol-1]
  in_dat = cbind(metadat, dat) #for lmer
  completeMod <- list();  testRes <- list();  fit <- metadat;
  fixeff.tab <- function(model){#format testTable
    tb = data.frame()
    for(i in 1:length(model$completeMod)){
      rw = format.row(model$completeMod[[i]],model$testRes[[i]])
      tb = rbind(tb,rw)
    }
    row.names(tb) = names(model$completeMod)
    pchi = ncol(tb)-(nrow(model$testRes[[i]])-1)+1
    value = list(table=tb,pindex=pchi)
    return(value)
  }
  format.row <- function(x,y){#format testTable
    coeff = t(data.frame(summary(x)$coefficients[,1]))
    pchi = t(data.frame(t(y)[4,-1]))
    colnames(pchi) = paste0('Pr(Chi).',rownames(y)[-1])
    cbind(coeff,pchi)
  }

  #Fitting linear mixed-effects models
  cat("\n\nPerforming LME analysis ...")
  fixterm = colnames(metadat)[fix] #get fix terms
  for(i in 1:ncol(dat)){#loop through variables
    val <- colnames(dat)[i]
    cat("Fitting model to",cnames[i],"\n")
    #fit a linear mixed model (LMM)
    completeMod[[val]] <- lme4::lmer(as.formula(paste(val,"~", paste(do.call(c,list(fixterm,random)), collapse = "+"))), data=in_dat)
    #test all fixed effects
    testRes[[val]] <- drop1(completeMod[[val]], test = "Chisq", trace = F)
    fit <- cbind(fit, predict(completeMod[[val]]))
    names(fit)[ncol(fit)] <- cnames[i]
  }
  names(completeMod) = names(testRes) = cnames #rename to original colnames
  methodls$testMethod = "Linear mixed model fit by REML"; methodls$fix = fixterm ; methodls$random = random;
  lme_result$completeMod=completeMod; lme_result$testRes=testRes; lme_result$testTable=fixeff.tab(lme_result)$table;
  lme_result$fittedDat=fit; lme_result$details = methodls;
  return(lme_result)
}
