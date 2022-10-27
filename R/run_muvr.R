#'Run MUVR
#'@description Perform ML analysis using MUVR package.
#'@usage run_muvr(METBObj, scale=FALSE, nRep=3, varRatio=0.75, partitionSize=0.7,
#'method="PLS", DA=TRUE, fitness="AUROC")
#'@param METBObj METBObj object contains list of data.
#'@param scale a logical variable indicating whether the variables should be scaled to have unit variance before PLS modelling. Default is FALSE.
#'@param nRep number of repetitions of double CV.
#'@param varRatio ratio of variables to include in inner loop iterations. Default is 0.75.
#'@param partitionSize number of data holding for training set. Default is 0.7.
#'@param method ML method. Choose one from the list: PLS, RF. Default is PLS.
#'@param DA a logical variable for performing classification. Default is TRUE. For performing regression, DA=FALSE.
#'@param fitness fitness function for model tuning and evaluation. Choose AUROC (default) or MISS for classification; or choose RMSEP for regression.
#'@details Based on MUVR version 0.0.975. (2021-09-21). The \code{partitionSize} is for calculating \code{nOuter}.
#'@return MVObject object.
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@references Lin Shi, et al. MUVR (2019) \url{https://doi.org/10.1093/bioinformatics/bty710}.
#'@seealso \code{\link[MUVR:MUVR]{MUVR::MUVR()}}, \code{\link[MUVR:vectSamp]{MUVR::vectSamp()}}
#'@examples
#'#res_muvr = run_muvr(METBObj, method="PLS")
#'#MUVR::plotVAL(res_muvr) #need more distinct colors
#'#MUVR::plotMV(res_muvr, model = 'min')
#'#MUVR::plotStability(res_muvr, model = 'min')
#'#MUVR::plotVIP(res_muvr, model = 'min') #error if nrep=1
#'#MUVR::getVIP(res_muvr, model = 'min')
#'@export
run_muvr <- function(METBObj, scale=FALSE, nRep=3, varRatio=0.75, partitionSize=0.7, method="PLS", DA=TRUE, fitness="AUROC"){
  #Check argument
  tmparg_method <- try(method <- match.arg(method, c("PLS","RF"), several.ok = FALSE), silent = TRUE)
  if (class(tmparg_method) == "try-error") {
    cat("\nERROR! Argument 'method' is not valid, choose one from the list: PLS,RF.\n")
    return(list())
  }
  tmparg_fitness <- try(fitness <- match.arg(fitness, c("AUROC","MISS","RMSEP"), several.ok = FALSE), silent = TRUE)
  if (class(tmparg_fitness) == "try-error") {
    cat("\nERROR! Argument 'method' is not valid, choose one from the list: AUROC,MISS,RMSEP.\n")
    return(list())
  }
  #Check argument
  if(sum(is.na(METBObj$X)) > 0){#Data contains missing values
    cat("\nThe data contains missing values. Data was not analyzed.\n")
    return(FALSE)
  }

  ###Parallel processing
  # library(foreach)
  # if (parallel) "%doVersion%"=get("%dopar%") else "%doVersion%"=get("%do%")


  ###INITIALIZE PARAMETERS
  scale=scale
  modReturn=TRUE
  methParam=list(compMax=5,robust=0.05,ntreeIn=150,ntreeOut=300,mtryMaxIn=150,oneHot=FALSE,NZV=FALSE,rfMethod="randomForest")
  #methParams
  if (method=='PLS') {
    methParam$NZV <- TRUE
  }
  ML=FALSE
  parallel=FALSE
  if (method == 'RF') library(randomForest)
  library(pROC)
  #files.sources = list.files("MUVR",pattern="*.R", full.names=TRUE)
  #sapply(files.sources, source)

  ###INITIALIZE VARIABLES
  start.time=proc.time()[3]
  modelReturn=list(call=match.call())
  ID = METBObj$ID #ID
  orgID = METBObj$orgID #orgID
  unik = METBObj$unik #unik
  unikID = METBObj$unikID #unikID
  classCol = METBObj$classCol #classCol
  classLs = METBObj$Y #classLs-list of classes for DA only
  X = METBObj$X #X
  Y = METBObj$Y #Y

  ###START CALCULATION
  #Find nOuter
  if (DA & identical(unikID,ID)) {##Independent samples => Classification
    cat('\nClassification; Assuming samples are independent.')
    nTest = (table(classLs))-(round(table(classLs) * partitionSize)) #finding test size
    nOuter = min(round((table(classLs))/nTest)) #finding no. of folds
  } else{##Repeated samples => Regression|Classification; Independent samples => Regression
    cat('\nClassification or Regression; Assuming samples are dependent.')
    nTest = (length(unikID))-(round(length(unikID) * partitionSize)) #finding test size
    nOuter = round(length(unikID)/nTest) #finding no. of folds
  }
  if(nOuter <= 1){
    nOuter = 2 #default no. of folds
    nInner = nOuter #no. of inner folds
  }else{
    nInner = nOuter-1 #default no. of inner folds
  }

  ###EXECUTE MUVR
  #Check input data
  if (length(dim(X))!=2) stop('\nWrong format of X matrix.\n')
  if (is.null(colnames(X))) stop('\nNo column names in X matrix.\n')
  if (nrow(X)!=length(Y)) stop('\nMust have same nSamp in X and Y.\n')
  if (any(is.na(X)) | any(is.na(Y))) stop('\nNo missing values allowed in X or Y data.\n')
  if (!is.null(dim(Y))) stop('\nY is not a vector.\n')

  #Format input data
  X = as.matrix(X)
  if (is.character(Y)) Y=factor(Y)
  if (is.factor(Y)) {
    cat('\nY is factor -> Classification (',length(unique(Y)),' classes)',sep='')
    DA=TRUE
  }
  if (is.numeric(Y) & DA) {
    Y=as.factor(Y)
    cat('\nDA=TRUE -> Y as factor -> Classification (',length(unique(Y)),' classes)',sep='')
  }
  #Remove nearZeroVariance variables
  if (methParam$NZV) {
    nzv=MUVR::nearZeroVar(X)
    if (length(nzv$Position)>0) {
      modelReturn$nzv=colnames(X)[nzv$Position]
      X=X[,-nzv$Position]
      cat('\n',length(nzv$Position),'variables with near zero variance detected -> removed from X and stored under $nzv')
    }
  }

  #Number of samples and variables
  nSamp=nrow(X)
  nVar=nVar0=ncol(X)

  #Internal modelling parameters
  if (method=='PLS') {
    if (nVar<methParam$compMax) methParam$compMax <- nVar # nCompMax cannot be larger than number of variables
  }
  #fitness
  if (missing(fitness)) {
    if (DA) {
      fitness='MISS'
      cat('\nMissing fitness -> MISS')
    } else {
      fitness='RMSEP'
      cat('\nMissing fitness -> RMSEP')
    }
  }

  #Store input data in list
  InData=list(X=X,Y=Y,ID=ID,scale=scale,nRep=nRep,nOuter=nOuter,nInner=nInner,varRatio=varRatio,DA=DA,fitness=fitness,method=method,methParam=methParam,ML=ML,parallel=parallel)

  #Allocate variables for outer loop
  if (DA) {
    if(nOuter>min(table(Y))) stop('\nnOuter is larger than the smallest group size. Consider lowering your nOuter to the smallest group size.',call.=TRUE)
    unikY=Y[unik]  #Counterintuitive, but needed for groupings by Ynames
    Ynames=sort(unique(Y))  #Find groups
    groups=length(Ynames) #Number of groups
    groupID=list()  #Allocate list for indices of groups
    for (g in 1:groups) {
      groupID[[g]]=unikID[unikY==Ynames[g]]  #Find indices per group
    }
    yPredMin=yPredMid=yPredMax=array(dim=c(length(Y),length(levels(Y)),nRep),dimnames=list(ID,levels(Y),paste('Rep',1:nRep,sep='')))
    yPredMinR=yPredMidR=yPredMaxR=matrix(nrow=length(Y),ncol=length(levels(Y)),dimnames=list(ID,levels(Y)))
  } else {
    yPredMin=yPredMid=yPredMax=matrix(nrow=length(Y),ncol=nRep,dimnames=list(ID,paste('Rep',1:nRep,sep='')))
    yPredMinR=yPredMidR=yPredMaxR=numeric(length(Y))
  }

  #Allocate response vectors and matrices for var's, nComp and VIP ranks over repetitions
  missRep=numeric(nRep)
  names(missRep)=paste(rep('rep',nRep),1:nRep,sep='')
  varRepMin=varRepMid=varRepMax=nCompRepMin=nCompRepMid=nCompRepMax=missRep
  nCompSegMin=nCompSegMid=nCompSegMax=matrix(nrow=nRep,ncol=nOuter,dimnames=list(paste('repetition',1:nRep,sep=''),paste('segment',1:nOuter,sep='')))
  VIPRepMin=VIPRepMid=VIPRepMax=matrix(data=nVar0,nrow=nVar0,ncol=nRep,dimnames=list(colnames(X),paste(rep('rep',nRep),1:nRep,sep='')))

  #Number of selected variables for iterations in inner loops
  var=numeric()
  cnt=0
  while (nVar>1) {
    cnt=cnt+1
    var=c(var,nVar)
    nVar=floor(varRatio*nVar)
  }

  #Allocate array for validation results
  VAL=array(dim=c(nOuter,cnt,nRep),dimnames=list(paste('outSeg',1:nOuter,paste=''),var,paste(rep('rep',nRep),1:nRep,sep='')))

  ##Start loops
  reps = list() #output from each repetition
  for (r in 1:nRep){##REPEATE DATA PARTITIONS
  #reps=foreach(r=1:nRep, .packages=packs, .export=exports) %doVersion% {
    if (modReturn) outMod=list()
    cat('\n','   Repetition ',r,' of ',nRep,':',sep='')
    #Random sampling outer sets
    if (DA & identical(unikID,ID)) {
      groupTest=list()  #Allocate list for samples within group
      for (gT in 1:groups) {
        groupTest[[gT]]=MUVR::vectSamp(groupID[[gT]],n=nOuter)  #Draw random samples within group
      }
      allTest=groupTest[[1]] #Add 1st groups to 'Master' sample of all groups
      for (gT in 2:groups) {  #Add subsequent groups
        allTest=allTest[order(sapply(allTest,length))]
        for (aT in 1:nOuter) {
          allTest[[aT]]=sort(c(allTest[[aT]],groupTest[[gT]][[aT]]))
        }
      }
    } else {
      allTest=MUVR::vectSamp(unikID,n=nOuter)
    }

    #Allocate variables for outer loop
    nCompOutMax=numeric(nOuter)
    names(nCompOutMax)=paste(rep('outSeg',nOuter),1:nOuter,sep='')
    varOutMin=varOutMid=varOutMax=nCompOutMin=nCompOutMid=nCompOutMax
    VIPOutMin=VIPOutMid=VIPOutMax=matrix(data=nVar0,nrow=nVar0,ncol=nOuter,dimnames=list(colnames(X),paste(rep('outSeg',nOuter),1:nOuter,sep='')))
    VALRep=matrix(nrow=nOuter,ncol=cnt)

    for (i in 1:nOuter){##OUTER LOOP
      cat('\n Segment ',i,' (variables):',sep='') # Counter
      #Draw out test set
      testID=allTest[[i]] #Draw out segment = holdout set BASED ON UNIQUE ID
      testIndex=ID%in%testID #Boolean for samples included
      xTest=X[testIndex,]
      yTest=Y[testIndex]
      inID=unikID[!unikID%in%testID]  #IDs not in test set
      if (DA & identical(unikID,ID)) inY=unikY[!unikID%in%testID]  #Counterintuitive, but needed for grouping by Ynames

      #Allocate variables for inner loop
      missIn=berIn=aucIn=rmsepIn=PRESSIn=nCompIn=matrix(nrow=nInner,ncol=cnt,dimnames=list(paste(rep('inSeg',nInner),1:nInner,sep=''),var))
      VIPInner=array(data=nVar0,dim=c(nVar0,cnt,nInner),dimnames=list(colnames(X),var,paste(rep('inSeg',nInner),1:nInner,sep='')))

      #Set variables for modeling
      incVar = colnames(X)

      for (count in 1:cnt) {##VARIABLE REDUCTION
        #Set parameters
        nVar = var[count]
        cat(nVar)
        if (method=='PLS') comp=min(c(nVar,methParam$compMax))
        if (method=='RF') {
          mtryIn=ifelse(DA,min(c(methParam$mtryMaxIn,floor(sqrt(nVar)))),min(c(methParam$mtryMaxIn,floor(nVar/3))))
          mtryIn=max(c(2,mtryIn))
        }
        #Random sampling inner sets
        if (DA & identical(unikID,ID)) {
          groupIDVal=list()
          for (g in 1:groups) {
            groupIDVal[[g]]=inID[inY==Ynames[g]]  #Find indices per group
          }
          groupVal=list()  #Allocate list for samples within group
          for (gV in 1:groups) {
            groupVal[[gV]]=MUVR::vectSamp(groupIDVal[[gV]],n=nInner)  #Draw random samples within group
          }
          allVal=groupVal[[1]] #Add 1st groups to 'Master' sample of all groups
          for (gV in 2:groups) {  #Add subsequent groups
            allVal=allVal[order(sapply(allVal,length))]
            for (aV in 1:nInner) {
              allVal[[aV]]=sort(c(allVal[[aV]],groupVal[[gV]][[aV]]))
            }
          }
        } else {
          allVal=MUVR::vectSamp(inID,n=nInner)
        }

        for (j in 1:nInner) {##INNER LOOP
          cat('.') # Counter
          valID=allVal[[j]] #Draw out segment = validation set
          valIndex=ID%in%valID
          xVal=X[valIndex,]
          xVal=subset(xVal,select=incVar)
          yVal=Y[valIndex]
          trainID=inID[!inID%in%valID]
          trainIndex=ID%in%trainID #Define Training segment
          xTrain=X[trainIndex,]
          xTrain=subset(xTrain,select=incVar)
          yTrain=Y[trainIndex]

          ##Run inner model
          if (method=='PLS') {
            inMod=MUVR::plsInner(xTrain,yTrain,xVal,yVal,DA,fitness,comp,scale=scale)
            nCompIn[j,count]=inMod$nComp
          } else {
            inMod=MUVR::rfInner(xTrain,yTrain,xVal,yVal,DA,fitness,mtry=mtryIn,ntree=methParam$ntreeIn,method=methParam$rfMethod)
          }
          #Store fitness metric
          if (fitness=='MISS') {
            missIn[j,count]=inMod$miss
          } else if (fitness=='BER') {
            berIn[j,count]=inMod$ber
          } else if (fitness=='AUROC') {
            aucIn[j,count]=inMod$auc
          } else {
            rmsepIn[j,count]=inMod$rmsep
            PRESSIn[j,count]=(inMod$rmsep^2)*length(yVal)
          }
          #Store VIs
          VIPInner[match(names(inMod$vi),rownames(VIPInner)),count,j]=inMod$vi

        }##END INNER LOOP

        ##Average inner VIP ranks
        VIPInAve=apply(VIPInner[,count,],1,mean)
        if (count<cnt) {
          #incVar=names(sort(VIPInAve))[1:var[count+1]] #Set variables for modeling
          incVar <- names(VIPInAve[order(VIPInAve)])[1:var[count + 1]] #Extract the names of the variables kept for the next iteration
        }
      }##END VARIABLE REDUCTION

      #Per outer loop: Average inner loop variables, nComp and VIP ranks
      if (fitness=='AUROC') {
        fitRank=colMeans(-aucIn)
        VALRep[i,]=colMeans(aucIn)
      } else if (fitness=='MISS') {
        fitRank=VALRep[i,]=colSums(missIn)
      } else if (fitness=='BER') {
        fitRank=VALRep[i,]=colMeans(berIn)
      }else {
        fitRank=colMeans(rmsepIn)
        VALRep[i,]=sqrt(colSums(PRESSIn)/sum(!testIndex))
      }
      fitRank=(fitRank-min(fitRank))/abs(diff(range(fitRank))) #Rescale fitRank to range 0-1
      if(all(is.nan(fitRank))) fitRank=rep(0,cnt) #If all VAL have the same value -> reset all fitRank to 0
      minIndex=max(which(fitRank <= methParam$robust)) #Get index of min number of val
      maxIndex=min(which(fitRank <= methParam$robust)) #Get index of max number of val
      varOutMin[i]=var[minIndex]
      varOutMax[i]=var[maxIndex]
      varOutMid[i]=round(exp(mean(log(c(var[minIndex],var[maxIndex])))))
      midIndex=which.min(abs(var-varOutMid[i])) #Get index of mid number of val
      if (method=='PLS') {#Average no. of comp
        nCompOutMin[i]=round(mean(nCompIn[,minIndex]))
        nCompOutMid[i]=round(mean(nCompIn[,midIndex]))
        nCompOutMax[i]=round(mean(nCompIn[,maxIndex]))
      }
      #Average inner VIP ranks from all count
      VIPOutMin[,i]=apply(VIPInner[,minIndex,],1,mean)
      VIPOutMid[,i]=apply(VIPInner[,midIndex,],1,mean)
      VIPOutMax[,i]=apply(VIPInner[,maxIndex,],1,mean)

      #Build outer model for min mid and max of nComp and predict YTEST
      xIn=X[!testIndex,]
      yIn=Y[!testIndex]
      incVarMin=rownames(VIPOutMin)[rank(VIPOutMin[,i])<=varOutMin[i]]
      incVarMid=rownames(VIPOutMid)[rank(VIPOutMid[,i])<=varOutMid[i]]
      incVarMax=rownames(VIPOutMax)[rank(VIPOutMax[,i])<=varOutMax[i]]

      if (method=='PLS'){
        #Min model
        if (DA) {
          plsOutMin=MUVR::plsda(subset(xIn,select=incVarMin),yIn,ncomp=nCompOutMin[i],near.zero.var=TRUE,scale=scale)
        }else{
          plsOutMin=MUVR::pls(subset(xIn,select=incVarMin),yIn,ncomp=nCompOutMin[i],near.zero.var=TRUE,scale=scale)
        }
        xTestMin=subset(xTest,select=incVarMin)
        yPredMinR[testIndex]=predict(plsOutMin,newdata=xTestMin,scale=scale)$predict[,,nCompOutMin[i]]  #
        #Mid model
        if (DA) {
          plsOutMid=MUVR::plsda(subset(xIn,select=incVarMid),yIn,ncomp=nCompOutMid[i],near.zero.var=TRUE,scale=scale)
        }else{
          plsOutMid=MUVR::pls(subset(xIn,select=incVarMid),yIn,ncomp=nCompOutMid[i],near.zero.var=TRUE,scale=scale)
        }
        xTestMid=subset(xTest,select=incVarMid)
        yPredMidR[testIndex]=predict(plsOutMid,newdata=xTestMid,scale=scale)$predict[,,nCompOutMid[i]]  #
        #Max model
        if (DA) {
          plsOutMax=MUVR::plsda(subset(xIn,select=incVarMax),yIn,ncomp=nCompOutMax[i],near.zero.var=TRUE,scale=scale)
        }else{
          plsOutMax=MUVR::pls(subset(xIn,select=incVarMax),yIn,ncomp=nCompOutMax[i],near.zero.var=TRUE,scale=scale)
        }
        xTestMax=subset(xTest,select=incVarMax)
        yPredMaxR[testIndex]=predict(plsOutMax,newdata=xTestMax,scale=scale)$predict[,,nCompOutMax[i]]  #
        if (modReturn) {
          outMod[[i]]=list(plsOutMin,plsOutMid,plsOutMax)
        }
      } else {
        rfOutMin=randomForest(x=subset(xIn,select=incVarMin),y=yIn,xtest=subset(xTest,select=incVarMin),ytest=yTest,ntree=methParam$ntreeOut,keep.forest=TRUE)
        if (DA) {
          yPredMinR[testIndex,]=rfOutMin$test$votes
        } else {
          yPredMinR[testIndex]=rfOutMin$test$predicted
        }
        rfOutMid=randomForest(x=subset(xIn,select=incVarMid),y=yIn,xtest=subset(xTest,select=incVarMid),ytest=yTest,ntree=methParam$ntreeOut,keep.forest=TRUE)
        if (DA) {
          yPredMidR[testIndex,]=rfOutMid$test$votes
        } else {
          yPredMidR[testIndex]=rfOutMid$test$predicted
        }
        rfOutMax=randomForest(x=subset(xIn,select=incVarMax),y=yIn,xtest=subset(xTest,select=incVarMax),ytest=yTest,ntree=methParam$ntreeOut,keep.forest=TRUE)
        if (DA) {
          yPredMaxR[testIndex,]=rfOutMax$test$votes
        } else {
          yPredMaxR[testIndex]=rfOutMax$test$predicted
        }
        if (modReturn) {
          outMod[[i]]=list(rfOutMin,rfOutMid,rfOutMax)
        }
      }
    }##END OUTER LOOP

    #Per repetition: Average outer loop predictions, VIP ranks and nComp for PLS
    parReturn=list(yPredMin=yPredMinR,yPredMid=yPredMidR,yPredMax=yPredMaxR)
    # parReturn$VIPRepMin=apply(VIPOutMin,1,mean)
    # parReturn$VIPRepMid=apply(VIPOutMid,1,mean)
    # parReturn$VIPRepMax=apply(VIPOutMax,1,mean)
    # VI ranks
    parReturn$VIPRepMin <- rowMeans(VIPOutMin)
    parReturn$VIPRepMid <- rowMeans(VIPOutMid)
    parReturn$VIPRepMax <- rowMeans(VIPOutMax)

    if (method=='PLS'){
      parReturn$nCompRepMin=round(mean(nCompOutMin))
      parReturn$nCompRepMid=round(mean(nCompOutMid))
      parReturn$nCompRepMax=round(mean(nCompOutMax))
      parReturn$nCompSegMin=nCompOutMin
      parReturn$nCompSegMid=nCompOutMid
      parReturn$nCompSegMax=nCompOutMax
    }
    #Model validation curves
    parReturn$VAL=VALRep
    #Calculate nVar per repetition
    fitRankRep=colSums(VALRep)
    if(fitness=='AUROC') fitRankRep=-fitRankRep
    fitRankRep=(fitRankRep-min(fitRankRep))/abs(diff(range(fitRankRep)))
    if(all(is.nan(fitRankRep))) fitRankRep=rep(0,cnt) # If all VAL have same value -> reset all fitRankRep to 0
    minIndex=max(which(fitRankRep<=methParam$robust))
    maxIndex=min(which(fitRankRep<=methParam$robust))
    parReturn$varRepMin=var[minIndex]
    parReturn$varRepMid=round(exp(mean(log(c(var[minIndex],var[maxIndex])))))
    parReturn$varRepMax=var[maxIndex]
    #Return underlying models
    if (modReturn) parReturn$outModel=outMod
    reps[[r]]=parReturn

  }##END REPEATE DATA PARTITIONS

  #Collect all repetitions
  if (modReturn) outMods=list()
  for (r in 1:nRep) {
    if (DA) yPredMin[,,r]=reps[[r]]$yPredMin else
      yPredMin[,r]=reps[[r]]$yPredMin
    if (DA) yPredMid[,,r]=reps[[r]]$yPredMid else
      yPredMid[,r]=reps[[r]]$yPredMid
    if (DA) yPredMax[,,r]=reps[[r]]$yPredMax else
      yPredMax[,r]=reps[[r]]$yPredMax
    varRepMin[r]=reps[[r]]$varRepMin
    varRepMid[r]=reps[[r]]$varRepMid
    varRepMax[r]=reps[[r]]$varRepMax
    VIPRepMin[,r]=reps[[r]]$VIPRepMin
    VIPRepMid[,r]=reps[[r]]$VIPRepMid
    VIPRepMax[,r]=reps[[r]]$VIPRepMax
    if (method=='PLS') {
      nCompRepMin[r]=reps[[r]]$nCompRepMin
      nCompRepMid[r]=reps[[r]]$nCompRepMid
      nCompRepMax[r]=reps[[r]]$nCompRepMax
      nCompSegMin[r,]=reps[[r]]$nCompSegMin
      nCompSegMid[r,]=reps[[r]]$nCompSegMid
      nCompSegMax[r,]=reps[[r]]$nCompSegMax
    }
    VAL[,,r]=reps[[r]]$VAL
    if (modReturn) outMods=c(outMods,reps[[r]]$outModel)
  }

  #Average predictions
  if (DA) {
    yPred=list()
    yPred[['min']]=apply(yPredMin,c(1,2),mean)
    yPred[['mid']]=apply(yPredMid,c(1,2),mean)
    yPred[['max']]=apply(yPredMax,c(1,2),mean)
  } else {
    yPred=apply(yPredMin,1,mean)[1:nrow(X)]
    yPred=cbind(yPred,apply(yPredMid,1,mean)[1:nrow(X)])
    yPred=cbind(yPred,apply(yPredMax,1,mean)[1:nrow(X)])
    colnames(yPred)=c('min','mid','max')
    rownames(yPred)=paste(1:nSamp,ID,orgID,sep='_ID')
  }
  modelReturn$yPred=yPred
  modelReturn$yPredPerRep <- list(minModel = yPredMin, # And predictions per repetition
                                  midModel = yPredMid,
                                  maxModel = yPredMax)
  if (DA) {
    auc=matrix(nrow=3,ncol=length(levels(Y)),dimnames=list(c('min','mid','max'),levels(Y)))
    for (cl in 1:length(levels(Y))) {
      auc[1,cl]=roc(Y==(levels(Y)[cl]),yPred[['min']][,cl], quiet = TRUE)$auc
      auc[2,cl]=roc(Y==(levels(Y)[cl]),yPred[['mid']][,cl], quiet = TRUE)$auc
      auc[3,cl]=roc(Y==(levels(Y)[cl]),yPred[['max']][,cl], quiet = TRUE)$auc
    }
    #Classify predictions
    miss=numeric(3) #miss classification
    yClass=data.frame(Y) #class prediction
    for (mo in 1:3) {
      classPred=factor(apply(yPred[[mo]],1,function(x) levels(Y)[which.max(x)]),levels=levels(Y))
      miss[mo]=sum(classPred!=Y)
      yClass[,mo]=classPred
    }
    names(miss)=colnames(yClass)=c('min','mid','max')
    rownames(yClass)=paste(1:nSamp,ID,orgID,sep='_ID')
    #Report
    modelReturn$yClass=yClass
    modelReturn$miss=miss
    modelReturn$auc=auc
  }

  #Average VIP ranks over repetitions
  VIP=apply(VIPRepMin,1,mean)
  VIP=cbind(VIP,apply(VIPRepMid,1,mean))
  VIP=cbind(VIP,apply(VIPRepMax,1,mean))
  colnames(VIP)=c('min','mid','max')
  modelReturn$VIP=VIP
  modelReturn$VIPPerRep=list(minModel=VIPRepMin,midModel=VIPRepMid,maxModel=VIPRepMax)
  #Calculate overall nVar
  fitRankAll=apply(VAL,2,mean)
  if(fitness=='AUROC') fitRankAll=-fitRankAll
  fitRankAll=(fitRankAll-min(fitRankAll))/abs(diff(range(fitRankAll))) #rescale to 0-1 range
  if(all(is.nan(fitRankAll))) fitRankAll=rep(0,cnt) #If all VAL have same value -> reset all fitRankAll to 0
  minIndex=max(which(fitRankAll<=methParam$robust))
  maxIndex=min(which(fitRankAll<=methParam$robust))
  nVar=c(var[minIndex],round(exp(mean(log(c(var[minIndex],var[maxIndex]))))),var[maxIndex])
  names(nVar)=c('min','mid','max')
  modelReturn$nVar=nVar
  modelReturn$nVarPerRep <- list(minModel = varRepMin,
                                 midModel = varRepMid,
                                 maxModel = varRepMax)
  if (method=='PLS') {
    #Average nComp over repetitions
    nComp=c(round(mean(nCompRepMin)),round(mean(nCompRepMid)),round(mean(nCompRepMax)))
    names(nComp)=c('min','mid','max')
    modelReturn$nComp=nComp
    modelReturn$nCompPerRep <- list(minModel = nCompRepMin,
                                    midModel = nCompRepMid,
                                    maxModel = nCompRepMax)
    modelReturn$nCompPerSeg <- list(minModel = nCompSegMin,
                                    midModel = nCompSegMid,
                                    maxModel = nCompSegMax)
  }
  modelReturn$VAL$metric=fitness
  modelReturn$VAL$VAL=VAL
  if (modReturn) modelReturn$outModels=outMods
  modelReturn$inData=InData

  #Build overall "Fit" method for calculating R2 and visualisations
  incVarMin=names(VIP[rank(VIP[,1])<=round(nVar[1]),1])
  incVarMid=names(VIP[rank(VIP[,2])<=round(nVar[2]),2])
  incVarMax=names(VIP[rank(VIP[,3])<=round(nVar[3]),3])
  if (method=='PLS'){
    #Min model
    if (DA) plsFitMin=MUVR::plsda(subset(X,select=incVarMin),Y,ncomp=round(nComp[1]),near.zero.var=TRUE,scale=scale) else
      plsFitMin=MUVR::pls(subset(X,select=incVarMin),Y,ncomp=round(nComp[1]),near.zero.var=TRUE,scale=scale)
    if (length(plsFitMin$nzv$Position)>0) incVarMin=incVarMin[!incVarMin%in%rownames(plsFitMin$nzv$Metrics)]
    yFitMin=predict(plsFitMin,newdata=subset(X,select=incVarMin),scale=scale)$predict[,,nComp[1]]  #
    #Mid model
    if (DA) plsFitMid=MUVR::plsda(subset(X,select=incVarMid),Y,ncomp=round(nComp[2]),near.zero.var=TRUE,scale=scale) else
      plsFitMid=MUVR::pls(subset(X,select=incVarMid),Y,ncomp=round(nComp[2]),near.zero.var=TRUE,scale=scale)
    if (length(plsFitMid$nzv$Position)>0) incVarMid=incVarMid[!incVarMid%in%rownames(plsFitMid$nzv$Metrics)]
    yFitMid=predict(plsFitMid,newdata=subset(X,select=incVarMid),scale=scale)$predict[,,nComp[2]]  #
    #Max model
    if (DA) plsFitMax=MUVR::plsda(subset(X,select=incVarMax),Y,ncomp=round(nComp[3]),near.zero.var=TRUE,scale=scale) else
      plsFitMax=MUVR::pls(subset(X,select=incVarMax),Y,ncomp=round(nComp[3]),near.zero.var=TRUE,scale=scale)
    if (length(plsFitMax$nzv$Position)>0) incVarMax=incVarMax[!incVarMax%in%rownames(plsFitMax$nzv$Metrics)]
    yFitMax=predict(plsFitMax,newdata=subset(X,select=incVarMax),scale=scale)$predict[,,nComp[3]]  #
    yFit=cbind(yFitMin,yFitMid,yFitMax)
    yRep=ncol(yFit)/3
    colnames(yFit)=rep(c('min','mid','max'),each=yRep)
    rownames(yFit)=ID
    modelReturn$Fit=list(yFit=yFit,plsFitMin=plsFitMin,plsFitMid=plsFitMid,plsFitMax=plsFitMax)
  } else {
    rfFitMin=suppressWarnings(randomForest(subset(X,select=incVarMin),Y)) # suppress warnings for ML (regression against fewer than 5 unique values)
    if (DA) {
      yFitMin=rfFitMin$votes
    } else {
      yFitMin=rfFitMin$predicted
    }
    rfFitMid=suppressWarnings(randomForest(subset(X,select=incVarMid),Y)) # suppress warnings for ML (regression against fewer than 5 unique values)
    if (DA) {
      yFitMid=rfFitMid$votes
    } else {
      yFitMid=rfFitMid$predicted
    }
    rfFitMax=suppressWarnings(randomForest(subset(X,select=incVarMax),Y)) # suppress warnings for ML (regression against fewer than 5 unique values)
    if (DA) {
      yFitMax=rfFitMax$votes
    } else {
      yFitMax=rfFitMax$predicted
    }
    yFit=cbind(yFitMin,yFitMid,yFitMax)
    yRep=ncol(yFit)/3
    colnames(yFit)=rep(c('min','mid','max'),each=yRep)
    rownames(yFit)=ID
    modelReturn$Fit=list(yFit=yFit,rfFitMin=rfFitMin,rfFitMid=rfFitMid,rfFitMax=rfFitMax)
  }

  #Calculate fit statistics
  if (!DA) {
    TSS <- sum((Y - mean(Y)) ^ 2)
    RSSMin <- sum((Y - yFitMin) ^ 2)
    RSSMid <- sum((Y - yFitMid) ^ 2)
    RSSMax <- sum((Y - yFitMax) ^ 2)
    PRESSMin <- sum((Y - yPred[,1]) ^ 2)
    PRESSMid <- sum((Y - yPred[,2]) ^ 2)
    PRESSMax <- sum((Y - yPred[,3]) ^ 2)
    R2 <- c(1 - (RSSMin / TSS),
            1 - (RSSMid / TSS),
            1 - (RSSMax / TSS))
    Q2 <- c(1 - (PRESSMin / TSS),
            1 - (PRESSMid / TSS),
            1 - (PRESSMax / TSS))
    names(R2) <- names(Q2) <- c('min', 'mid', 'max')
    # Report
    modelReturn$fitMetric <- list(R2 = R2,
                                  Q2 = Q2)
  } else {
    modelReturn$fitMetric <- list(CR = 1 - (miss / length(Y)))
  }

  #Stop timer
  end.time=proc.time()[3]
  modelReturn$calcMins=(end.time-start.time)/60
  modelReturn$scale=scale; modelReturn$varRatio=varRatio; modelReturn$partitionSize=partitionSize;
  modelReturn$nRep=nRep; modelReturn$DA=DA; modelReturn$nOuter=nOuter; modelReturn$nInner=nInner;
  cat('\n Elapsed time',(end.time-start.time)/60,'mins \n')
  class(modelReturn)=c('MVObject',ifelse(DA,'Classification',ifelse(ML,'Multilevel','Regression')),method)
  return(modelReturn)
}
