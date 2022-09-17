#'Perform enrichment analysis
#'@description perform enrichment analysis using entity-level statistics e.g. p-values. The function wraps around the main functions of \pkg{\link{piano}}.
#'@usage enrichment_analyze(input_data, pcol=2, fccol=NULL, method="reporter",
#'nodetype="compound", settype="pathway", organism="hsa", size=3)
#'@param input_data a data frame of input data containing list of entities, statistics (e.g. p-values, t-values or F-values) and directions (e.g. fold-changes). See details below.
#'@param pcol a number specifying the column number that contains statistical values (e.g. p-values, t-values or F-values). Default is 2.
#'@param fccol a number specifying the column number that contains direction values (e.g. fold-changes). Default is NULL and directionality classes will be ignored.
#'@param method a string specifying the enrichment analysis method. It can be one of methods from \pkg{\link{piano}} e.g. reporter (default), fisher, median, mean, stouffer. See \code{\link{runGSA}}.
#'@param nodetype a string specifying a node type. It can be one of compound (default), protein, gene.
#'@param settype a string specifying a set type. It can be one of pathway (default), chemicalclass.
#'@param organism a string specifying organism code from KEGG database, the parameter will not be used for \code{settype = "chemicalclass"}. Choose one from hsa (default), tdc.
#'@param size a number specifying the minimum number of members in each annotation term to be used in the analysis. Default = 3.
#'@details For pathway analysis, Metabox uses KEGG ID (e.g.C12078) for compounds, UniProt entry (e.g.P0C9J6) for proteins, Ensembl (e.g.ENSG00000139618) for genes.
#'For chemical class analysis, HMDB ID (e.g.HMDB0000001) is used for compounds.
#'For \code{input_data}, 1st column = list of variable/feature IDs, 2nd column = statistics and 3rd column = directions (if applicable).
#'@return a list of the following components:
#'
#'enrichment = enrichment results.
#'
#'network = network format for \code{\link[piano:networkPlot]{piano::networkPlot()}}.
#'
#'details = a list of analysis details: univsize, testMethod, pAdjusted, inputsize, numsets, minsize.
#'
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@references Patil K. and Nielsen J. (2005) Uncovering transcriptional regulation of metabolism by using metabolic network topology. Proceedings of the National Academy of Sciences of the United States of America 102(8), 2685.
#'@references Oliveira A., Patil K., and Nielsen J. (2008) Architecture of transcriptional regulatory circuits is knitted over the topology of bio-molecular interaction networks. BMC Systems Biology 2(1), 17.
#'@references Väremo L., Nielsen J., and Nookaew I. (2013) Enriching the gene set analysis of genome-wide data by incorporating directionality of gene expression and combining statistical hypotheses and methods. Nucleic Acids Research, 41(8), pp. 4378-4391.
#'@seealso \code{\link[piano:loadGSC]{piano::loadGSC()}}, \code{\link[piano:runGSAhyper]{piano::runGSAhyper()}}, \code{\link[piano:GSAsummaryTable]{piano::GSAsummaryTable()}}
#'@examples
#'#out=enrichment_analyze(fnanal_data$compound_data, pcol=5, fccol=6)
#'#piano::networkPlot(out$network,class="non",significance=0.05)
#' @export
enrichment_analyze <- function(input_data, pcol=2, fccol=NULL, method="reporter", nodetype="compound", settype="pathway", organism="hsa", size=3){
  #Check argument
  tmparg <- try(nodetype <- match.arg(tolower(nodetype), c("compound","protein","gene"), several.ok = FALSE), silent = TRUE)
  if (class(tmparg) == "try-error") {
    stop("argument 'nodetype' is not valid, choose one from the list: compound,protein,gene")
  }
  tmparg <- try(settype <- match.arg(tolower(settype), c("pathway","chemicalclass"), several.ok = FALSE), silent = TRUE)
  if (class(tmparg) == "try-error") {
    stop("argument 'settype' is not valid, choose one from the list: pathway,chemicalclass")
  }
  tmparg <- try(method <- match.arg(tolower(method), c("reporter","fisher","median","mean","stouffer"), several.ok = FALSE), silent = TRUE)
  if (class(tmparg) == "try-error") {
    stop("argument 'method' is not valid, choose one from the list: reporter,fisher,median,mean,stouffer")
  }
  if (class(input_data) != "data.frame"){
    stop("argument 'input_data' is not valid, a data frame is required.")
  }
  tmparg <- try(organism <- match.arg(tolower(organism), c("hsa","tdc"), several.ok = FALSE), silent = TRUE)
  if (class(tmparg) == "try-error") {
    stop("argument 'organism' is not valid, choose one from the list: hsa,tdc")
  }
  #Initialize parameters
  enrich_result = list(); methodls = list(); methodls$univsize = 0; #Working data
  if(settype=="chemicalclass"){
    anno_dat = chemclass_databox  #load annotation data
    pair_dat = anno_dat[['pairData']] #load pair data
    meta_dat = anno_dat[['metaData']] #load meta data
    gsc = pair_dat[['compound_class']] #use only chem_class for this version
    nodeid = "hmdb"
  }else if(nodetype=="compound" && settype=="pathway"){
    anno_dat = anno_databox[[organism]]  #load annotation data
    pair_dat = anno_dat[['pairData']] #load pair data
    meta_dat = anno_dat[['metaData']] #load meta data
    gsc = pair_dat[[paste0(tolower(nodetype),'_annotation_',tolower(settype))]] #select data set
    nodeid = "kegg" #should support pubchem
  }else{#gene, protein pathways
    anno_dat = anno_databox[[organism]]  #load annotation data
    pair_dat = anno_dat[['pairData']] #load pair data
    meta_dat = anno_dat[['metaData']] #load meta data
    gsc = pair_dat[[paste0(tolower(nodetype),'_annotation_',tolower(settype))]] #select data set
    nodeid = "GID"
  }
  formatMembername <- function(rtable, nodetype){
    if(nodetype == "compound"){
      mg = merge(data.frame(GID=unlist(rtable$member)), compound_databox, by.x="GID", by.y=nodeid, sort=FALSE)
      mg$name
    }else{
      mg = merge(data.frame(GID=unlist(rtable$member)), meta_dat[[tolower(nodetype)]], sort=FALSE)
      mg$name
    }
  }

  cat("\nExecuting function ...\n")
  inst_pkg = NULL
  if(!requireNamespace("piano", quietly = TRUE)){#check and install required package
    cat("\nMissing the required package 'piano', trying to install the package ...\n")
    inst_pkg = install_pkgs('piano')
  }
  if(length(unlist(inst_pkg))){
    enr_res = list(enrichment = data.frame(), network = data.frame())
    cat("\nERROR! Could not install the required package 'piano'. Data was not analyzed.\n")
  }else{
    enr_res = tryCatch({
      gs = piano::loadGSC(gsc, type="data.frame")
      pval = input_data[,pcol]
      pval = as.numeric(pval)
      names(pval) = input_data[,1]
      fc = NULL
      if(!is.null(fccol) && fccol !=''){
        fc = input_data[,fccol]
        fc = as.numeric(fc)
        names(fc) = input_data[,1]
      }
      gsaRes = piano::runGSA(geneLevelStats=pval, directions=fc, gsc=gs, geneSetStat=method, gsSizeLim=c(size,Inf))
      methodls$univsize = length(unique(gsc$from))
      resTab = piano::GSAsummaryTable(gsaRes)
      if(nrow(resTab)>0){
        ## output
        cat("Formatting output ...\n")
        colnames(resTab) = gsub("Name","id",colnames(resTab))
        resTab = merge(resTab,meta_dat[[tolower(settype)]], by.x = 'id', by.y = 'GID')
        sumz = data.frame(table(gsc['to']))
        colnames(sumz) = c("to","freq")
        resTab = merge(resTab,sumz, by.x='id', by.y='to')
        resTab$member = lapply(resTab$id, function(x) names(piano::geneSetSummary(gsaRes, geneSet=x)$geneLevelStats)) #get members
        colnames(resTab) = gsub("Genes \\(tot\\)","hits",colnames(resTab))
        colnames(resTab) = gsub("Genes \\(up\\)","up",colnames(resTab))
        colnames(resTab) = gsub("Genes \\(down\\)","down",colnames(resTab))
        colnames(resTab) = gsub("freq","total",colnames(resTab))
        drops = c("id","name","total","hits","up","down")
        resTab = resTab[ ,c(drops, colnames(resTab)[!(colnames(resTab) %in% drops)])] #rearrange columns
        resTab = resTab[ ,c(1:(ncol(resTab)-3),ncol(resTab))]
        colnames(resTab) = gsub("\\.","_",colnames(resTab))
        colnames(resTab) = gsub("p adj","p_adj",colnames(resTab))
        resTab$membername = apply(resTab,1,formatMembername, nodetype=nodetype) #get member names
        resTab[is.na(resTab)] = '' #replace NA
        resTab$member = vapply(resTab$member, paste, collapse = ", ", character(1L)) #format list into a character vector
        resTab$membername = vapply(resTab$membername, paste, collapse = " | ", character(1L)) #format list into a character vector
        cat("Returning ",nrow(resTab)," enrichment sets ...\n")
        list(enrichment = resTab, network = gsaRes[-32])
      }else{
        list(enrichment = data.frame(), network = data.frame())
      }
    },
    error=function(e){
      cat(e$message)
      #message(e)
      cat("\nERROR! Data was not analyzed.\n")
      list(enrichment = data.frame(), network = data.frame())
    })
    methodls$testMethod = paste(settype, "enrichment analysis by", method); methodls$inputsize=nrow(input_data);
    methodls$numsets=nrow(enr_res$enrichment); methodls$minsize=size;  methodls$pAdjusted = "Benjamini & Hochberg or FDR";
  }
  enrich_result$enrichment = enr_res$enrichment; enrich_result$network = enr_res$network; enrich_result$details = methodls;
  return(enrich_result)
}
