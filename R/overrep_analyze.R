#'Perform overrepresentation analysis
#'@description perform overrepresentation analysis to retrieve overrepresented annotation terms of the input entities. The function wraps around the main functions of \pkg{\link{piano}}.
#'@usage overrep_analyze(txtinput, nodetype="compound", settype="pathway",
#'organism="hsa", size=3, universe)
#'@param txtinput a character vector of variable/feature IDs e.g. c('C12078', 'C02273'). See details below.
#'@param nodetype a string specifying a node type. It can be one of compound (default), protein, gene.
#'@param settype a string specifying a set type. It can be one of pathway (default), chemicalclass.
#'@param organism a string specifying organism code from KEGG database, the parameter will not be used for \code{settype = "chemicalclass"}. Choose one from hsa (default), tdc.
#'@param size a number specifying the minimum number of members in each annotation term to be used in the analysis. Default = 3.
#'@param universe a character vector of variable/feature IDs that represent the universe. Default's to all unique IDs in a set collection.
#'@details For pathway analysis, Metabox uses KEGG ID (e.g.C12078) for compounds, UniProt entry (e.g.P0C9J6) for proteins, NCBI-GeneID (e.g.10327) for genes.
#'For chemical class analysis, HMDB ID (e.g.HMDB0000001) is used for compounds.
#'@return a list of the following components:
#'
#'enrichment = enrichment results.
#'
#'network = network format for \code{\link[piano:networkPlot]{piano::networkPlot()}}.
#'
#'details = a list of analysis details: univsize, bgsize, testMethod, pAdjusted, inputsize, numsets, minsize.
#'
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@references Fisher R. (1932) Statistical methods for research workers. Oliver and Boyd, Edinburgh.
#'@references Väremo L., Nielsen J., and Nookaew I. (2013) Enriching the gene set analysis of genome-wide data by incorporating directionality of gene expression and combining statistical hypotheses and methods. Nucleic Acids Research, 41(8), pp. 4378-4391.
#'@seealso \pkg{\link{piano}}, \code{\link[piano:loadGSC]{piano::loadGSC()}}, \code{\link[piano:runGSAhyper]{piano::runGSAhyper()}}, \code{\link[piano:GSAsummaryTable]{piano::GSAsummaryTable()}}
#'@examples
#'#out=overrep_analyze(fnanal_data$compound_data$kegg, organism = "tdc", size = 5) #pathway ORA
#'#out=overrep_analyze(fnanal_data$combined_data$id[1:9], nodetype = "protein") #pathway ORA
#'#out=overrep_analyze(fnanal_data$compound_data$hmdb, settype="chemicalclass") #chemical class ORA
#'@export
overrep_analyze <- function (txtinput, nodetype="compound", settype="pathway", organism="hsa", size=3, universe){
  #Check argument
  tmparg <- try(nodetype <- match.arg(tolower(nodetype), c("compound","protein","gene"), several.ok = FALSE), silent = TRUE)
  if (class(tmparg) == "try-error") {
    stop("argument 'nodetype' is not valid, choose one from the list: compound,protein,gene")
  }
  tmparg <- try(settype <- match.arg(tolower(settype), c("pathway","chemicalclass"), several.ok = FALSE), silent = TRUE)
  if (class(tmparg) == "try-error") {
    stop("argument 'settype' is not valid, choose one from the list: pathway,chemicalclass")
  }
  tmparg <- try(organism <- match.arg(tolower(organism), c("hsa","tdc"), several.ok = FALSE), silent = TRUE)
  if (class(tmparg) == "try-error") {
    stop("argument 'organism' is not valid, choose one from the list: hsa,tdc")
  }
  #Initialize parameters
  overrep_result = list(); methodls = list(); methodls$univsize = 0; methodls$bgsize = 0; #Working data
  txtinput = unique(stringr::str_trim(unlist(txtinput))) #remove whiteline, duplicate
  txtinput = txtinput[!is.na(txtinput)]
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
    nodeid = "kegg"
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
    }else if(nodetype == "gene"){#use ncbi-geneid
      mg = merge(data.frame(GID=unlist(rtable$member)), meta_dat[[tolower(nodetype)]], by.x="GID", by.y="xref", sort=FALSE)
      mg$name
    }else{
      mg = merge(data.frame(GID=unlist(rtable$member)), meta_dat[[tolower(nodetype)]], sort=FALSE)
      mg$name
    }
  }

  cat("\nExecuting function ...\n")
  ovr_res = tryCatch({
    gs = piano::loadGSC(gsc, type="data.frame")
    if(missing(universe)){
      gsaRes = piano::runGSAhyper(genes = txtinput, gsc = gs, gsSizeLim=c(size,Inf))
    }else{
      gsaRes = piano::runGSAhyper(genes = txtinput, gsc = gs, gsSizeLim=c(size,Inf), universe=universe)
    }
    methodls$univsize = sum(gsaRes$contingencyTable[[1]])
    methodls$bgsize = sum(gsaRes$contingencyTable[[1]])-length(txtinput)
    resTab = data.frame(gsaRes$resTab)
    if(nrow(resTab)>0){
      ## output
      cat("Formatting output ...\n")
      colnames(resTab) = c('p_val','p_adj','hits','not_input','not_hits','not_in_set')
      resTab$id = row.names(resTab)
      resTab$set_size = resTab$hits + resTab$not_input
      resTab = merge(resTab,meta_dat[[tolower(settype)]], by.x = 'id', by.y = 'GID')
      cols = c("id","name","set_size","hits","not_hits","not_in_set","p_val","p_adj")
      resTab = resTab[ ,cols] #rearrange columns
      resTab$member = lapply(resTab$id, function(x) txtinput[txtinput %in% gsaRes$gsc[[x]]]) #get members
      resTab$membername = apply(resTab,1,formatMembername, nodetype=nodetype) #get member names
      resTab$member = vapply(resTab$member, paste, collapse = ", ", character(1L)) #format list into a character vector
      resTab$membername = vapply(resTab$membername, paste, collapse = " | ", character(1L)) #format list into a character vector
      cat("Returning ",nrow(resTab)," enrichment sets ...\n")
      list(enrichment = resTab, network = gsaRes)
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
  methodls$testMethod = paste(settype, "overrepresentation analysis with Fisher's exact test"); methodls$inputsize=length(txtinput);
  methodls$numsets=nrow(ovr_res$enrichment); methodls$minsize=size;  methodls$pAdjusted = "Benjamini & Hochberg or FDR";
  overrep_result$enrichment = ovr_res$enrichment; overrep_result$network = ovr_res$network; overrep_result$details = methodls;
  return(overrep_result)
}
