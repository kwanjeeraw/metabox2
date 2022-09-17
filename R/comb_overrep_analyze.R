#'Perform integrated pathway overrepresentation analysis
#'@description perform integrated pathway overrepresentation analysis using hypergeometric tests.
#'If there are more than one entity types, Fisher's method is used to combine p-values before adjusted.
#'@usage comb_overrep_analyze(input_data, organism="hsa", size=3)
#'@param input_data a data frame of input data. The 1st column contains list of entities and the 2nd indicates entity types (i.e. gene, protein, compound). See details below.
#'@param organism a string specifying organism code from KEGG database. Choose one from hsa (default), tdc.
#'@param size a number specifying the minimum number of members in each annotation term to be used in the analysis. Default = 3.
#'@details Metabox uses KEGG ID (e.g.C12078) for compounds, UniProt entry (e.g.P0C9J6) for proteins, Ensembl (e.g.ENSG00000139618) for genes.
#'@return a list of the following components:
#'
#'enrichment = enrichment results.
#'
#'annotations = a data frame of annotation pairs.
#'
#'details = a list of analysis details: univsize, testMethod, pAdjusted, pCombined, inputsize, inputtype, numsets, minsize.
#'
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@references Johnson NL., Kotz S., and Kemp AW. (1992) Univariate Discrete Distributions, Second Edition. New York: Wiley.
#'@references Fisher R. (1932) Statistical methods for research workers. Oliver and Boyd, Edinburgh.
#'@seealso \code{\link{phyper}}, \code{\link{p.adjust}},\code{\link{pchisq}}
#'@examples
#'#out = comb_overrep_analyze(fnanal_data$combined_data, organism = "hsa", size = 3)
#'@import dplyr
#'@export
comb_overrep_analyze <- function (input_data, organism="hsa", size=3){
  #Check argument
  if (class(input_data) != "data.frame"){
    stop("argument 'input_data' is not valid, a data frame is required.")
  }
  tmparg <- try(organism <- match.arg(tolower(organism), c("hsa","tdc"), several.ok = FALSE), silent = TRUE)
  if (class(tmparg) == "try-error") {
    stop("argument 'organism' is not valid, choose one from the list: hsa,tdc")
  }
  #Initialize parameters
  overrep_result = list(); methodls = list(); methodls$univsize = 0; #Working data
  colnames(input_data) = c("id","type") #set colnames
  inptypels = unique(input_data$type)
  anno_dat = anno_databox[[organism]]  #load annotation data
  pair_dat = anno_dat[['pairData']] #load pair data
  meta_dat = anno_dat[['metaData']] #load meta data
  edgesls = do.call(rbind, lapply(inptypels, function(x){
    gsc = pair_dat[[paste0(tolower(x),'_annotation_pathway')]] #select data set
    comp_meta = dplyr::inner_join(x=input_data, y=gsc, by=c("id"="from"))[,c(1,3,2)]
    colnames(comp_meta) = c("from","to","type")
    comp_meta
  }))
  unqnodes = c(unique(edgesls$from),unique(edgesls$to))
  ntypels = c(unique(edgesls$type),"pathway")
  nodesls = do.call(rbind, lapply(ntypels, function(x){
    if(tolower(x) != "compound"){
      gsc = meta_dat[[tolower(x)]] #select data set
      dplyr::inner_join(x=data.frame(id=unqnodes), y=gsc, by=c("id"="GID"))
    }else{
      gsc = compound_databox #select data set
      comp_meta = dplyr::inner_join(x=data.frame(id=unqnodes), y=gsc, by=c("id"="kegg"))[,c(1,3,2)]
      colnames(comp_meta) = c("id","name","xref")
      comp_meta$type = "Compound"
      comp_meta
    }
  }))
  annonws = list(nodes=nodesls, edges=edgesls)
  annomem = plyr::ddply(annonws$edges,c("to"),plyr::summarise,member=list(from)) #get members
  getattb = dplyr::left_join(annonws$edges[,2:1],annonws$nodes,by=c("from"="id")) #get from attributes
  totalhit = getattb %>% dplyr::group_by(type) %>% dplyr::summarise(numuniq = n_distinct(from)) #get total no. of hit for each entity type
  subanno = getattb %>% dplyr::group_by(to,type) %>% dplyr::tally() #get no. of entities in each annotation term
  colnames(subanno) = c('id','nodelabel','count')
  subanno = subanno[subanno$count >= size,] #filter by annotation size
  if(nrow(subanno) > 0){
    cat("\nPerforming overrepresentation analysis ...\n")
    #get annotation details
    annostat = do.call(list, lapply(ntypels[ntypels!="pathway"], function (x){
      list(ntypesize=length(unique(pair_dat[[paste0(tolower(x),'_annotation_pathway')]]$from)), #get total no. of annotated entities (universe)
           annsize=pair_dat[[paste0(tolower(x),'_annotation_pathway')]] %>% dplyr::count(to) #get total size of each annotation term
      )
    }))
    names(annostat) = ntypels[ntypels!="pathway"]
    methodls$univsize = paste(names(annostat), unlist(annostat,recursive = F)[c(1,3)], collapse = ', ');
    overDF = data.frame(stringsAsFactors = FALSE)
    for(i in 1:nrow(subanno)){#overrepresentation analysis
      nodetyp = tolower(subanno$nodelabel[i])
      toid = subanno$id[i] #annotation term
      hitnum = subanno$count[i] #hit number
      annosize = annostat[[nodetyp]]$annsize #annotation size
      blackAnno = annostat[[nodetyp]]$ntypesize - annosize[annosize$to==toid,2] #no. of entities not in the annotation term
      kdrawn = totalhit[tolower(totalhit$type)==nodetyp,]$numuniq #nrow(nodelist)
      pval = phyper(hitnum-1, annosize[annosize$to==toid,2], blackAnno, kdrawn, lower.tail = F) #hypergeometric test
      hyp = data.frame(id=as.character(toid), p=pval, no_of_entities=hitnum,
                       annotation_size=annosize[annosize$to==toid,2], universe_size=annostat[[nodetyp]]$ntypesize, nodelabel=nodetyp, stringsAsFactors = FALSE)
      overDF = rbind(overDF,hyp)
    }
    overDF$p_adj = p.adjust(overDF$p, method = "BH") #adjust raw p-values
    #format output table
    pcombine = plyr::ddply(overDF,c('id'),plyr::summarise,p_combine=combine_pvals(p)) #Fisher combines raw p-values
    pcombine$p_combine_adj = p.adjust(pcombine$p_combine, method = "BH") #adjust raw p_combines
    pcombine$p = plyr::ddply(overDF,c('id'),plyr::summarise,p=list(paste0(nodelabel,' (',p,')')))$p
    pcombine$p_adj = plyr::ddply(overDF,c('id'),plyr::summarise,p_adj=list(paste0(nodelabel,' (',p_adj,')')))$p_adj
    pcombine$no_of_entities = plyr::ddply(overDF,c('id'),plyr::summarise,no_of_entities=list(paste0(nodelabel,' (',no_of_entities,')')))$no_of_entities
    pcombine$annotation_size = plyr::ddply(overDF,c('id'),plyr::summarise,annotation_size=list(paste0(nodelabel,' (',annotation_size,')')))$annotation_size
    pcombine$universe_size = plyr::ddply(overDF,c('id'),plyr::summarise,universe_size=list(paste0(nodelabel,' (',universe_size,')')))$universe_size
    pcombine = pcombine[order(pcombine$p_combine_adj),]
    pcombine$rank = seq(1:nrow(pcombine))
    pcombine = merge(annonws$nodes, pcombine, by='id') #merge annotation attributes and results
    pcombine = pcombine[,c(ncol(pcombine),1:(ncol(pcombine)-1))] #rearrange columns
    pcombine = dplyr::left_join(pcombine, annomem, by=c('id'='to'))
    overDF = pcombine #output
    memname = plyr::ddply(getattb,c('to'),plyr::summarise,membername=list(name))
    overDF = dplyr::left_join(overDF, memname, by=c('id'='to'))
    overDF$p = vapply(overDF$p, paste, collapse = ", ", character(1L)) #format list into a character vector
    overDF$p_adj = vapply(overDF$p_adj, paste, collapse = ", ", character(1L)) #format list into a character vector
    overDF$no_of_entities = vapply(overDF$no_of_entities, paste, collapse = ", ", character(1L)) #format list into a character vector
    overDF$annotation_size = vapply(overDF$annotation_size, paste, collapse = ", ", character(1L)) #format list into a character vector
    overDF$universe_size = vapply(overDF$universe_size, paste, collapse = ", ", character(1L)) #format list into a character vector
    overDF$member = vapply(overDF$member, paste, collapse = ", ", character(1L)) #format list into a character vector
    overDF$membername = vapply(overDF$membername, paste, collapse = " | ", character(1L)) #format list into a character vector
    cat("Returning ",nrow(overDF)," enrichment sets ...\n")
    ovr_res = list(enrichment=overDF, annotations=annonws$edges) #output
  }else{
    ovr_res = list(enrichment=data.frame(), annotations=data.frame()) #output
  }
  methodls$testMethod = paste("Integrated pathway overrepresentation analysis with Hypergeometric test");
  methodls$pAdjusted = "Benjamini & Hochberg or FDR"; methodls$pCombined = "Fisher's method";
  methodls$inputsize=length(input_data$id); methodls$inputtype=length(unique(input_data$type)); methodls$numsets=nrow(ovr_res$enrichment); methodls$minsize=size;
  overrep_result$enrichment = ovr_res$enrichment; overrep_result$annotations = ovr_res$annotations; overrep_result$details = methodls;
  return(overrep_result)
}
