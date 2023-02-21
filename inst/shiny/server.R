server = function(input, output, session) {

  options(DT.options = list(pageLength = 20)) #Global option for dt-table output

  metabox_report_dir <- paste0(tempdir(),"/Metabox_report/")
  if(!dir.exists(metabox_report_dir)){
    dir.create(metabox_report_dir)
  }
  if(dir.exists(metabox_report_dir)){
    setwd(metabox_report_dir)
    tempdir <- getwd()
  }

  ##ClearValues###
  clear_all <- function(Inte = NULL, Inter = FALSE){

    if(!is.null(Inte)){
    if(Inte){
      if(!is.null(input$INfileMBPL1)){ ## Shiny retains input$INfile$datapath(file.csv) in tmpdir !!!!!, it must be deleted to avoid self-triggered reactive function in upload page
        if(file.exists(input$INfileMBPL1$datapath)){ file.remove(input$INfileMBPL1$datapath)}
      }
      if(!is.null(input$INfileMBPL2)){
        if(file.exists(input$INfileMBPL2$datapath)){ file.remove(input$INfileMBPL2$datapath)}
      }
      if(!is.null(input$INfileMBPL3)){
        if(file.exists(input$INfileMBPL3$datapath)){ file.remove(input$INfileMBPL3$datapath)}
      }
      if(!is.null(input$INfileMBPL4)){
        if(file.exists(input$INfileMBPL4$datapath)){ file.remove(input$INfileMBPL4$datapath)}
      }
      if(!is.null(input$INfileMBPL5)){
        if(file.exists(input$INfileMBPL5$datapath)){ file.remove(input$INfileMBPL5$datapath)}
      }
      if(!is.null(metboshow$metbo_read_mbpl5)){metboshow$metbo_read_mbpl5 <- NULL}
      if(!is.null(metboshow$metbo_read_mbpl4)){metboshow$metbo_read_mbpl4 <- NULL}
      if(!is.null(metboshow$metbo_read_mbpl3)){metboshow$metbo_read_mbpl3 <- NULL}
      if(!is.null(metboshow$metbo_read_mbpl2)){metboshow$metbo_read_mbpl2 <- NULL}
      if(!is.null(metboshow$metbo_read_mbpl1)){metboshow$metbo_read_mbpl1 <- NULL}

      if(!is.null(metboshow$metbo_data_mbpl1)){metboshow$metbo_data_mbpl1 <- NULL}
      if(!is.null(metboshow$metbo_data_mbpl2)){metboshow$metbo_data_mbpl2 <- NULL}
      if(!is.null(metboshow$metbo_data_mbpl3)){metboshow$metbo_data_mbpl3 <- NULL}
      if(!is.null(metboshow$metbo_data_mbpl4)){metboshow$metbo_data_mbpl4 <- NULL}
      if(!is.null(metboshow$metbo_data_mbpl5)){metboshow$metbo_data_mbpl5 <- NULL}

      if(!is.null(metboshow$metbo_data_mbpl1_m)){metboshow$metbo_data_mbpl1_m <- NULL}
      if(!is.null(metboshow$metbo_data_mbpl2_m)){metboshow$metbo_data_mbpl2_m <- NULL}
      if(!is.null(metboshow$metbo_data_mbpl3_m)){metboshow$metbo_data_mbpl3_m <- NULL}
      if(!is.null(metboshow$metbo_data_mbpl4_m)){metboshow$metbo_data_mbpl4_m <- NULL}
      if(!is.null(metboshow$metbo_data_mbpl5_m)){metboshow$metbo_data_mbpl5_m <- NULL}

      if(!is.null(metboshow$metbo_data_mbpl1_ana)){metboshow$metbo_data_mbpl1_ana <- NULL}
      if(!is.null(metboshow$metbo_data_mbpl2_ana)){metboshow$metbo_data_mbpl2_ana <- NULL}
      if(!is.null(metboshow$metbo_data_mbpl3_ana)){metboshow$metbo_data_mbpl3_ana <- NULL}
      if(!is.null(metboshow$metbo_data_mbpl4_ana)){metboshow$metbo_data_mbpl4_ana <- NULL}
      if(!is.null(metboshow$metbo_data_mbpl5_ana)){metboshow$metbo_data_mbpl5_ana <- NULL}
      if(!is.null(metboshow$metbo_mbpl)){metboshow$metbo_mbpl <- NULL}
      output$txtbox.missing_m <- renderPrint({cat("")})
      shinyjs::hide("next_m")
      shinyjs::hide("MIsummary_m")
      shinyjs::hide("MIoutput_m")
      hideTab(inputId = "MIoutputT_m", target = "MBPLdatS1_m")
      hideTab(inputId = "MIoutputT_m", target = "MBPLdatS2_m")
      hideTab(inputId = "MIoutputT_m", target = "MBPLdatS3_m")
      hideTab(inputId = "MIoutputT_m", target = "MBPLdatS4_m")
      hideTab(inputId = "MIoutputT_m", target = "MBPLdatS5_m")

      shinyjs::hide("idMBPL1")
      shinyjs::hide("idMBPL2")
      shinyjs::hide("idMBPL3")
      shinyjs::hide("idMBPL4")
      shinyjs::hide("idMBPL5")
      shinyjs::hide("classMBPL1")
      shinyjs::hide("firstMBPL1")
      shinyjs::hide("classMBPL2")
      shinyjs::hide("firstMBPL2")
      shinyjs::hide("classMBPL3")
      shinyjs::hide("firstMBPL3")
      shinyjs::hide("classMBPL4")
      shinyjs::hide("firstMBPL4")
      shinyjs::hide("classMBPL5")
      shinyjs::hide("firstMBPL5")
      shinyjs::hide("INfileMBPL1")
      shinyjs::hide("INfileMBPL2")
      shinyjs::hide("INfileMBPL3")
      shinyjs::hide("INfileMBPL4")
      shinyjs::hide("INfileMBPL5")
      shinyjs::hide("nameMBPL1")
      shinyjs::hide("nameMBPL2")
      shinyjs::hide("nameMBPL3")
      shinyjs::hide("nameMBPL4")
      shinyjs::hide("nameMBPL5")
      shinyjs::hide('MBPLsummary')
      shinyjs::hide("NextMBPL")
      shinyjs::disable("UploadMBPL")
      shinyjs::hide("next_m")
      shinyjs::hide("MBPLplotTab")

      shinyjs::reset("MBPLdmun")
      shinyjs::reset("idMBPL1")
      shinyjs::reset("idMBPL2")
      shinyjs::reset("idMBPL3")
      shinyjs::reset("idMBPL4")
      shinyjs::reset("idMBPL5")
      shinyjs::reset("classMBPL1")
      shinyjs::reset("firstMBPL1")
      shinyjs::reset("classMBPL2")
      shinyjs::reset("firstMBPL2")
      shinyjs::reset("classMBPL3")
      shinyjs::reset("firstMBPL3")
      shinyjs::reset("classMBPL4")
      shinyjs::reset("firstMBPL4")
      shinyjs::reset("classMBPL5")
      shinyjs::reset("firstMBPL5")
      shinyjs::reset("INfileMBPL1")
      shinyjs::reset("INfileMBPL2")
      shinyjs::reset("INfileMBPL3")
      shinyjs::reset("INfileMBPL4")
      shinyjs::reset("INfileMBPL5")
      shinyjs::reset("nameMBPL1")
      shinyjs::reset("nameMBPL2")
      shinyjs::reset("nameMBPL3")
      shinyjs::reset("nameMBPL4")
      shinyjs::reset("nameMBPL5")

      shinyjs::reset("missCheck_m")
      shinyjs::reset("reall_m")
      shinyjs::reset("remPercent_m")
      shinyjs::reset("MBPLnf")
      shinyjs::reset("MBPLoptdim")
      shinyjs::reset("MBPLnboot")
      shinyjs::reset("MBPLnrepet")
      shinyjs::reset("MBPLtest")
      shinyjs::reset("MBPLnpermut")
      shinyjs::reset("MBPLthreshold")
    }else{
      metboshow$S_missing_M1 <- FALSE
      metboshow$B_missing <- FALSE

      reset("INfile") #reset does not clear csv in tmp dir
      if(!is.null(input$INfile)){ ## Shiny retains input$INfile$datapath(file.csv) in tmpdir !!!!!, it must be deleted to avoid self-triggered reactive function in upload page
        if(file.exists(input$INfile$datapath)){ file.remove(input$INfile$datapath)}
      }
      if(!is.null(metboshow$metbo_read)){metboshow$metbo_read <- NULL}
      if(!is.null(metboshow$metbo_input)){metboshow$metbo_input <- NULL}
      if(!is.null(metboshow$metbo_QCnorm)){metboshow$metbo_QCnorm <- NULL}
      if(!is.null(metboshow$metbo_missing)){metboshow$metbo_missing <- NULL}
      if(!is.null(metboshow$metbo_norm)){metboshow$metbo_norm <- NULL}
      if(!is.null(metboshow$metbo_tran)){metboshow$metbo_tran <- NULL}
      if(!is.null(metboshow$metbo_scal)){metboshow$metbo_scal <- NULL}
      if(!is.null(metboshow$keepValueM)){metboshow$keepValueM <- NULL}
      if(!is.null(metboshow$keepValueN)){metboshow$keepValueN <- NULL}
      if(!is.null(metboshow$keepValueDP)){metboshow$keepValueDP <- NULL}
      if(!is.null(metboshow$keepValueMP)){metboshow$keepValueMP <- NULL}

      reset("idcol")
      reset("classS")
      reset("firstV")
      reset("next1")
      shinyjs::hide("UPdata")
      output$txtbox.missing <-  renderPrint({cat("")})
      reset("reall")
      reset("remPercent")
      reset("missCheck")
      reset("QCnormSelect")
      reset("exdesign")
      reset("istd")
      reset("sampleType")
      reset("injectionOrder")
      shinyjs::show("skip2")
      shinyjs::hide("next3")
      shinyjs::show("skip1")
      shinyjs::hide("next2")
      reset("DATAnormSelect")
      reset("ref")
      reset("tranM")
      reset("scalM")
      shinyjs::show("skip3")
      reset("var.equal")
      reset("ispara")
      reset("factor2Col")
      reset("doposthoc")
      if(!is.null(metboshow$metbo_uni)){metboshow$metbo_uni <- NULL}
      reset("multiM")
      reset("Mscale")
      reset("MULScorePC1")
      reset("MULScorePC2")
      reset("MULselectPC")
      if(!is.null(metboshow$metbo_multi)){metboshow$metbo_multi <- NULL}
      reset("corrM")
      if(!is.null(metboshow$metbo_corr)){metboshow$metbo_corr <- NULL}
      reset("LMEfix")
      reset("LMErandom")
      if(!is.null(metboshow$LMEfix)){metboshow$LMEfix <- NULL}
      reset("muvr")
      reset("nRep")
      reset("varRatio")
      reset("partitionSize")
      reset("da")
      reset("fitn")
      if(!is.null(metboshow$metbo_ml)){metboshow$metbo_ml <- NULL}
      reset("MBPLselectvipprop")
      reset("MBPLselectPC")
      reset("MBPLselectloadprop")
    }
    }

    if(is.null(Inte) && Inter == TRUE){
    if(!is.null(input$infileOver)){
      if(file.exists(input$infileOver$datapath)){ file.remove(input$infileOver$datapath)}
    }
    metboshow$EntityIDs <- NULL
    metboshow$metbo_overr <- NULL
    shinyjs::enable("EntityIDs")
    reset('EntityIDs')
    reset('INfileOVER')
    reset('OVERnode')
    reset('OVERset')
    reset('OVERsize')
    shinyjs::hide('OVERsummary')
    shinyjs::hide('OVERoutput')

    if(!is.null(input$metbo_ENRfile)){
      if(file.exists(input$metbo_ENRfile$datapath)){ file.remove(input$metbo_ENRfile$datapath)}
    }
    metboshow$metbo_ENRfile <- NULL
    metboshow$metbo_enr <- NULL
    reset('INfileENR')
    reset('ENRfccol')
    reset('ENRpcol')
    reset('ENRset')
    reset('ENRsize')
    reset('ENRmed')
    reset('ENRnode')
    shinyjs::hide('ENRsummary')
    shinyjs::hide('ENRoutput')

    if(!is.null(input$INfileORA)){
      if(file.exists(input$INfileORA$datapath)){ file.remove(input$INfileORA$datapath)}
    }
    metboshow$metbo_cORA_data <- NULL
    metboshow$metbo_cORA <- NULL
    reset('INfileORA')
    reset('ORAsize')
    shinyjs::hide('ORAsummary')
    shinyjs::hide('ORAoutput')
    }
  }

  metboshow <- reactiveValues()
  metboshow$S_missing_M1 <- FALSE
  metboshow$B_missing <- FALSE
  #Next button######
  observeEvent(input$start,{updateTabItems(session, "tabs", "upload")
    metboshow$full <- TRUE
    metboshow$stat <- FALSE
    metboshow$bio <- FALSE
    metboshow$nextMBPL_M <- FALSE
    clear_all(Inte = FALSE)
    shinyjs::js$resetSelected()
  })
  observeEvent(input$start2,{updateTabItems(session, "tabs", "upload")
    metboshow$full <- FALSE
    metboshow$stat <- TRUE
    metboshow$bio <- FALSE
    metboshow$nextMBPL_M <- FALSE
    clear_all(Inte = FALSE)
    shinyjs::js$resetSelected()

  })
  observeEvent(input$start3,{updateTabItems(session, "tabs", "upload")
    metboshow$full <- FALSE
    metboshow$stat <- FALSE
    metboshow$bio <- TRUE
    metboshow$nextMBPL_M <- FALSE
    clear_all(Inte = FALSE)
    shinyjs::js$resetSelected()

  })

  observeEvent(input$start4,{updateTabItems(session, "tabs", "DPtabs")
    updateTabItems(session, "tabs", "DIP")
    metboshow$full <- FALSE
    metboshow$stat <- FALSE
    metboshow$bio <- FALSE
    metboshow$nextMBPL_M <- FALSE
    clear_all(Inter = TRUE)
    shinyjs::js$resetSelected()
  })


  observeEvent(input$start5,{updateTabItems(session, "tabs", "mbpl")
    updateTabItems(session, "MBPLTabset", "mbplUP")
    metboshow$full <- FALSE
    metboshow$stat <- FALSE
    metboshow$bio <- FALSE
    metboshow$nextMBPL_M <- TRUE
    clear_all(Inte = TRUE)
    shinyjs::js$resetSelected()

  })
  observeEvent(input$next1,{
    if(metboshow$full == TRUE){
      updateTabItems(session, "tabs", "DPtabs")
      updateTabItems(session, "DPtabset", "missing")
    }
    if(metboshow$stat == TRUE){
      updateTabItems(session, "tabs", "statTabs") }
    if(metboshow$bio == TRUE){
      updateTabItems(session, "tabs", "bio") }
  })

  observeEvent(input$next2,{
    if(metboshow$full == TRUE){
      updateTabItems(session, "DPtabset", "QCnorm") }
    if(metboshow$stat == TRUE){
      updateTabItems(session, "tabs", "statTabs") }
    if(metboshow$bio == TRUE){
      updateTabItems(session, "tabs", "bio") }

  })

  observeEvent(input$next3,{updateTabItems(session, "DPtabset", "norm")})
  observeEvent(input$next4,{updateTabItems(session, "tabs", "statTabs")})
  observeEvent(input$next5,{updateTabItems(session, "tabs", "bio") })
  observeEvent(input$NextMBPL,{
    if(metboshow$nextMBPL_M == FALSE){
      updateTabItems(session, "MBPLTabset", "mbplAna")
    }else{
      shinyjs::hide("MIsummary_m")
      shinyjs::hide("MIoutput_m")
      shinyjs::hide("next_m")
      shinyjs::hide("txtbox.missing_m")
      updateTabItems(session, "MBPLTabset", "missing_m")
    }
  })
  observeEvent(input$next_m,{
    updateTabItems(session, "MBPLTabset", "mbplAna")
  })

  observeEvent(input$nextB_S,{updateTabItems(session, "tabs", "statTabs") })
  observeEvent(input$nextB_M,{updateTabItems(session, "tabs", "DPtabs")
    updateTabItems(session, "DPtabset", "missing")  })

  shinyjs::hide("nextB_M")
  observe({
    req(metboshow$bio); req(metboshow$B_missing);
    if(metboshow$bio == TRUE & metboshow$B_missing == TRUE){
      shinyjs::show("nextB_M")
    }else{
      shinyjs::hide("nextB_M")
    }
  })

  observeEvent(input$nextS_B1,{updateTabItems(session, "tabs", "bio") })
  observeEvent(input$nextS_B2,{updateTabItems(session, "tabs", "bio") })
  observeEvent(input$nextS_B3,{updateTabItems(session, "tabs", "bio") })
  observeEvent(input$nextS_B4,{updateTabItems(session, "tabs", "bio") })

  shinyjs::hide("nextS_M1")
  shinyjs::hide("nextS_M2")
  shinyjs::hide("nextS_M3")
  shinyjs::hide("nextS_M4")
  observe({
    req(metboshow$stat); req(metboshow$S_missing_M1);

    if(metboshow$stat == TRUE & metboshow$S_missing_M1 == TRUE){
      shinyjs::show("nextS_M1")
      shinyjs::show("nextS_M2")
      shinyjs::show("nextS_M3")
      shinyjs::show("nextS_M4")
    }
    else{
      shinyjs::hide("nextS_M1")
      shinyjs::hide("nextS_M2")
      shinyjs::hide("nextS_M3")
      shinyjs::hide("nextS_M4")
    }
  })

  observeEvent(input$nextS_M1,{updateTabItems(session, "tabs", "DPtabs")
    updateTabItems(session, "DPtabset", "missing") })
  observeEvent(input$nextS_M2,{updateTabItems(session, "tabs", "DPtabs")
    updateTabItems(session, "DPtabset", "missing") })
  observeEvent(input$nextS_M3,{updateTabItems(session, "tabs", "DPtabs")
    updateTabItems(session, "DPtabset", "missing") })
  observeEvent(input$nextS_M4,{updateTabItems(session, "tabs", "DPtabs")
    updateTabItems(session, "DPtabset", "missing") })

  observeEvent(input$skip1,{
    if(metboshow$full == TRUE){
      updateTabItems(session, "DPtabset", "QCnorm") }
    if(metboshow$stat == TRUE){
      updateTabItems(session, "tabs", "statTabs") }
    if(metboshow$bio == TRUE){
      updateTabItems(session, "tabs", "bio") }
  })

  observeEvent(input$skip2,{updateTabItems(session, "DPtabset", "norm") })

  #Back button######
  observeEvent(input$back1,{
    # if(metboshow$nextMBPL_M == TRUE){
    #   updateTabItems(session, "tabs", "home")
    # }else{
    updateTabItems(session, "tabs", "upload")
    # }
  })

  observe({
    if( metboshow$B_missing == TRUE || metboshow$S_missing_M1 == TRUE){
      shinyjs::hide("back1")
      shinyjs::hide("skip1")
    }
  })

  observeEvent(input$back2,{updateTabItems(session, "DPtabset", "missing")})
  observeEvent(input$back3,{updateTabItems(session, "DPtabset", "QCnorm") })

  observeEvent(input$back5,{
    if(metboshow$full == TRUE || metboshow$stat == TRUE){
      if(metboshow$full == TRUE){
        updateTabItems(session, "tabs", "DPtabs")
        updateTabItems(session, "DPtabset", "norm")
      }else{
        updateTabItems(session, "tabs", "statTabs")
      }
    }else{
      updateTabItems(session, "tabs", "upload")
    }
  })

  observeEvent(input$back6,{
    if(metboshow$full == TRUE || metboshow$bio == TRUE){
      if(metboshow$full == TRUE){
        updateTabItems(session, "tabs", "DPtabs")
        updateTabItems(session, "DPtabset", "norm")
      }else{
        updateTabItems(session, "tabs", "bio")
      }
    }else{
      updateTabItems(session, "tabs", "upload")
    }
  })
  observeEvent(input$back7,{
    if(metboshow$full == TRUE || metboshow$bio == TRUE){
      if(metboshow$full == TRUE){
        updateTabItems(session, "tabs", "DPtabs")
        updateTabItems(session, "DPtabset", "norm")
      }else{
        updateTabItems(session, "tabs", "bio")
      }
    }else{
      updateTabItems(session, "tabs", "upload")
    }
  })
  observeEvent(input$back8,{
    if(metboshow$full == TRUE || metboshow$bio == TRUE){
      if(metboshow$full == TRUE){
        updateTabItems(session, "tabs", "DPtabs")
        updateTabItems(session, "DPtabset", "norm")
      }else{
        updateTabItems(session, "tabs", "bio")
      }
    }else{
      updateTabItems(session, "tabs", "upload")
    }
  })
  observeEvent(input$back9,{
    if(metboshow$full == TRUE || metboshow$bio == TRUE){
      if(metboshow$full == TRUE){
        updateTabItems(session, "tabs", "DPtabs")
        updateTabItems(session, "DPtabset", "norm")
      }else{
        updateTabItems(session, "tabs", "bio")
      }
    }else{
      updateTabItems(session, "tabs", "upload")
    }
  })
  observeEvent(input$BackMBPL,{updateTabItems(session, "tabs", "home") })
  observeEvent(input$BackMBPL2,{
    if(metboshow$nextMBPL_M == TRUE){
      updateTabItems(session, "MBPLTabset", "missing_m")
    }else{
      updateTabItems(session, "MBPLTabset", "mbplUP") }
  })

  observeEvent(input$back_m,{updateTabItems(session, "MBPLTabset", "mbplUP")  })

  observeEvent(input$backOVER,{updateTabItems(session, "tabs", "home") })
  observeEvent(input$backENR,{updateTabItems(session, "tabs", "home") })
  observeEvent(input$backORA,{updateTabItems(session, "tabs", "home") })

  #####disable navbar####
  shinyjs::disable(selector = '.navbar-nav a[data-value="DPtabs"')
  shinyjs::disable(selector = '.navbar-nav a[data-value="statTabs"')
  shinyjs::disable(selector = '.navbar-nav a[data-value="bio"')
  shinyjs::js$disableTab("upload")
  shinyjs::js$disableTab("missing")
  shinyjs::js$disableTab("QCnorm")
  shinyjs::js$disableTab("norm")
  shinyjs::js$disableTab("DIP")
  shinyjs::js$disableTab("mbpl")
  shinyjs::js$disableTab("mbplUP")
  shinyjs::js$disableTab("missing_m")
  shinyjs::js$disableTab("mbplAna")


  ########################################################################

  ##first page, hide options before data uploading
  observe({
    shinyjs::hide("idcol")
    shinyjs::hide("classS")
    shinyjs::hide("firstV")
    shinyjs::hide("next1")
    req(metboshow$metbo_read)
    if (!is.null(metboshow$metbo_read)) {
      shinyjs::show("idcol")
      shinyjs::show("classS")
      shinyjs::show("firstV")
      shinyjs::show("next1")
    }
  })

  ##set input column before set metbo obj####
  observe({
    req(metboshow$metbo_read)
    if(ncol(metboshow$metbo_read) >= 3){
      updateSelectInput(session, "idcol",
                        choices=c(colnames(metboshow$metbo_read[1]),colnames(metboshow$metbo_read)))

      updateSelectInput(session, "classS",
                        choices=c(colnames(metboshow$metbo_read[2]),colnames(metboshow$metbo_read)))

      updateSelectInput(session, "firstV",
                        choices=c(colnames(metboshow$metbo_read[3]),colnames(metboshow$metbo_read)))
    }

  })

  ##option in impute missing value page#####
  observe({
    shinyjs::hide("remPercent")
    req(input$reall)
    if (input$reall == FALSE){
      observeEvent(input$reall,{shinyjs::show("remPercent")
        shinyjs::show("missCheck")} )}
    if (input$reall == TRUE){
      observeEvent(input$reall,{shinyjs::hide("remPercent")
        shinyjs::hide("missCheck")} )}
  })

  ##IS and QC page  page options######
  observe({
    req(metboshow$firstV)
    req(input$QCnormSelect)
    if (input$QCnormSelect != 'none'){
      if(input$QCnormSelect == 'ccmn' || input$QCnormSelect == 'ruv2'){shinyjs::show("exdesign")
        updateSelectInput(session, "exdesign",
                          choices=c("none",colnames(metboshow$keepValueN$inputdata[1:metboshow$firstV-1])))
      }else{shinyjs::hide("exdesign")}

      if(input$QCnormSelect == 'ccmn' || input$QCnormSelect == 'nomis' || input$QCnormSelect == 'ruvrand' || input$QCnormSelect == 'ruv2'){shinyjs::show("istd")
        updateSelectInput(session, "istd",
                          choices=c(NULL,colnames(metboshow$keepValueN$X)))
      }else{shinyjs::hide("istd")}

      if(input$QCnormSelect == 'serrf' || input$QCnormSelect == 'loess'){shinyjs::show("sampleType")
        updateSelectInput(session, "sampleType",
                          choices=c("none",colnames(metboshow$keepValueN$inputdata[1:metboshow$firstV-1])))
      }else{shinyjs::hide("sampleType")}

      if(input$QCnormSelect == 'serrf' || input$QCnormSelect == 'loess'){shinyjs::show("injectionOrder")
        updateSelectInput(session, "injectionOrder",
                          choices=c("none",colnames(metboshow$keepValueN$inputdata[1:metboshow$firstV-1])))
      }else{shinyjs::hide("injectionOrder")}

      if(input$QCnormSelect == 'serrf' || input$QCnormSelect == 'loess'){shinyjs::show("batch")
        updateSelectInput(session, "batch",
                          choices=c("none",colnames(metboshow$keepValueN$inputdata[1:metboshow$firstV-1])))
      }else{shinyjs::hide("batch")}
    }

    if (input$QCnormSelect == 'none'){
      shinyjs::hide("exdesign")
      shinyjs::hide("istd")
      shinyjs::hide("sampleType")
      shinyjs::hide("injectionOrder")
      shinyjs::hide("batch")
    }

    observe({
      req(input$DATAnormSelect)
      if(input$DATAnormSelect == "pqn"){shinyjs::show("ref")
        updateSelectInput(session, "ref",
                          choices=c(NULL,row.names(metboshow$keepValueDP$X),row.names(metboshow$keepValueDP$X)))
      }else{shinyjs::hide("ref")}
    })

    #option in data preprocessing######
    req(input$DATAnormSelect)
    if(input$DATAnormSelect == 'vsn' ){
      shinyjs::disable("scalM")
      shinyjs::disable("tranM")
      updateRadioButtons(session,"scalM" ,selected = "none")
      updateRadioButtons(session,"tranM" ,selected = "none")
    }else{
      shinyjs::enable("scalM")
      shinyjs::enable("tranM")}

    # if(input$norm == 'none' ){
    #   shinyjs::enable("scalM")
    #   shinyjs::enable("tranM")
    # }
  })

  #option in ML biomarker analysis
  observe({
    shinyjs::hide("fitn")
    req(input$da)
    if(input$da == 'tr'){shinyjs::show("fitn")}else{shinyjs::hide("fitn")}
  })

  #set value for all method#####
  shinyjs::hide("UPdata")
  metbo.input <- reactive({
    req(input$INfile)
    if(file.exists(input$INfile$datapath)){
      validate(need(tools::file_ext(input$INfile$datapath) == "csv", "File type is not supported, please use 'csv' file as input"))
      metboshow$metbo_read <- read_input_data(input$INfile$datapath)
      if(ncol(metboshow$metbo_read) < 3){
        cat("\nThe input requires at least 3 cloumns: ID column, class column and numeric value","\nPlease check the input format...")
        shinyjs::disable("next1")
        return(NULL)
      }else{
        shinyjs::enable("next1")
      }
      #cat("done reading")
      metboshow$idcol <- which(input$idcol==colnames(metboshow$metbo_read))
      metboshow$classS <- which(input$classS==colnames(metboshow$metbo_read))
      metboshow$firstV <- which(input$firstV==colnames(metboshow$metbo_read))

      req(metboshow$metbo_read);req(metboshow$idcol);req(metboshow$classS);req(metboshow$firstV);
      metboshow$metbo_input <- set_input_obj(metboshow$metbo_read,metboshow$idcol,metboshow$classS,
                                             metboshow$firstV)

      if(!is.null(metboshow$metbo_missing)){metboshow$metbo_missing <- NULL}
      if(!is.null(metboshow$metbo_norm)){metboshow$metbo_norm <- NULL}
      if(!is.null(metboshow$metbo_tran)){metboshow$metbo_tran <- NULL}
      if(!is.null(metboshow$metbo_scal)){metboshow$metbo_scal <- NULL}
      if(!is.null(metboshow$metbo_QCnorm)){metboshow$metbo_QCnorm <- NULL}
      #if(!is.null(metboshow$metbo_IStran)){metboshow$metbo_IStran <- NULL}
      if (!is.null(metboshow$keepValueMP)){metboshow$keepValueMP <- NULL}
      output$txtbox.missing <-  renderPrint({cat("")})
    }

  })

  output$txtbox <- renderPrint({
    if(is.null(metbo.input()) || is.null(metboshow$metbo_input)){ cat("")}
  })

  observe({
    #req(metboshow$metbo_input)
    req(metboshow$metbo_read)
    isolate({
      shinyjs::show("UPdata")
      output$INtable <- DT::renderDataTable({
        #metboshow$metbo_input$inputdata
        metboshow$metbo_read
      })
    })
  })

  #run impute missing value#####
  # hide("MIoutput")
  shinyjs::hide("next2")
  shinyjs::hide("MIsummary")
  shinyjs::hide("MIoutput")
  observeEvent(input$run1, {
    shinyjs::show("MIsummary")
    shinyjs::show("MIoutput")
    shinyjs::show("next2")
    shinyjs::show("INtable3")
    shinyjs::show("txtbox.missing")
    shinyjs::hide("skip1")
    shinyjs::hide("nextS_M1")
    shinyjs::hide("nextS_M2")
    shinyjs::hide("nextS_M3")
    shinyjs::hide("nextS_M4")
    shinyjs::hide("nextB_M")
    output$txtbox.QCnorm <- renderPrint({cat("")})
    if(!is.null(metboshow$metbo_norm)){metboshow$metbo_norm <- NULL}
    if(!is.null(metboshow$metbo_tran)){metboshow$metbo_tran <- NULL}
    if(!is.null(metboshow$metbo_scal)){metboshow$metbo_scal <- NULL}
    if(!is.null(metboshow$metbo_QCnorm)){metboshow$metbo_QCnorm <- NULL}
    #if(!is.null(metboshow$metbo_IStran)){metboshow$metbo_IStran <- NULL}
    if (!is.null(metboshow$keepValueMP)){metboshow$keepValueMP <- NULL}

    if(!is.null(metboshow$metbo_input)){
      if(sum(is.na(metboshow$metbo_input$X)) == 0){#No missing value to impute
        output$txtbox.missing <- renderPrint({cat("\nNo missing value to impute.\n")})
        shinyjs::hide("MIoutput")
        return(NULL)
      }

      shinybusy::show_modal_spinner(color = "#3CB371", text = "Computing ...")
      metboshow$metbo_missing <- impute_missing_data(metboshow$metbo_input, method = input$missCheck, removeall=input$reall, cutoff=input$remPercent )
      # shinyjs::hide("txtbox.QCnorm")
      if (input$reall == TRUE){
        output$txtbox.missing <- renderPrint({cat("\nRemove all variables with missing values.\n")})}else{
          if(input$missCheck=="min"){
            output$txtbox.missing <- renderPrint({
              cat("\nImpute missing values with a chosen method.\n")
              cat("\nRemove variables with >",isolate({input$remPercent}),"% missing values.\n")
              cat("\nImpute missing values with half of 'min' value of each column.\n")})
          }else{
            output$txtbox.missing <- renderPrint({
              cat("\nImpute missing values with a chosen method.\n")
              cat("\nRemove variables with >",isolate({input$remPercent}),"% missing values.\n")
              cat("\nImpute missing values with",isolate({input$missCheck}),"of each column.\n")})
          }
        }
      shinybusy::remove_modal_spinner()
    }else{
      output$txtbox.missing <- renderPrint({cat("\nNot found data to impute missing values.\n")})
      shinyjs::hide("MIoutput")
      return(NULL)
    }

  })

  observe({
    if(!is.null(metboshow$metbo_missing)){
      output$INtable3 <- DT::renderDataTable({ as.data.frame(metboshow$metbo_missing$X) })
    }else{
      output$txtbox.missing <-  renderPrint({cat("")})
      shinyjs::hide("MIsummary")
      shinyjs::hide("MIoutput")
      output$INtable3 <- NULL
    }
  })

  #set current value for is-qc normalization
  observe({
    if (is.null(metboshow$metbo_missing)){metboshow$keepValueN <- metboshow$metbo_input}
    if (!is.null(metboshow$metbo_missing)){metboshow$keepValueN <- metboshow$metbo_missing}
  })

  #is-qc normalization######
  shinyjs::hide("next3")
  observe({
    #add conditions for each QCnormSelect
    req(input$QCnormSelect);
    if(input$QCnormSelect == "none"){shinyjs::disable("run2")}else{shinyjs::enable("run2")}
    if(input$QCnormSelect == "ccmn"){
      if(input$exdesign != "none" && !is.null(input$istd)){shinyjs::enable("run2")}else{shinyjs::disable("run2")}
    }
    if(input$QCnormSelect == "nomis"){
      if(!is.null(input$istd)){shinyjs::enable("run2")}else{shinyjs::disable("run2")}
    }
    if(input$QCnormSelect == "ruvrand"){
      if(!is.null(input$istd)){shinyjs::enable("run2")}else{shinyjs::disable("run2")}
    }
    if(input$QCnormSelect == "ruv2"){
      if(input$exdesign != "none" && !is.null(input$istd)){shinyjs::enable("run2")}else{shinyjs::disable("run2")}
    }
    if(input$QCnormSelect == "serrf"){
      if(input$sampleType != "none" && input$injectionOrder != "none" && input$batch != "none"){shinyjs::enable("run2")}else{shinyjs::disable("run2")}
    }
    if(input$QCnormSelect == "loess"){
      if(input$sampleType != "none" && input$injectionOrder != "none" && input$batch != "none"){shinyjs::enable("run2")}else{shinyjs::disable("run2")}
    }
  })

  shinyjs::hide("QCsummary")
  shinyjs::hide("QCplotTab")

  observeEvent(input$run2, {
    shinyjs::show("next3")
    shinyjs::hide("skip2")
    shinyjs::show("QCsummary")
    shinyjs::show("QCplotTab")
    shinyjs::show("txtbox.QCnorm")
    shinyjs::show("INtable4")
    shinyjs::hide("Nsummary")
    shinyjs::hide("NplotTab")
    shinyjs::hide("txtbox.DATnorm")
    # shinyjs::show("txtbox.QCnorm")
    if(input$QCnormSelect == "none"){return(NULL)}
    if(input$QCnormSelect != "none"){

      if(!is.null(metboshow$metbo_QCnorm)){metboshow$metbo_QCnorm <- NULL}
      if(!is.null(metboshow$metbo_norm)){metboshow$metbo_norm <- NULL}
      if(!is.null(metboshow$metbo_tran)){metboshow$metbo_tran <- NULL}
      if(!is.null(metboshow$metbo_scal)){metboshow$metbo_scal <- NULL}
      if(!is.null(metboshow$keepValueMP)){metboshow$keepValueMP <- NULL}

      if(sum(is.na(metboshow$keepValueN$X)) > 0){#Data contains missing values
        output$txtbox.QCnorm <- renderPrint({
          shinyjs::hide("QCplotTab")
          cat("\nThe data contains missing values.\n")
          cat("\nThe data was not normalized.\n")
          cat("\nPlease impute missing values.\n")
        })
        shinyjs::hide("next3")
        shinyjs::show("skip2")
        return(NULL)
      }

      shinybusy::show_modal_spinner(color = "#3CB371", text = "Computing ...")
      metboshow$metbo_QCnorm <- normalize_input_data_byqc(metboshow$keepValueN,method = input$QCnormSelect, factorCol=which(input$exdesign == colnames(metboshow$keepValueN$inputdata)),
                                                          istd=match(input$istd,colnames(metboshow$keepValueN$X)),
                                                          sampleType=which(input$sampleType == colnames(metboshow$keepValueN$inputdata)),
                                                          injectionOrder=which(input$injectionOrder == colnames(metboshow$keepValueN$inputdata)),
                                                          batch=which(input$batch == colnames(metboshow$keepValueN$inputdata)))
      output$txtbox.QCnorm <- renderPrint({
        cat("\nExecuting function ...\n")
        cat("\nData normalization with\n",isolate({input$QCnormSelect}),"\nDone!\n")
      })
      shinybusy::remove_modal_spinner()

      # if(input$QCnormSelect == "none" & input$IStranM != "none"){
      #
      #   withProgress(message = 'Analysis in progress', min = 1,max = 100,value =10,
      #                {
      #                  metboshow$metbo_IStran <- transform_input_data(metboshow$keepValueN, method = input$IStranM)
      #
      #                  output$txtbox.QCnorm <- renderPrint({
      #                    cat("\nExecuting function ...\n")
      #                    cat("\nData transformation with\n",isolate({input$IStranM}),"\nDone!\n")
      #
      #                  })
      #
      #                })
      # }
      # if(input$QCnormSelect != "none" & input$IStranM == "none"){
      #   withProgress(message = 'Analysis in progress', min = 1,max = 100,value =10,
      #                {
      #                  metboshow$metbo_QCnorm <- normalize_input_data_byqc(metboshow$keepValueN,method = input$QCnormSelect, exdesign=which(input$exdesign == colnames(metboshow$keepValueN$inputdata)),
      #                                                                      istd=match(input$istd,colnames(metboshow$keepValueN$X)),
      #                                                                      sampleType=which(input$sampleType == colnames(metboshow$keepValueN$inputdata)),
      #                                                                      injectionOrder=which(input$injectionOrder == colnames(metboshow$keepValueN$inputdata)),
      #                                                                      batch=which(input$batch == colnames(metboshow$keepValueN$inputdata)))
      #                  output$txtbox.QCnorm <- renderPrint({
      #                    cat("\nExecuting function ...\n")
      #                    cat("\nData normalization with\n",isolate({input$QCnormSelect}),"\nDone!\n")
      #
      #                  })
      #
      #                })
      # }
      # if(input$QCnormSelect != "none" & input$IStranM != "none"){
      #   withProgress(message = 'Analysis in progress', min = 1,max = 100,value =10,
      #                {
      #                  metboshow$metbo_IStran <- transform_input_data(metboshow$keepValueN, method = input$IStranM)
      #
      #                  metboshow$metbo_QCnorm <- normalize_input_data_byqc(metboshow$metbo_IStran,method = input$QCnormSelect, exdesign=which(input$exdesign == colnames(metboshow$keepValueN$inputdata)),
      #                                                                      istd=match(input$istd,colnames(metboshow$keepValueN$X)),
      #                                                                      sampleType=which(input$sampleType == colnames(metboshow$keepValueN$inputdata)),
      #                                                                      injectionOrder=which(input$injectionOrder == colnames(metboshow$keepValueN$inputdata)),
      #                                                                      batch=which(input$batch == colnames(metboshow$keepValueN$inputdata)))
      #
      #                  output$txtbox.QCnorm <- renderPrint({
      #                    cat("\nExecuting function ...\n")
      #                    cat("\nData transformation with\n",isolate({input$IStranM}),"\nDone!\n")
      #                    cat("\nData normalization with\n",isolate({input$QCnormSelect}),"\nDone!\n")
      #
      #                  })
      #
      #                })
      # }

    }

  })

  shinyjs::hide("rev_is")
  observe({
    if(!is.null(metboshow$metbo_QCnorm)){
      shinyjs::show("deleteC_qc")
      output$INtable4 <- DT::renderDataTable({
        datatable(metboshow$metbo_QCnorm$X, selection = list(target = "column"))
      })

      output$plotQCoverOri <- renderPlotly({
        p_isO <- ggplotly(pcaplot_overview(metboshow$keepValueN,scale=FALSE, plot_title="Before normalization"), source = 'pca_isO')
        p_isO %>%
          layout(dragmode = "select") %>%
          event_register("plotly_selected")
      })

      output$txtbox.QCoverOri <- renderPrint({
        if(!is.null(event_data("plotly_selected", source = 'pca_isO'))){
          event_data("plotly_selected", source = 'pca_isO')
        }
      })

      #output$plotQCboxOri <- renderPlotly({ggplotly(boxplot_overview(metboshow$keepValueN$X[,1:ncol(metboshow$keepValueN$X)<=50]))})
      output$plotQCboxOriS <- renderPlot({rlaplot_overview(metboshow$keepValueN, limitx = TRUE, dolog = TRUE, plot_title="Before normalization")})
      output$plotQCdenOri <- renderPlot({densityplot_overview(metboshow$keepValueN$X, plot_title="Before normalization densityplot plot (variables)")})
      output$plotQCdenOriS <- renderPlot({densityplot_overview(metboshow$keepValueN$X, plotvar=FALSE, plot_title="Before normalization densityplot plot (samples)")})

      output$plotQCover <- renderPlotly({
        isolate({
          progress <- shiny::Progress$new()
          on.exit(progress$close())
          progress$set(message = "Drawing plot ", value = 0)
          for (i in 1:10) {
            progress$inc(0.1, detail = "...")
          }
          Sys.sleep(1)
        })
        p_is <- ggplotly(pcaplot_overview(metboshow$metbo_QCnorm,scale=FALSE, plot_title="After normalization"), source = 'pca_is')
        p_is %>%
          layout(dragmode = "select") %>%
          event_register("plotly_selected")
      })

      output$txtbox.QCover <- renderPrint({
        if(!is.null(event_data("plotly_selected", source = 'pca_is'))){
          shinyjs::show("rev_is")
          event_data("plotly_selected", source = 'pca_is')
        }else{shinyjs::hide("rev_is")}
      })

      # output$plotQCbox <- renderPlotly({ggplotly(boxplot_overview(metboshow$metbo_QCnorm$X[,1:ncol(metboshow$keepValueN$X)<=50]))})
      output$plotQCboxS <- renderPlot({
        isolate({
          progress <- shiny::Progress$new()
          on.exit(progress$close())
          progress$set(message = "Drawing plot ", value = 0)
          for (i in 1:10) {
            progress$inc(0.1, detail = "...")
          }
          Sys.sleep(1)
        })
        rlaplot_overview(metboshow$metbo_QCnorm, limitx = TRUE, dolog = TRUE, plot_title="After normalization")})
      output$plotQCden <- renderPlot({
        isolate({
          progress <- shiny::Progress$new()
          on.exit(progress$close())
          progress$set(message = "Drawing plot ", value = 0)
          for (i in 1:10) {
            progress$inc(0.1, detail = "...")
          }
          Sys.sleep(1)
        })
        densityplot_overview(metboshow$metbo_QCnorm$X, plot_title="After normalization densityplot plot (variables)")})
      output$plotQCdenS <- renderPlot({densityplot_overview(metboshow$metbo_QCnorm$X, plotvar=FALSE, plot_title="After normalization densityplot plot (samples)")})
    }else{
      shinyjs::hide("QCsummary")
      shinyjs::hide("QCplotTab")
      shinyjs::hide("txtbox.QCnorm")
      shinyjs::hide("INtable4")
      shinyjs::hide("deleteC_qc")
      output$INtable4 <- NULL
      output$plotQCoverOri <- renderPlotly({plotly_empty()})
      output$plotQCboxOri <- renderPlotly({plotly_empty()})
      output$plotQCdenOri <- renderPlotly({plotly_empty()})
      output$plotQCdenOriS <- renderPlotly({plotly_empty()})
      output$plotQCover <- renderPlotly({plotly_empty()})
      output$plotQCbox <- renderPlotly({plotly_empty()})
      output$plotQCden <- renderPlotly({plotly_empty()})
      output$plotQCdenS <- renderPlotly({plotly_empty()})
    }
  })

  observeEvent(input$deleteC_qc,{
    if (!is.null(input$INtable4_columns_selected)) {
      metboshow$metbo_QCnorm$X <- metboshow$metbo_QCnorm$X[,-as.numeric(input$INtable4_columns_selected)]
    }
  })

  observeEvent(input$rev_is,{

    if (!is.null(event_data("plotly_selected", source = 'pca_is'))){
      row_rem_is <- as.numeric(which(metboshow$metbo_QCnorm$orgID %in% event_data("plotly_selected", source = 'pca_is')$key))
      if (length(row_rem_is) != 0){

        metboshow$metbo_QCnorm$orgID <- metboshow$metbo_QCnorm$orgID[-row_rem_is]
        metboshow$metbo_QCnorm$ID <- metboshow$metbo_QCnorm$ID[-row_rem_is]
        metboshow$metbo_QCnorm$Y <- metboshow$metbo_QCnorm$Y[-row_rem_is]
        rm_unikID_IS <- which(metboshow$metbo_QCnorm$unikID %in% unique(row_rem_is))
        metboshow$metbo_QCnorm$unikID <- metboshow$metbo_QCnorm$unikID[-rm_unikID_IS]
        metboshow$metbo_QCnorm$unik <- metboshow$metbo_QCnorm$unik[-row_rem_is]
        metboshow$metbo_QCnorm$X <- metboshow$metbo_QCnorm$X[-row_rem_is,]
        metboshow$metbo_QCnorm$inputdata <- metboshow$metbo_QCnorm$inputdata[-row_rem_is,]

        output$txtbox.QCnorm <- renderPrint({
          cat("Removed sample ID:",event_data("plotly_selected", source = 'pca_is')$key, "\n")
        })
      }else{
        output$txtbox.QCnorm <-  renderPrint({
          cat("Please select value to remove", "\n")
        })
        return(NULL)}
    }else{
      output$txtbox.QCnorm <-  renderPrint({
        cat("Please select value to remove", "\n")
      })
      return(NULL)}
  })

  #set current value for data preprocessing
  observe({
    if (is.null(metboshow$metbo_missing)){metboshow$keepValueDP <- metboshow$metbo_input}
    if (!is.null(metboshow$metbo_missing)){metboshow$keepValueDP <- metboshow$metbo_missing}
    if (!is.null(metboshow$metbo_QCnorm)){metboshow$keepValueDP <- metboshow$metbo_QCnorm}
  })

  #data preprocessing#####
  observe({
    req(input$scalM)
    req(input$tranM)
    req(input$DATAnormSelect)
    if(input$scalM == "none" & input$tranM == "none" & input$DATAnormSelect == "none"){shinyjs::disable("run3")}else{shinyjs::enable("run3")}
  })

  shinyjs::hide("Nsummary")
  shinyjs::hide("NplotTab")

  observeEvent(input$run3, {
    shinyjs::show("Nsummary")
    shinyjs::show("NplotTab")
    shinyjs::show("txtbox.DATnorm")
    shinyjs::show("INtable5")
    if(!is.null(metboshow$metbo_norm)){metboshow$metbo_norm <- NULL}
    if(!is.null(metboshow$metbo_tran)){metboshow$metbo_tran <- NULL}
    if(!is.null(metboshow$metbo_scal)){metboshow$metbo_scal <- NULL}
    if(sum(is.na(metboshow$keepValueDP$X)) > 0){#Data contains missing values
      output$txtbox.DATnorm <- renderPrint({
        shinyjs::hide("Ntabset")
        cat("\nThe data contains missing values.\n")
        cat("\nThe data was not normalized.\n")
        cat("\nPlease impute missing values.\n")
      })
      return(NULL)
    }
    shinybusy::show_modal_spinner(color = "#3CB371", text = "Computing ...")

    ########1#SxTxNx
    if(input$scalM == "none" && input$tranM == "none" && input$DATAnormSelect == "none"){return(NULL)}
    ########2#SxTxNo
    if(input$scalM == "none" && input$tranM == "none" && input$DATAnormSelect != "none"){
      metboshow$metbo_norm <- normalize_input_data_bydata(metboshow$keepValueDP,method = input$DATAnormSelect,
                                                          ref = input$ref)

      output$txtbox.DATnorm <- renderPrint({
        cat("\nExecuting function ...\n")
        cat("\nData normalization with\n",isolate({input$DATAnormSelect}),"\nDone!\n")
      })

    }
    ########3#SxToNo
    if(input$scalM == "none" && input$tranM != "none" && input$DATAnormSelect != "none"){

      metboshow$metbo_norm <- normalize_input_data_bydata(metboshow$keepValueDP,method = input$DATAnormSelect,
                                                          ref = input$ref)
      metboshow$metbo_tran <- transform_input_data(metboshow$metbo_norm, method = input$tranM)

      output$txtbox.DATnorm <- renderPrint({
        cat("\nExecuting function ...\n")
        cat("\nData normalization with",isolate({input$DATAnormSelect}),"\n")
        cat("\nData transformation with",isolate({input$tranM}),"transformation.\n","\nDone!\n")
      })

    }

    ########4#SxToNx
    if(input$scalM == "none" && input$tranM != "none" && input$DATAnormSelect == "none"){
      metboshow$metbo_tran <- transform_input_data(metboshow$keepValueDP, method = input$tranM)
      output$txtbox.DATnorm <- renderPrint({
        cat("\nExecuting function ...\n")
        cat("\nData transformation with",isolate({input$tranM}),"transformation.\n","\nDone!\n")
      })
    }

    ########5#SoToNx
    if(input$scalM != "none" && input$tranM != "none" && input$DATAnormSelect == "none"){
      metboshow$metbo_tran <- transform_input_data(metboshow$keepValueDP, method = input$tranM)
      metboshow$metbo_scal <- scale_input_data(metboshow$metbo_tran, method = input$scalM)
      output$txtbox.DATnorm <- renderPrint({
        cat("\nExecuting function ...\n")
        cat("\nData transformation with",isolate({input$tranM}),"transformation.\n")
        cat("\nData scaling with",isolate({input$scalM}),"scaling.\n","\nDone!\n")
      })
    }
    ########6#SoTxNo
    if(input$scalM != "none" && input$tranM == "none" && input$DATAnormSelect != "none"){

      metboshow$metbo_norm <- normalize_input_data_bydata(metboshow$keepValueDP,method = input$DATAnormSelect,
                                                          ref = input$ref)
      metboshow$metbo_scal <- scale_input_data(metboshow$metbo_norm, method = input$scalM)
      output$txtbox.DATnorm <- renderPrint({
        cat("\nExecuting function ...\n")
        cat("\nData normalization with",isolate({input$DATAnormSelect}),"\n")
        cat("\nData scaling with",isolate({input$scalM}),"scaling.\n","\nDone!\n")
      })
    }
    ########7#SoTxNx
    if(input$scalM != "none" && input$tranM == "none" && input$DATAnormSelect == "none"){
      metboshow$metbo_scal <- scale_input_data(metboshow$keepValueDP, method = input$scalM)
      output$txtbox.DATnorm <- renderPrint({
        cat("\nExecuting function ...\n")
        cat("\nData scaling with",isolate({input$scalM}),"scaling.\n","\nDone!\n")
      })
    }
    ########8#SoToNo
    if(input$scalM != "none" && input$tranM != "none" && input$DATAnormSelect != "none"){

      metboshow$metbo_norm <- normalize_input_data_bydata(metboshow$keepValueDP,method = input$DATAnormSelect,
                                                          ref = input$ref)
      metboshow$metbo_tran <- transform_input_data(metboshow$metbo_norm, method = input$tranM)
      metboshow$metbo_scal <- scale_input_data(metboshow$metbo_tran, method = input$scalM)
      output$txtbox.DATnorm <- renderPrint({
        cat("\nExecuting function ...\n")
        cat("\nData normalization with",isolate({input$DATAnormSelect}),"\n")
        cat("\nData transformation with",isolate({input$tranM}),"transformation.\n")
        cat("\nData scaling with",isolate({input$scalM}),"scaling.\n","\nDone!\n")
      })
    }

    if (!is.null(metboshow$metbo_norm)){metboshow$keepValueMP <- metboshow$metbo_norm}
    if (!is.null(metboshow$metbo_tran)){metboshow$keepValueMP <- metboshow$metbo_tran}
    if (!is.null(metboshow$metbo_scal)){metboshow$keepValueMP <- metboshow$metbo_scal}

    shinybusy::remove_modal_spinner()

  })

  # output$txtbox.norm <- renderPrint({
  #   req(metbo.norm())
  #   if(is.null(metbo.norm())){ cat("")}
  # })

  #set current value for M-value for down stream analysis #####

  observe({
    if (is.null(metboshow$metbo_missing)){metboshow$keepValueM <- metboshow$metbo_input}
    if (!is.null(metboshow$metbo_missing)){metboshow$keepValueM <- metboshow$metbo_missing}
    if (!is.null(metboshow$metbo_QCnorm)){metboshow$keepValueM <- metboshow$metbo_QCnorm}
    if (!is.null(metboshow$keepValueMP)){metboshow$keepValueM <- metboshow$keepValueMP}

  })

  shinyjs::hide("rev_N")
  observe({
    if(!is.null(metboshow$keepValueMP)){
      shinyjs::show("deleteC_norm")
      output$INtable5 <- DT::renderDataTable({
        datatable(metboshow$keepValueMP$X,selection = list(target = "column"))
      })


      output$plotNoverOri <- renderPlotly({
        p_NO <-ggplotly(pcaplot_overview(metboshow$keepValueDP, scale=FALSE, plot_title="Before normalization"), source = 'pca_NO')
        p_NO %>%
          layout(dragmode = "select") %>%
          event_register("plotly_selected")
      })

      output$txtbox.NoverOri <- renderPrint({
        if(!is.null(event_data("plotly_selected", source = 'pca_NO'))){
          event_data("plotly_selected", source = 'pca_NO')
        }
      })

      if(!is.null(metboshow$metbo_scal)){
        output$plotNover <- renderPlotly({
          isolate({
            progress <- shiny::Progress$new()
            on.exit(progress$close())
            progress$set(message = "Drawing plot ", value = 0)
            for (i in 1:10) {
              progress$inc(0.1, detail = "...")
            }
            Sys.sleep(1)
          })
          p_N <- ggplotly(pcaplot_overview(metboshow$keepValueMP, scale=FALSE, plot_title="After normalization"), source = 'pca_N')
          p_N %>%
            layout(dragmode = "select") %>%
            event_register("plotly_selected")
        })
      }

      if(is.null(metboshow$metbo_scal)){
        output$plotNover <- renderPlotly({
          isolate({
            progress <- shiny::Progress$new()
            on.exit(progress$close())
            progress$set(message = "Drawing plot ", value = 0)
            for (i in 1:10) {
              progress$inc(0.1, detail = "...")
            }
            Sys.sleep(1)
          })
          p_N <- ggplotly(pcaplot_overview(metboshow$keepValueMP, scale = FALSE, plot_title="After normalization"), source = 'pca_N')
          p_N %>%
            layout(dragmode = "select") %>%
            event_register("plotly_selected")
        })
      }

      output$txtbox.Nover <- renderPrint({
        if(!is.null(event_data("plotly_selected", source = 'pca_N'))){
          shinyjs::show("rev_N")
          event_data("plotly_selected", source = 'pca_N')
        }else{shinyjs::hide("rev_N")}
      })

      #output$plotNboxOri <- renderPlotly({ggplotly(boxplot_overview(metboshow$keepValueDP$X[,1:ncol(metboshow$keepValueDP$X)<=50]))})
      output$plotNboxOriS <- renderPlot({rlaplot_overview(metboshow$keepValueDP, limitx = TRUE, dolog = TRUE, plot_title="Before normalization")})
      output$plotNdenOri <- renderPlot({densityplot_overview(metboshow$keepValueDP$X, plot_title="Before normalization densityplot plot (variables)")})
      output$plotNdenOriS <- renderPlot({densityplot_overview(metboshow$keepValueDP$X, plotvar=FALSE, plot_title="Before normalization densityplot plot (samples)")})
      # output$plotNbox <- renderPlotly({ggplotly(boxplot_overview(metboshow$keepValueMP$X[,1:ncol(metboshow$keepValueDP$X)<=50]))})
      output$plotNboxS <- renderPlot({
        isolate({
          progress <- shiny::Progress$new()
          on.exit(progress$close())
          progress$set(message = "Drawing plot ", value = 0)
          for (i in 1:10) {
            progress$inc(0.1, detail = "...")
          }
          Sys.sleep(1)
          if(!is.null(metboshow$metbo_tran)|!is.null(metboshow$metbo_scal)){
            rlaplot_overview(metboshow$keepValueMP, limitx = TRUE, dolog = FALSE, plot_title="After normalization")
          }else{
            rlaplot_overview(metboshow$keepValueMP, limitx = TRUE, dolog = TRUE, plot_title="After normalization")#dolog for normalization
          }
          })
      })

      output$plotNden <- renderPlot({
        isolate({
          progress <- shiny::Progress$new()
          on.exit(progress$close())
          progress$set(message = "Drawing plot ", value = 0)
          for (i in 1:10) {
            progress$inc(0.1, detail = "...")
          }
          Sys.sleep(1)
        })
        densityplot_overview(metboshow$keepValueMP$X, plot_title="After normalization densityplot plot (variables)")})
      output$plotNdenS <- renderPlot({densityplot_overview(metboshow$keepValueMP$X, plotvar=FALSE, plot_title="After normalization densityplot plot (samples)")})
    }else{
      shinyjs::hide("Nsummary")
      shinyjs::hide("NplotTab")
      shinyjs::hide("txtbox.DATnorm")
      shinyjs::hide("INtable5")
      shinyjs::hide("deleteC_norm")

      output$plotNoverOri <- renderPlotly({plotly_empty()})
      output$plotNboxOri <- renderPlotly({plotly_empty()})
      output$plotNdenOri <- renderPlotly({plotly_empty()})
      output$plotNdenOriS <- renderPlotly({plotly_empty()})
      output$plotNover <- renderPlotly({plotly_empty()})
      output$plotNbox <- renderPlotly({plotly_empty()})
      output$plotNden <- renderPlotly({plotly_empty()})
      output$plotNdenS <- renderPlotly({plotly_empty()})
    }

  })

  observeEvent(input$deleteC_norm,{
    if (!is.null(input$INtable5_columns_selected)) {
      metboshow$keepValueMP$X <- metboshow$keepValueMP$X[,-as.numeric(input$INtable5_columns_selected)]
    }
  })

  observeEvent(input$rev_N,{

    if (!is.null(event_data("plotly_selected", source = 'pca_N'))){
      row_rem_N <- as.numeric(which(metboshow$keepValueMP$orgID %in% event_data("plotly_selected", source = 'pca_N')$key))
      if (length(row_rem_N) != 0){

        metboshow$keepValueMP$ID <- metboshow$keepValueMP$ID[-row_rem_N]
        metboshow$keepValueMP$orgID <- metboshow$keepValueMP$orgID[-row_rem_N]
        metboshow$keepValueMP$Y <- metboshow$keepValueMP$Y[-row_rem_N]
        rm_unikID_N <- which(metboshow$keepValueMP$unikID %in% unique(row_rem_N))
        metboshow$keepValueMP$unikID <- metboshow$keepValueMP$unikID[-rm_unikID_N]
        metboshow$keepValueMP$unik <- metboshow$keepValueMP$unik[-row_rem_N]
        metboshow$keepValueMP$inputdata <- metboshow$keepValueMP$inputdata[-row_rem_N,]
        metboshow$keepValueMP$X <- metboshow$keepValueMP$X[-row_rem_N,]

        output$txtbox.DATnorm <- renderPrint({
          cat("Removed sample ID:",event_data("plotly_selected", source = 'pca_N')$key, "\n")
        })
      }else{
        output$txtbox.DATnorm <-  renderPrint({
          cat("Please select value to remove", "\n")
        })
        return(NULL)}
    }else{
      output$txtbox.DATnorm <-  renderPrint({
        cat("Please select value to remove", "\n")
      })
      return(NULL)}
  })

  #run statistical analysis####

  #run univariate analysis####
  observe({
    req(metboshow$firstV)
    updateSelectInput(session, "factor2Col",
                      choices=c("None" = "FALSE",colnames(metboshow$keepValueM$inputdata[which(colnames(metboshow$keepValueM$inputdata[1:metboshow$firstV-1]) != input$idcol &
                                                                                                 colnames(metboshow$keepValueM$inputdata[1:metboshow$firstV-1]) != input$classS
                      )]
                      )))
  })

  observe({
    metboshow$var.equal <- as.logical(input$var.equal)
    metboshow$ispara <- as.logical(input$ispara)
    metboshow$doposthoc<- as.logical(input$doposthoc)
    if (input$factor2Col == "FALSE"){
      metboshow$factor2Col <- as.logical(input$factor2Col)
      enable("var.equal")
    }else{
      metboshow$factor2Col <- which(input$factor2Col==colnames(metboshow$keepValueM$inputdata))
      disable("var.equal")
      updateRadioButtons(session,"var.equal" ,selected = "FALSE")

    }
  })

  #shinyjs::hide("UNIsummary")
  shinyjs::hide("UNIplotTab")
  shinyjs::hide("UNIoutput")
  metbo.uni <- eventReactive(input$run6, {
    #shinyjs::show("UNIsummary")
    output$INtableUNI <- NULL
    if(!is.null(metboshow$metbo_uni)){metboshow$metbo_uni <- NULL}
    if(sum(is.na(metboshow$keepValueM$X)) > 0){#Data contains missing values
      cat("\nThe data contains missing values.\n")
      cat("\nThe analysis was not performed.\n")
      cat("\nPlease impute missing values.\n")
      metboshow$metbo_uni <- NULL
      metboshow$S_missing_M1 <- TRUE
      shinyjs::hide("UNIplotTab")
      shinyjs::hide("UNIoutput")
      return(NULL)
    }
    #shinybusy::show_modal_spinner(color = "#3CB371", text = "Computing ...")

    if(!is.null(metboshow$keepValueM)){
      metboshow$metbo_uni = tryCatch({
        univ_analyze(metboshow$keepValueM, var.equal = metboshow$var.equal, ispara = metboshow$ispara, doposthoc = metboshow$doposthoc, factor2Col = metboshow$factor2Col)
      }
      ,error=function(e){
        shinybusy::remove_modal_spinner()
        message(e)

        #   # return(NULL)
        # showNotification("Error:",e$message,type = "error")
        #   # return(NULL)
      }
      )
    }else{return(NULL)}
    shinybusy::remove_modal_spinner()
    #return(metboshow$metbo_uni)
  })

  observe({
    req(metboshow$metbo_uni)
    isolate({
      if(!is.null(metboshow$metbo_uni)){

        shinyjs::show("UNIplotTab")
        shinyjs::show("UNIoutput")
        fcinx = grepl("_fc", colnames(metboshow$metbo_uni$stat_summ)) #get fold change column
        fcdf = as.data.frame(metboshow$metbo_uni$stat_sum[,fcinx],row.names = row.names(metboshow$metbo_uni$stat_sum)) #set fold change data frame
        colnames(fcdf) = colnames(metboshow$metbo_uni$stat_sum)[fcinx]
        if(input$factor2Col == "FALSE"){
          pval_dt = cbind(pvalue = metboshow$metbo_uni$p_value, adjusted_pvalue = metboshow$metbo_uni$p_adj) #set pval data frame
        }else{
          pval_dt = merge(t(metboshow$metbo_uni$p_value), t(metboshow$metbo_uni$p_adj), by=0,sort=F, suffixes=c(".pvalue",".adjusted_pvalue"))[,-1]
          row.names(pval_dt) = row.names(metboshow$metbo_uni$stat_summ)
        }
        metboshow$univ_merge_table = as.data.frame(merge(pval_dt,fcdf, by=0, sort=FALSE)) #set univariate output table
        colnames(metboshow$univ_merge_table)[1] = "variable"
        if(metboshow$metbo_uni$details$posthocTest == "No test"){
          output$INtableUNI <- DT::renderDataTable({metboshow$univ_merge_table})
        }else{
          poshoc_dt = get_posthoc_summary(metboshow$metbo_uni) #set post hoc table
          poshocdf = data.frame(variable=rownames(poshoc_dt),significant=poshoc_dt$significant)
          metboshow$poshoc_merge_table <-as.data.frame(merge(metboshow$univ_merge_table,poshocdf,sort=FALSE))
          metboshow$poshoc_merge_table$significant <- case_when(nchar(metboshow$poshoc_merge_table$significant) > 100 ~
                                                                  paste(str_sub(metboshow$poshoc_merge_table$significant, 1, 100), "..."),
                                                                nchar(metboshow$poshoc_merge_table$significant) <= 100 ~
                                                                  metboshow$poshoc_merge_table$significant)

          output$INtableUNI <- DT::renderDataTable({
            metboshow$poshoc_merge_table
          })
        }
      }else{
        shinyjs::hide("UNIplotTab")
        shinyjs::hide("UNIoutput")
        output$INtableUNI <- NULL
      }
    })

  })

  observe({
    if(!is.null(metboshow$metbo_uni)){
      # shinyjs::show("UNIsummary")
      # shinyjs::show("UNIplotTab")
      # shinyjs::show("UNIoutput")
      #metboshow$plotUNI <- data.frame(metboshow$metbo_uni$p_adj[-which(metboshow$metbo_uni$p_adj >= 0.05)])
      #output$INtablePoshocUNI <- DT::renderDataTable({get_posthoc_table(metboshow$metbo_uni)})

      output$plotPvalUNI <- renderPlotly({
        isolate({
          progress <- shiny::Progress$new()
          on.exit(progress$close())
          progress$set(message = "Drawing plot ", value = 0)
          for (i in 1:10) {
            progress$inc(0.1, detail = "...")
          }
          Sys.sleep(1)
          if(input$factor2Col == "FALSE"){
            filtpv = data.frame(sort(metboshow$metbo_uni$p_adj)[(1:length(metboshow$metbo_uni$p_adj))<=100])
            p_uni <- ggplotly(pvalplot_overview(filtpv,plot_title = "Statistical significance plot"))

          }else{
            f1inx = grep("adjusted_pvalue", colnames(metboshow$univ_merge_table))[1]
            pvtb = data.frame(metboshow$univ_merge_table[,f1inx], row.names = metboshow$univ_merge_table$variable)
            pvtb = pvtb[order(pvtb),, drop=FALSE]
            filtpv = head(pvtb,n=100)
            p_uni <- ggplotly(pvalplot_overview(filtpv,plot_title = "Statistical significance plot"))
          }
          p_uni %>%
            event_register("plotly_click")
        })
      })

      output$plotBoxUNI <- renderPlotly({
        req(event_data("plotly_click")$key)
        isolate({
          row_click <- as.numeric(which(colnames(metboshow$keepValueM$X) == event_data("plotly_click")$key))
          if(input$factor2Col != "FALSE"){
            ggplotly(boxplot_byF1F2(metboshow$keepValueM,row_click,metboshow$keepValueM$Y,metboshow$keepValueM$inputdata[,metboshow$factor2Col])) %>% layout(boxmode = "group")
          }else{
            ggplotly(boxplot_byF1(metboshow$keepValueM,row_click,metboshow$keepValueM$Y))
          }
        })
      })
    }else{
      shinyjs::hide("UNIplotTab")
      shinyjs::hide("UNIoutput")
      output$plotPvalUNI <- renderPlotly({plotly_empty()})
      output$plotBoxUNI <- renderPlotly({plotly_empty()})
    }
  })

  output$txtbox.uni <- renderPrint({
    if(is.null(metbo.uni())){
      cat("")
    }
    if(is.null(metboshow$metbo_uni)){
      cat("")
    }
  })


  #run multivariate analysis####
  shinyjs::hide("MULplotTab")

  metbo.multi <- eventReactive(input$run7, {
    #observeEvent(input$run7, {
    if(!is.null(metboshow$metbo_multi)){metboshow$metbo_multi <- NULL}
    shinyjs::show("MULplotTab")
    if(sum(is.na(metboshow$keepValueM$X)) > 0){#Data contains missing values
      cat("\nThe data contains missing values.\n")
      cat("\nThe analysis was not performed.\n")
      cat("\nPlease impute missing values.\n")
      metboshow$S_missing_M1 <- TRUE
      shinyjs::hide("MULplotTab")
      return(NULL)
    }
    shinybusy::show_modal_spinner(color = "#3CB371", text = "Computing ...")

    if(!is.null(metboshow$keepValueM)){
      metboshow$metbo_multi = tryCatch({
        multiv_analyze(metboshow$keepValueM,method=input$multiM, scale=input$Mscale)
      },error=function(e){
        shinybusy::remove_modal_spinner()
        message(e)
        return()
      })
    }else{return(NULL)}

    shinybusy::remove_modal_spinner()

    if(!is.null(metboshow$metbo_multi)){
      updateSelectInput(session, "MULselectPC",
                        choices=c(colnames(metboshow$metbo_multi$loading_val)[1],colnames(metboshow$metbo_multi$loading_val)))

      updateSelectInput(session, "MULScorePC1",
                        choices=c(colnames(metboshow$metbo_multi$score_val)[1],colnames(metboshow$metbo_multi$score_val)))

      updateSelectInput(session, "MULScorePC2",
                        choices=c(colnames(metboshow$metbo_multi$score_val)[2],colnames(metboshow$metbo_multi$score_val)))
    }else{return(NULL)}

  })

  observe({
    if(!is.null(metboshow$metbo_multi)){
      req(input$MULScorePC1);req(input$MULScorePC2);req(metboshow$metbo_multi)
      isolate({
        if((input$multiM == "pca" || input$multiM == "pls") && ncol(metboshow$metbo_multi$score_val) > 1 && input$MULScorePC2 != "NA"){#pca or pls
          shinyjs::show("MULScorePC1")
          shinyjs::show("MULScorePC2")
          metboshow$MULScoreCol1 <- which(input$MULScorePC1 == colnames(metboshow$metbo_multi$score_val))
          metboshow$MULScoreCol2 <- which(input$MULScorePC2 == colnames(metboshow$metbo_multi$score_val))
          output$plotPCloadMUL <- renderPlotly({#loading by pc
            progress <- shiny::Progress$new()
            on.exit(progress$close())
            progress$set(message = "Drawing plot ", value = 0)
            for (i in 1:10) {
              progress$inc(0.1, detail = "...")
            }
            Sys.sleep(1)
            filtld = metboshow$metbo_multi$loading_val[order(abs(metboshow$metbo_multi$loading_val[,input$MULselectPC]),decreasing = TRUE),]
            ggplotly(multiv_loadingplot_bypc(head(filtld,n=100), pc = input$MULselectPC))
          })
          output$plotScoreMUL <- renderPlotly({#score
            progress <- shiny::Progress$new()
            on.exit(progress$close())
            progress$set(message = "Drawing plot ", value = 0)
            for (i in 1:10) {
              progress$inc(0.1, detail = "...")
            }
            Sys.sleep(1)
            ggplotly(multiv_scoreplot(metboshow$keepValueM,cbind(PCx=metboshow$metbo_multi$score_val[,input$MULScorePC1],PCY=metboshow$metbo_multi$score_val[,input$MULScorePC2]),
                                      pcx=metboshow$metbo_multi$model_summary$R2X[metboshow$MULScoreCol1],pcy=metboshow$metbo_multi$model_summary$R2X[metboshow$MULScoreCol2], plot_title="Score plot"))
          })
          output$plotLoadMUL <- renderPlotly({#loading
            ggplotly(multiv_loadingplot(cbind(PCX=metboshow$metbo_multi$loading_val[,input$MULScorePC1],PCY=metboshow$metbo_multi$loading_val[,input$MULScorePC2]), plot_title="Loading plot"))
          })
          if(input$multiM == "pca"){#vip pca
            output$plotvipMUL <- renderPlotly({ggplotly(ggplot()+ggtitle("VIP wasn't calculated for PCA, showing no plot."))})
          }else{#vip pls
            output$plotvipMUL <- renderPlotly({
              progress <- shiny::Progress$new()
              on.exit(progress$close())
              progress$set(message = "Drawing plot ", value = 0)
              for (i in 1:10) {
                progress$inc(0.1, detail = "...")
              }
              Sys.sleep(1)
              filtvp = sort(metboshow$metbo_multi$vip_val, decreasing = TRUE)
              ggplotly(multiv_vipplot(head(filtvp,n=100)))
            })
          }
          metboshow$mul_table <- data.frame(cbind(metboshow$metbo_multi$loading_val,VIP=metboshow$metbo_multi$vip_val))
          output$INtableMUL <- DT::renderDataTable({
            metboshow$mul_table
          })
        }else if(input$multiM == "opls"){#opls
          shinyjs::hide("MULScorePC1")
          shinyjs::hide("MULScorePC2")
          if(!is.na(metboshow$metbo_multi$score_val)[,1] && !is.na(metboshow$metbo_multi$oscore_val)[,1]){
            output$plotPCloadMUL <- renderPlotly({#loading by pc
              progress <- shiny::Progress$new()
              on.exit(progress$close())
              progress$set(message = "Drawing plot ", value = 0)
              for (i in 1:10) {
                progress$inc(0.1, detail = "...")
              }
              Sys.sleep(1)
              filtld = metboshow$metbo_multi$loading_val[order(abs(metboshow$metbo_multi$loading_val[,1]),decreasing = TRUE),,drop=FALSE]
              ggplotly(multiv_loadingplot_bypc(head(filtld,n=100), pc = 1))
            })
            output$plotScoreMUL <- renderPlotly({#score
              progress <- shiny::Progress$new()
              on.exit(progress$close())
              progress$set(message = "Drawing plot ", value = 0)
              for (i in 1:10) {
                progress$inc(0.1, detail = "...")
              }
              Sys.sleep(1)
              ggplotly(multiv_scoreplot(metboshow$keepValueM,score_data=metboshow$metbo_multi$score_val,
                                        pcx=metboshow$metbo_multi$model_summary$R2X[1],pcy=metboshow$metbo_multi$model_summary$R2X[2],
                                        oscore_data=metboshow$metbo_multi$oscore_val, plot_title="Score plot"))
            })
            output$plotLoadMUL <- renderPlotly({#loading
              ggplotly(multiv_loadingplot(metboshow$metbo_multi$loading_val,oloading_data=metboshow$metbo_multi$oloading_val,plot_title="Loading plot"))
            })
            output$plotvipMUL <- renderPlotly({#vip
              progress <- shiny::Progress$new()
              on.exit(progress$close())
              progress$set(message = "Drawing plot ", value = 0)
              for (i in 1:10) {
                progress$inc(0.1, detail = "...")
              }
              Sys.sleep(1)
              filtvp = sort(metboshow$metbo_multi$vip_val, decreasing = TRUE)
              ggplotly(multiv_vipplot(head(filtvp,n=100)))
            })
            metboshow$mul_table <- data.frame(cbind(metboshow$metbo_multi$loading_val,metboshow$metbo_multi$oloading_val,VIP=metboshow$metbo_multi$vip_val,OVIP=metboshow$metbo_multi$ovip_val))
            output$INtableMUL <- DT::renderDataTable({
              metboshow$mul_table
            })
          }
        }else{#no model
          shinyjs::hide("MULScorePC1")
          shinyjs::hide("MULScorePC2")
          output$plotScoreMUL <- renderPlotly({ggplotly(ggplot()+ggtitle("Found only 1 PC, showing no plot."))})
          output$plotPCloadMUL <- renderPlotly({ggplotly(ggplot()+ggtitle("Found only 1 PC, showing no plot."))})
          output$plotvipMUL <- renderPlotly({ggplotly(ggplot()+ggtitle("Found only 1 PC, showing no plot."))})
          output$plotLoadMUL <- renderPlotly({plotly_empty()})
          output$INtableMUL <- DT::renderDataTable({data.frame()})
        }
      })
    }
    if(is.null(metboshow$metbo_multi)){
      shinyjs::hide("MULplotTab")
      shinyjs::hide("MULScorePC1")
      shinyjs::hide("MULScorePC2")
      output$plotScoreMUL <- renderPlotly({plotly_empty()})
      output$plotLoadMUL <- renderPlotly({plotly_empty()})
      output$plotvipMUL <- renderPlotly({plotly_empty()})
      output$plotPCloadMUL <- renderPlotly({plotly_empty()})
      output$INtableMUL <- DT::renderDataTable({data.frame()})
    }
  })

  output$txtbox.multi <- renderPrint({
    if(is.null(metbo.multi())){ cat("")}
    if(is.null(metboshow$metbo_multi)){
      cat("")
    }
  })

  #run correlation analysis####
  shinyjs::hide("CORplotTab")
  metbo.corr <- eventReactive(input$run8, {
    if(!is.null(metboshow$metbo_corr)){metboshow$metbo_corr <- NULL}
    shinyjs::show("CORplotTab")
    if(sum(is.na(metboshow$keepValueM$X)) > 0){#Data contains missing values
      cat("\nThe data contains missing values.\n")
      cat("\nThe analysis was not performed.\n")
      cat("\nPlease impute missing values.\n")
      metboshow$S_missing_M1 <- TRUE
      shinyjs::hide("CORplotTab")
      return(NULL)
    }

    # if(!is.null(input$INfileCorr)){
    #   if(tools::file_ext(input$INfileCorr$datapath) != "csv"){
    #     cat("File type is not supported, please use 'csv' file as input")
    #     return(metboshow$infileCorr <- NULL)
    #   }
    #   metboshow$infileCorr <- read.csv(input$INfileCorr$datapath, header = T)
    #   metboshow$infileCorr <- as.data.frame(metboshow$infileCorr)
    #   if (!is.numeric(metboshow$infileCorr)){
    #     cat("\nThe data contains non-numeric type values.\n")
    #     cat("\nThe analysis was not performed.\n\n")
    #     return(metboshow$infileCorr <- NULL)
    #   }
    # }else{metboshow$infileCorr <- NULL}
    shinybusy::show_modal_spinner(color = "#3CB371", text = "Computing ...")
    # metboshow$metbo_corr <- correlation_analyze(metboshow$keepValueM, y=metboshow$infileCorr, method=input$corrM)
    metboshow$metbo_corr <- correlation_analyze(metboshow$keepValueM, method=input$corrM)
    shinybusy::remove_modal_spinner()

  })

  output$txtbox.corr <- renderPrint({
    if(is.null(metbo.corr())){ cat("")}
    if(is.null(metboshow$metbo_corr)){
      cat("")
    }
  })

  observe({
    if(!is.null(metboshow$metbo_corr)){
      shinyjs::show("CORplotTab")
      isolate({
        output$INtableCOR <- DT::renderDataTable({data.frame(metboshow$metbo_corr$corr_data)})
        output$plotHeatCOR <- renderPlotly({ggplotly(corrplot_heatmap(metboshow$metbo_corr$corr_data, plot_title="Correlation heatmap"))})
      })
    }else{
      shinyjs::hide("CORplotTab")
      output$INtableCOR <- DT::renderDataTable({data.frame()})
      output$plotHeatCOR <-renderPlotly({plotly_empty()})
    }
  })

  #run Linear mixed-effects modeling####
  shinyjs::hide("LMEplotTab")
  #shinyjs::hide("LMEsummaryT")
  shinyjs::hide("LMEoutputT")

  observe({

    req(metboshow$firstV)
    updateSelectInput(session, "LMEfix",
                      choices=c(colnames(metboshow$keepValueM$inputdata[1]),colnames(metboshow$keepValueM$inputdata[1:metboshow$firstV-1])))
  })

  metbo.lme <- eventReactive(input$run9, {
    if(!is.null(metboshow$LMEfix)){metboshow$LMEfix <- NULL}

    metboshow$LMEfix <- match(input$LMEfix,colnames(metboshow$keepValueM$inputdata[1:metboshow$firstV-1]))

    if(sum(is.na(metboshow$keepValueM$X)) > 0){#Data contains missing values
      cat("\nThe data contains missing values.\n")
      cat("\nThe analysis was not performed.\n")
      cat("\nPlease impute missing values.\n")
      metboshow$S_missing_M1 <- TRUE
      metboshow$LMEfix <- NULL
      return(NULL)
    }
    shinybusy::show_modal_spinner(color = "#3CB371", text = "Computing ...")

    if(!is.null(metboshow$keepValueM)){
      metboshow$metbo_lme = tryCatch({
        lme_analyze(metboshow$keepValueM, fix = metboshow$LMEfix, random = input$LMErandom)
      },error=function(e){
        shinybusy::remove_modal_spinner()
        message(e)
        return()
      })
    }else{return(NULL)}
    shinybusy::remove_modal_spinner()
  })

  output$txtbox.lme <- renderPrint({
    if(is.null(metbo.lme())){ cat("")}
    if(is.null(metboshow$metbo_lme)){
      cat("")
    }
  })

  observe({
    if(!is.null(metboshow$metbo_lme)){
      output$INtableLME <- DT::renderDataTable({metboshow$metbo_lme$testTable})
      shinyjs::show("LMEplotTab")
      #shinyjs::show("LMEsummaryT")
      shinyjs::show("LMEoutputT")
      #output$plotPvalLME <- renderPlotly({ggplotly(pvalplot_overview(data.frame(metboshow$metbo_lme$table),plot_title = "Adjusted p-value plot"))})
    }else{
      output$INtableLME <- DT::renderDataTable({data.frame()})
      shinyjs::hide("LMEplotTab")
      shinyjs::hide("LMEoutputT")
    }
  })

  #run biomaker analysis####
  shinyjs::hide("MUVplotTab")

  metbo_ml <- eventReactive(input$run5, {
    if(!is.null(metboshow$metbo_ml)){metboshow$metbo_ml <- NULL}
    if(sum(is.na(metboshow$keepValueM$X)) > 0){#Data contains missing values
      shinyjs::hide("MUVplotTab")
      cat("\nThe data contains missing values.\n")
      cat("\nThe analysis was not performed.\n")
      cat("\nPlease impute missing values.\n")
      metboshow$B_missing <- TRUE
      return(NULL)
    }
    #metboshow$keepValueM$Y <- as.factor(metboshow$keepValueM$Y)
    shinybusy::show_modal_spinner(color = "#3CB371", text = "Computing ...")
    if(input$varRatio > 1 || input$partitionSize > 1 ){
      if(input$varRatio > 1){
        metboshow$varRatio <- 0.75
        cat("Var Ratio paramerter was set higher than 1.","\tThe default paramerter is used: 'Var Ratio = 0.75'")

      }
      if(input$partitionSize > 1){
        metboshow$partitionSize <- 0.7
        cat("Partition size paramerter was set higher than 1.","\tThe default paramerter is used: 'patition size = 0.7'")
      }
    }else{
      metboshow$varRatio <- as.numeric(input$varRatio)
      metboshow$partitionSize <- as.numeric(input$partitionSize)
    }
    if(!is.null(metboshow$keepValueM)){
      if(input$da == 'tr'){#classificcation
        DA = TRUE; fitness = input$fitn;
      }
      if(input$da == 'fa'){#regresion
        DA = FALSE; fitness = "RMSEP";
      }
      metboshow$metbo_ml = tryCatch({
        run_muvr(metboshow$keepValueM, method = input$muvr, nRep=input$nRep, varRatio=metboshow$varRatio, partitionSize=metboshow$partitionSize, DA=DA, fitness=fitness)
      },error=function(e){
        shinybusy::remove_modal_spinner()
        message(e)
        return()
      })
    }else{return(NULL)}
    shinybusy::remove_modal_spinner()
  })

  observe({
    if (!is.null(metboshow$metbo_ml)){
      shinyjs::show("MUVplotTab")
      output$plotVAL1 <- renderPlotly({ggplotly(MUVR_plotval(metboshow$metbo_ml,plot_title="Performnce evaluation plot"))})
      output$INtableVIP <- DT::renderDataTable({MUVR_getvip(metboshow$metbo_ml)})
      output$plotVIP <- renderPlotly({ggplotly(MUVR_plotvip(metboshow$metbo_ml,plot_title="Variable importance (VIP) plot of minimal-optimal model"))})
    }else{
      shinyjs::hide("MUVplotTab")
      output$plotVAL1 <- renderPlotly({plotly_empty()})
      output$INtableVIP <- DT::renderDataTable({data.frame()})
      output$plotVIP <- renderPlotly({plotly_empty()})
    }

  })

  output$txtbox.muv <- renderPrint({
    if(is.null(metbo_ml())){ cat("")}
    if(is.null(metboshow$metbo_ml)){cat("")}
  })

  ####upload files for mbpls-da ####
  shinyjs::hide('NextMBPL1')
  shinyjs::hide('NextMBPL2')
  shinyjs::hide('NextMBPL3')
  shinyjs::hide('NextMBPL4')
  shinyjs::hide('NextMBPL5')
  shinyjs::hide("idMBPL1")
  shinyjs::hide("idMBPL2")
  shinyjs::hide("idMBPL3")
  shinyjs::hide("idMBPL4")
  shinyjs::hide("idMBPL5")
  shinyjs::hide("classMBPL1")
  shinyjs::hide("firstMBPL1")
  shinyjs::hide("classMBPL2")
  shinyjs::hide("firstMBPL2")
  shinyjs::hide("classMBPL3")
  shinyjs::hide("firstMBPL3")
  shinyjs::hide("classMBPL4")
  shinyjs::hide("firstMBPL4")
  shinyjs::hide("classMBPL5")
  shinyjs::hide("firstMBPL5")
  shinyjs::hide("INfileMBPL1")
  shinyjs::hide("INfileMBPL2")
  shinyjs::hide("INfileMBPL3")
  shinyjs::hide("INfileMBPL4")
  shinyjs::hide("INfileMBPL5")
  shinyjs::hide("nameMBPL1")
  shinyjs::hide("nameMBPL2")
  shinyjs::hide("nameMBPL3")
  shinyjs::hide("nameMBPL4")
  shinyjs::hide("nameMBPL5")
  shinyjs::hide('MBPLsummary')
  shinyjs::hide("NextMBPL")
  shinyjs::disable("UploadMBPL")
  shinyjs::hide("next_m")
  shinyjs::hide("MBPLplotTab")
  hideTab(inputId = "MBPLtabDat", target = "MBPLdatS1")
  hideTab(inputId = "MBPLtabDat", target = "MBPLdatS2")
  hideTab(inputId = "MBPLtabDat", target = "MBPLdatS3")
  hideTab(inputId = "MBPLtabDat", target = "MBPLdatS4")
  hideTab(inputId = "MBPLtabDat", target = "MBPLdatS5")


  observeEvent(input$NextMBPL1,{
    shinyjs::show("INfileMBPL1")
    if(!is.null(metboshow$metbo_read_mbpl1)){
      shinyjs::show("nameMBPL1")
      shinyjs::show("idMBPL1")
      shinyjs::show("classMBPL1")
      shinyjs::show("firstMBPL1")
    }else{
      shinyjs::hide("nameMBPL1")
      shinyjs::hide("idMBPL1")
      shinyjs::hide("classMBPL1")
      shinyjs::hide("firstMBPL1")
    }
    # if(!is.null(metboshow$metbo_data_mbpl1)){
    #   if(sum(is.na(metboshow$metbo_data_mbpl1$X)) > 0){
    #     shinyjs::show("nextMBPL_M1")
    #   }else{
    #     shinyjs::hide("nextMBPL_M1")
    #   }
    # }else{
    #   shinyjs::hide("nextMBPL_M1")
    # }
    shinyjs::hide("idMBPL2")
    shinyjs::hide("idMBPL3")
    shinyjs::hide("idMBPL4")
    shinyjs::hide("idMBPL5")
    shinyjs::hide("nameMBPL2")
    shinyjs::hide("nameMBPL3")
    shinyjs::hide("nameMBPL4")
    shinyjs::hide("nameMBPL5")
    shinyjs::hide("classMBPL2")
    shinyjs::hide("firstMBPL2")
    shinyjs::hide("classMBPL3")
    shinyjs::hide("firstMBPL3")
    shinyjs::hide("classMBPL4")
    shinyjs::hide("firstMBPL4")
    shinyjs::hide("classMBPL5")
    shinyjs::hide("firstMBPL5")
    shinyjs::hide("INfileMBPL2")
    shinyjs::hide("INfileMBPL3")
    shinyjs::hide("INfileMBPL4")
    shinyjs::hide("INfileMBPL5")
  })
  observeEvent(input$NextMBPL2,{

    shinyjs::show("INfileMBPL2")
    if(!is.null(metboshow$metbo_read_mbpl2)){
      shinyjs::show("nameMBPL2")
      shinyjs::show("idMBPL2")
      shinyjs::show("classMBPL2")
      shinyjs::show("firstMBPL2")
    }else{
      shinyjs::hide("nameMBPL2")
      shinyjs::hide("idMBPL2")
      shinyjs::hide("classMBPL2")
      shinyjs::hide("firstMBPL2")
    }
    shinyjs::hide("idMBPL1")
    shinyjs::hide("idMBPL3")
    shinyjs::hide("idMBPL4")
    shinyjs::hide("idMBPL5")
    shinyjs::hide("nameMBPL1")
    shinyjs::hide("nameMBPL3")
    shinyjs::hide("nameMBPL4")
    shinyjs::hide("nameMBPL5")
    shinyjs::hide("classMBPL1")
    shinyjs::hide("firstMBPL1")
    shinyjs::hide("classMBPL3")
    shinyjs::hide("firstMBPL3")
    shinyjs::hide("classMBPL4")
    shinyjs::hide("firstMBPL4")
    shinyjs::hide("classMBPL5")
    shinyjs::hide("firstMBPL5")
    shinyjs::hide("INfileMBPL1")
    shinyjs::hide("INfileMBPL3")
    shinyjs::hide("INfileMBPL4")
    shinyjs::hide("INfileMBPL5")

  })
  observeEvent(input$NextMBPL3,{

    shinyjs::show("INfileMBPL3")
    if(!is.null(metboshow$metbo_read_mbpl3)){
      shinyjs::show("nameMBPL3")
      shinyjs::show("idMBPL3")
      shinyjs::show("classMBPL3")
      shinyjs::show("firstMBPL3")
    }else{
      shinyjs::hide("nameMBPL3")
      shinyjs::hide("idMBPL3")
      shinyjs::hide("classMBPL3")
      shinyjs::hide("firstMBPL3")
    }
    shinyjs::hide("idMBPL1")
    shinyjs::hide("idMBPL2")
    shinyjs::hide("idMBPL4")
    shinyjs::hide("idMBPL5")
    shinyjs::hide("nameMBPL1")
    shinyjs::hide("nameMBPL2")
    shinyjs::hide("nameMBPL4")
    shinyjs::hide("nameMBPL5")
    shinyjs::hide("classMBPL1")
    shinyjs::hide("firstMBPL1")
    shinyjs::hide("classMBPL2")
    shinyjs::hide("firstMBPL2")
    shinyjs::hide("classMBPL4")
    shinyjs::hide("firstMBPL4")
    shinyjs::hide("classMBPL5")
    shinyjs::hide("firstMBPL5")
    shinyjs::hide("INfileMBPL1")
    shinyjs::hide("INfileMBPL2")
    shinyjs::hide("INfileMBPL4")
    shinyjs::hide("INfileMBPL5")

  })
  observeEvent(input$NextMBPL4,{

    shinyjs::show("INfileMBPL4")
    if(!is.null(metboshow$metbo_read_mbpl4)){
      shinyjs::show("idMBPL4")
      shinyjs::show("nameMBPL4")
      shinyjs::show("classMBPL4")
      shinyjs::show("firstMBPL4")
    }else{
      shinyjs::hide("nameMBPL4")
      shinyjs::hide("idMBPL4")
      shinyjs::hide("classMBPL4")
      shinyjs::hide("firstMBPL4")
    }
    shinyjs::hide("idMBPL1")
    shinyjs::hide("idMBPL2")
    shinyjs::hide("idMBPL3")
    shinyjs::hide("idMBPL5")
    shinyjs::hide("nameMBPL1")
    shinyjs::hide("nameMBPL2")
    shinyjs::hide("nameMBPL3")
    shinyjs::hide("nameMBPL5")
    shinyjs::hide("classMBPL1")
    shinyjs::hide("firstMBPL1")
    shinyjs::hide("classMBPL3")
    shinyjs::hide("firstMBPL3")
    shinyjs::hide("classMBPL2")
    shinyjs::hide("firstMBPL2")
    shinyjs::hide("classMBPL5")
    shinyjs::hide("firstMBPL5")
    shinyjs::hide("INfileMBPL1")
    shinyjs::hide("INfileMBPL2")
    shinyjs::hide("INfileMBPL3")
    shinyjs::hide("INfileMBPL5")

  })
  observeEvent(input$NextMBPL5,{

    shinyjs::show("INfileMBPL5")
    if(!is.null(metboshow$metbo_read_mbpl5)){
      shinyjs::show("nameMBPL5")
      shinyjs::show("idMBPL5")
      shinyjs::show("classMBPL5")
      shinyjs::show("firstMBPL5")
    }else{
      shinyjs::hide("nameMBPL5")
      shinyjs::hide("idMBPL5")
      shinyjs::hide("classMBPL5")
      shinyjs::hide("firstMBPL5")
    }
    shinyjs::hide("idMBPL1")
    shinyjs::hide("idMBPL2")
    shinyjs::hide("idMBPL3")
    shinyjs::hide("idMBPL4")
    shinyjs::hide("nameMBPL1")
    shinyjs::hide("nameMBPL2")
    shinyjs::hide("nameMBPL3")
    shinyjs::hide("nameMBPL4")
    shinyjs::hide("classMBPL1")
    shinyjs::hide("firstMBPL1")
    shinyjs::hide("classMBPL3")
    shinyjs::hide("firstMBPL3")
    shinyjs::hide("classMBPL2")
    shinyjs::hide("firstMBPL2")
    shinyjs::hide("classMBPL4")
    shinyjs::hide("firstMBPL4")
    shinyjs::hide("INfileMBPL1")
    shinyjs::hide("INfileMBPL2")
    shinyjs::hide("INfileMBPL3")
    shinyjs::hide("INfileMBPL4")

  })

  ###read MBPL dataset 1###
  observe({
    req(input$INfileMBPL1)
    if(file.exists(input$INfileMBPL1$datapath)){
      validate(need(tools::file_ext(input$INfileMBPL1$datapath) == "csv", "File type is not supported, please use 'csv' file as input"))
      metboshow$metbo_read_mbpl1 <- tryCatch(
        {
          read.csv(file = input$INfileMBPL1$datapath, stringsAsFactors=TRUE, check.names=FALSE)
        }, error=function(e){
          print(e)
          return(NULL)
        })

      if(any(dim(metboshow$metbo_read_mbpl1) == 0) || ncol(metboshow$metbo_read_mbpl1) < 2) {
        shinyjs::show("MBPLsummary")
        output$txtbox.mbpl <- renderPrint({
          cat("\nERROR! Fail to read the input data set 1.\n"
              ,"Check the following guidelines:\n"
              ,"- the input data is a table of comma separated values with .csv file format.\n"
              ,"- samples are in rows and variables/features are in columns.\n"
              ,"- The input data requires at least 2 cloumns: class column and numeric value\n"
              ,"- missing values should be blank or NA.\n")
        })
        return(NULL)
      }
      return(metboshow$metbo_read_mbpl1)
    }else{return()}
  })
  observe({
    if(!is.null(metboshow$metbo_read_mbpl1)){
      shinyjs::show("nameMBPL1")
      shinyjs::show("idMBPL1")
      shinyjs::show("classMBPL1")
      shinyjs::show("firstMBPL1")
      if(ncol(metboshow$metbo_read_mbpl1) >= 2){
        updateTextInput(session, "nameMBPL1",
                        value= as.character(gsub(".csv","",input$INfileMBPL1$name)))

        updateSelectInput(session, "idMBPL1",
                          choices=c("none",colnames(metboshow$metbo_read_mbpl1)))

        updateSelectInput(session, "classMBPL1",
                          choices=c(colnames(metboshow$metbo_read_mbpl1[2]),colnames(metboshow$metbo_read_mbpl1)))

        updateSelectInput(session, "firstMBPL1",
                          choices=c(colnames(metboshow$metbo_read_mbpl1[3]),colnames(metboshow$metbo_read_mbpl1)))
      }
    }else{
      shinyjs::hide("nameMBPL1")
      shinyjs::hide("idMBPL1")
      shinyjs::hide("classMBPL1")
      shinyjs::hide("firstMBPL1")
      return(NULL)
    }
  })
  ###read MBPL dataset 2###
  observe({
    req(input$INfileMBPL2)
    if(file.exists(input$INfileMBPL2$datapath)){
      validate(need(tools::file_ext(input$INfileMBPL2$datapath) == "csv", "File type is not supported, please use 'csv' file as input"))
      metboshow$metbo_read_mbpl2 <- tryCatch(
        {
          read.csv(file = input$INfileMBPL2$datapath, stringsAsFactors=TRUE, check.names=FALSE)
        }, error=function(e){
          print(e)
          return(NULL)
        })

      if(any(dim(metboshow$metbo_read_mbpl2) == 0) || ncol(metboshow$metbo_read_mbpl2) < 2) {
        shinyjs::show("MBPLsummary")
        output$txtbox.mbpl <- renderPrint({
          cat("\nERROR! Fail to read the input data set 2.\n"
              ,"Check the following guidelines:\n"
              ,"- the input data is a table of comma separated values with .csv file format.\n"
              ,"- samples are in rows and variables/features are in columns.\n"
              ,"- The input data requires at least 2 cloumns: class column and numeric value\n"
              ,"- missing values should be blank or NA.\n")
        })
        return(NULL)
      }
      return(metboshow$metbo_read_mbpl2)
    }else(return())
  })
  observe({
    if(!is.null(metboshow$metbo_read_mbpl2)){
      shinyjs::show("nameMBPL2")
      shinyjs::show("idMBPL2")
      shinyjs::show("classMBPL2")
      shinyjs::show("firstMBPL2")
      if(ncol(metboshow$metbo_read_mbpl2) >= 2){
        updateTextInput(session, "nameMBPL2",
                        value= as.character(gsub(".csv","",input$INfileMBPL2$name)))

        updateSelectInput(session, "idMBPL2",
                          choices=c("none",colnames(metboshow$metbo_read_mbpl2)))

        updateSelectInput(session, "classMBPL2",
                          choices=c(colnames(metboshow$metbo_read_mbpl2[2]),colnames(metboshow$metbo_read_mbpl2)))

        updateSelectInput(session, "firstMBPL2",
                          choices=c(colnames(metboshow$metbo_read_mbpl2[3]),colnames(metboshow$metbo_read_mbpl2)))
      }
    }else{
      shinyjs::hide("nameMBPL2")
      shinyjs::hide("idMBPL2")
      shinyjs::hide("classMBPL2")
      shinyjs::hide("firstMBPL2")
      return(NULL)
    }
  })
  ###read MBPL dataset 3###
  observe({
    req(input$INfileMBPL3)
    if(file.exists(input$INfileMBPL3$datapath)){
      validate(need(tools::file_ext(input$INfileMBPL3$datapath) == "csv", "File type is not supported, please use 'csv' file as input"))
      metboshow$metbo_read_mbpl3 <- tryCatch(
        {
          read.csv(file = input$INfileMBPL3$datapath, stringsAsFactors=TRUE, check.names=FALSE)
        }, error=function(e){
          print(e)
          return(NULL)
        })

      if(any(dim(metboshow$metbo_read_mbpl3) == 0) || ncol(metboshow$metbo_read_mbpl3) < 2) {
        shinyjs::show("MBPLsummary")
        output$txtbox.mbpl <- renderPrint({
          cat("\nERROR! Fail to read the input data set 3.\n"
              ,"Check the following guidelines:\n"
              ,"- the input data is a table of comma separated values with .csv file format.\n"
              ,"- samples are in rows and variables/features are in columns.\n"
              ,"- The input data requires at least 2 cloumns: class column and numeric value\n"
              ,"- missing values should be blank or NA.\n")
        })
        return(NULL)
      }
      return(metboshow$metbo_read_mbpl3)
    }else(return())
  })

  observe({
    if(!is.null(metboshow$metbo_read_mbpl3)){
      shinyjs::show("nameMBPL3")
      shinyjs::show("idMBPL3")
      shinyjs::show("classMBPL3")
      shinyjs::show("firstMBPL3")
      if(ncol(metboshow$metbo_read_mbpl3) >= 2){
        updateTextInput(session, "nameMBPL3",
                        value= as.character(gsub(".csv","",input$INfileMBPL3$name)))

        updateSelectInput(session, "idMBPL3",
                          choices=c("none",colnames(metboshow$metbo_read_mbpl3)))

        updateSelectInput(session, "classMBPL3",
                          choices=c(colnames(metboshow$metbo_read_mbpl3[2]),colnames(metboshow$metbo_read_mbpl3)))

        updateSelectInput(session, "firstMBPL3",
                          choices=c(colnames(metboshow$metbo_read_mbpl3[3]),colnames(metboshow$metbo_read_mbpl3)))
      }
    }else{
      shinyjs::hide("nameMBPL3")
      shinyjs::hide("idMBPL3")
      shinyjs::hide("classMBPL3")
      shinyjs::hide("firstMBPL3")
      return(NULL)
    }
  })
  ###read MBPL dataset 4###
  observe({
    req(input$INfileMBPL4)
    if(file.exists(input$INfileMBPL4$datapath)){
      validate(need(tools::file_ext(input$INfileMBPL4$datapath) == "csv", "File type is not supported, please use 'csv' file as input"))
      metboshow$metbo_read_mbpl4 <- tryCatch(
        {
          read.csv(file = input$INfileMBPL4$datapath, stringsAsFactors=TRUE, check.names=FALSE)
        }, error=function(e){
          print(e)
          return(NULL)
        })

      if(any(dim(metboshow$metbo_read_mbpl4) == 0) || ncol(metboshow$metbo_read_mbpl4) < 2) {
        shinyjs::show("MBPLsummary")
        output$txtbox.mbpl <- renderPrint({
          shinyjs::show("MBPLsummary")
          cat("\nERROR! Fail to read the input data set 4.\n"
              ,"Check the following guidelines:\n"
              ,"- the input data is a table of comma separated values with .csv file format.\n"
              ,"- samples are in rows and variables/features are in columns.\n"
              ,"- The input data requires at least 2 cloumns: class column and numeric value\n"
              ,"- missing values should be blank or NA.\n")
        })
        return(NULL)
      }
      return(metboshow$metbo_read_mbpl4)
    }else{retun()}
  })
  observe({
    if(!is.null(metboshow$metbo_read_mbpl4)){
      shinyjs::show("nameMBPL4")
      shinyjs::show("idMBPL4")
      shinyjs::show("classMBPL4")
      shinyjs::show("firstMBPL4")
      if(ncol(metboshow$metbo_read_mbpl4) >= 2){
        updateTextInput(session, "nameMBPL4",
                        value= as.character(gsub(".csv","",input$INfileMBPL4$name)))

        updateSelectInput(session, "idMBPL4",
                          choices=c("none",colnames(metboshow$metbo_read_mbpl1)))

        updateSelectInput(session, "classMBPL4",
                          choices=c(colnames(metboshow$metbo_read_mbpl4[2]),colnames(metboshow$metbo_read_mbpl4)))

        updateSelectInput(session, "firstMBPL4",
                          choices=c(colnames(metboshow$metbo_read_mbpl4[3]),colnames(metboshow$metbo_read_mbpl4)))
      }
    }else{
      shinyjs::hide("nameMBPL4")
      shinyjs::hide("idMBPL4")
      shinyjs::hide("classMBPL4")
      shinyjs::hide("firstMBPL4")
      return(NULL)
    }
  })
  ###read MBPL dataset 5###
  observe({
    req(input$INfileMBPL5)
    if(file.exists(input$INfileMBPL5$datapath)){
      validate(need(tools::file_ext(input$INfileMBPL5$datapath) == "csv", "File type is not supported, please use 'csv' file as input"))
      metboshow$metbo_read_mbpl5 <- tryCatch(
        {
          read.csv(file = input$INfileMBPL5$datapath, stringsAsFactors=TRUE, check.names=FALSE)
        }, error=function(e){
          print(e)
          return(NULL)
        })

      if(any(dim(metboshow$metbo_read_mbpl5) == 0) || ncol(metboshow$metbo_read_mbpl5) < 2) {
        shinyjs::show("MBPLsummary")
        output$txtbox.mbpl <- renderPrint({
          cat("\nERROR! Fail to read the input data set 5.\n"
              ,"Check the following guidelines:\n"
              ,"- the input data is a table of comma separated values with .csv file format.\n"
              ,"- samples are in rows and variables/features are in columns.\n"
              ,"- The input data requires at least 2 cloumns: class column and numeric value\n"
              ,"- missing values should be blank or NA.\n")
        })
        return(NULL)
      }
      return(metboshow$metbo_read_mbpl5)
    }else{return()}
  })
  observe({
    if(!is.null(metboshow$metbo_read_mbpl5)){
      shinyjs::show("nameMBPL5")
      shinyjs::show("idMBPL5")
      shinyjs::show("classMBPL5")
      shinyjs::show("firstMBPL5")
      if(ncol(metboshow$metbo_read_mbpl5) >= 2){
        updateTextInput(session, "nameMBPL5",
                        value= as.character(gsub(".csv","",input$INfileMBPL5$name)))

        updateSelectInput(session, "idMBPL5",
                          choices=c("none",colnames(metboshow$metbo_read_mbpl5)))

        updateSelectInput(session, "classMBPL5",
                          choices=c(colnames(metboshow$metbo_read_mbpl5[2]),colnames(metboshow$metbo_read_mbpl5)))

        updateSelectInput(session, "firstMBPL5",
                          choices=c(colnames(metboshow$metbo_read_mbpl5[3]),colnames(metboshow$metbo_read_mbpl5)))
      }
    }else{
      shinyjs::hide("nameMBPL5")
      shinyjs::hide("idMBPL5")
      shinyjs::hide("classMBPL5")
      shinyjs::hide("firstMBPL5")
      return(NULL)
    }
  })

  ###MBPL etc###
  set_input_obj_decat <- function(inputdata, idCol=1, classCol=2, xCol=3){
    #Initialize parameters
    input_obj = list(); #Working data
    #Set ID; ID is numeric; for MUVR
    if (missing(idCol)) {
      ID = 1:nrow(inputdata) #list of samples
      orgID = ID #keep original IDs, for statistical analyses
    }else{
      ID = inputdata[,idCol] #list of ids
      orgID = ID #keep original IDs, for statistical analyses
      if (is.character(ID)) ID=factor(ID)
      if (!is.numeric(ID)) ID=as.numeric(ID) #convert IDs to numbers, for MUVR
    }
    # #Set classCol; classCol is numeric
    # if (missing(classCol)) {
    #   cat("\nMissing class/response column; Assuming it is next to sample ID column.\n")
    #   classCol = max(which(sapply(inputdata, is.factor)))
    # }
    # #Set xCol; xCol is numeric
    # if (missing(xCol)) {
    #   cat("\nMissing the first variable/feature column; Assuming it is next to class/response column.\n")
    #   xCol = 3
    # }

    ###INITIALIZE VARIABLES
    X = inputdata[,-1:-(xCol-1)] #X;
    Y = inputdata[,classCol] #Y;
    unik = !duplicated(ID)  #boolean of unique IDs; for MUVR
    unikID = ID[unik] #list of ids; for MUVR
    isRepeated = (!identical(unikID,ID)) #boolean of repeated experiment

    input_obj$inputdata = inputdata; input_obj$ID = ID;  input_obj$orgID = orgID; input_obj$unik = unik; input_obj$unikID = unikID; input_obj$idCol = idCol;
    input_obj$xCol = xCol; input_obj$classCol = classCol; input_obj$X = X; input_obj$Y = Y; input_obj$isRepeated = isRepeated; input_obj$details = list();
    class(input_obj)=c('METBObj')
    return(input_obj)
  }
  observe({
    if(input$MBPLdmun == 2){
      shinyjs::hide("NextMBPL")
      output$txtbox.mbpl <- renderPrint({
        cat("")
      })
      shinyjs::show("NextMBPL1")
      shinyjs::show("NextMBPL2")
      shinyjs::hide("NextMBPL3")
      shinyjs::hide("NextMBPL4")
      shinyjs::hide("NextMBPL5")
      shinyjs::hide("nameMBPL3")
      shinyjs::hide("nameMBPL4")
      shinyjs::hide("nameMBPL5")
      shinyjs::hide("idMBPL3")
      shinyjs::hide("idMBPL4")
      shinyjs::hide("idMBPL5")
      shinyjs::hide("classMBPL3")
      shinyjs::hide("firstMBPL3")
      shinyjs::hide("classMBPL4")
      shinyjs::hide("firstMBPL4")
      shinyjs::hide("classMBPL5")
      shinyjs::hide("firstMBPL5")
      shinyjs::hide("INfileMBPL3")
      shinyjs::hide("INfileMBPL4")
      shinyjs::hide("INfileMBPL5")
      shinyjs::hide("MBPLplotTab")
      shinyjs::disable("UploadMBPL")

      req(metboshow$metbo_read_mbpl1);req(metboshow$metbo_read_mbpl2);
      if(!is.null(metboshow$metbo_read_mbpl1) && !is.null(metboshow$metbo_read_mbpl2) ){
        shinyjs::enable("UploadMBPL")
      }else{
        shinyjs::disable("UploadMBPL")
      }
    }
    if(input$MBPLdmun == 3){
      shinyjs::hide("NextMBPL")
      output$txtbox.mbpl <- renderPrint({
        cat("")
      })
      shinyjs::show("NextMBPL1")
      shinyjs::show("NextMBPL2")
      shinyjs::show("NextMBPL3")
      shinyjs::hide("NextMBPL4")
      shinyjs::hide("NextMBPL5")
      shinyjs::hide("nameMBPL4")
      shinyjs::hide("nameMBPL5")
      shinyjs::hide("idMBPL4")
      shinyjs::hide("idMBPL5")
      shinyjs::hide("classMBPL4")
      shinyjs::hide("firstMBPL4")
      shinyjs::hide("classMBPL5")
      shinyjs::hide("firstMBPL5")
      shinyjs::hide("INfileMBPL4")
      shinyjs::hide("INfileMBPL5")
      shinyjs::hide("MBPLplotTab")
      shinyjs::disable("UploadMBPL")

      req(metboshow$metbo_read_mbpl1);req(metboshow$metbo_read_mbpl2);req(metboshow$metbo_read_mbpl3);
      if(!is.null(metboshow$metbo_read_mbpl1) && !is.null(metboshow$metbo_read_mbpl2) && !is.null(metboshow$metbo_read_mbpl3) ){
        shinyjs::enable("UploadMBPL")
      }else{
        shinyjs::disable("UploadMBPL")
      }
    }
    if(input$MBPLdmun == 4){
      shinyjs::hide("NextMBPL")
      output$txtbox.mbpl <- renderPrint({
        cat("")
      })
      shinyjs::show("NextMBPL1")
      shinyjs::show("NextMBPL2")
      shinyjs::show("NextMBPL3")
      shinyjs::show("NextMBPL4")
      shinyjs::hide("idMBPL5")
      shinyjs::hide("nameMBPL5")
      shinyjs::hide("NextMBPL5")
      shinyjs::hide("classMBPL5")
      shinyjs::hide("firstMBPL5")
      shinyjs::hide("INfileMBPL5")
      shinyjs::hide("MBPLplotTab")
      shinyjs::disable("UploadMBPL")

      req(metboshow$metbo_read_mbpl1);req(metboshow$metbo_read_mbpl2);req(metboshow$metbo_read_mbpl3);req(metboshow$metbo_read_mbpl4);
      if(!is.null(metboshow$metbo_read_mbpl1) && !is.null(metboshow$metbo_read_mbpl2) && !is.null(metboshow$metbo_read_mbpl3)
         && !is.null(metboshow$metbo_read_mbpl4)){
        shinyjs::enable("UploadMBPL")
      }else{
        shinyjs::disable("UploadMBPL")
      }
    }
    if(input$MBPLdmun == 5){
      shinyjs::hide("NextMBPL")
      output$txtbox.mbpl <- renderPrint({
        cat("")
      })
      shinyjs::show("NextMBPL1")
      shinyjs::show("NextMBPL2")
      shinyjs::show("NextMBPL3")
      shinyjs::show("NextMBPL4")
      shinyjs::show("NextMBPL5")
      shinyjs::hide("MBPLplotTab")
      shinyjs::disable("UploadMBPL")
      req(metboshow$metbo_read_mbpl1);req(metboshow$metbo_read_mbpl2);req(metboshow$metbo_read_mbpl3);req(metboshow$metbo_read_mbpl4);
      req(metboshow$metbo_read_mbpl5);
      if(!is.null(metboshow$metbo_read_mbpl1) && !is.null(metboshow$metbo_read_mbpl2) && !is.null(metboshow$metbo_read_mbpl3)
         && !is.null(metboshow$metbo_read_mbpl4) && !is.null(metboshow$metbo_read_mbpl5)){
        shinyjs::enable("UploadMBPL")
      }else{
        shinyjs::disable("UploadMBPL")
      }
    }

  })

  ###click o upload button befor MBPLSDA-analysis
  #more than 1000 wraning
  observeEvent(input$UploadMBPL,{
    if(!is.null(metboshow$metbo_data_mbpl1)){metboshow$metbo_data_mbpl1 <- NULL}
    if(!is.null(metboshow$metbo_data_mbpl2)){metboshow$metbo_data_mbpl2 <- NULL}
    if(!is.null(metboshow$metbo_data_mbpl3)){metboshow$metbo_data_mbpl3 <- NULL}
    if(!is.null(metboshow$metbo_data_mbpl4)){metboshow$metbo_data_mbpl4 <- NULL}
    if(!is.null(metboshow$metbo_data_mbpl5)){metboshow$metbo_data_mbpl5 <- NULL}
    if(!is.null(metboshow$metbo_data_mbpl1_m)){metboshow$metbo_data_mbpl1_m <- NULL}
    if(!is.null(metboshow$metbo_data_mbpl2_m)){metboshow$metbo_data_mbpl2_m <- NULL}
    if(!is.null(metboshow$metbo_data_mbpl3_m)){metboshow$metbo_data_mbpl3_m <- NULL}
    if(!is.null(metboshow$metbo_data_mbpl4_m)){metboshow$metbo_data_mbpl4_m <- NULL}
    if(!is.null(metboshow$metbo_data_mbpl5_m)){metboshow$metbo_data_mbpl5_m <- NULL}
    shinyjs::show("MBPLsummary")
    shinyjs::show("MBPLplotTab")
    ###set input object for MBPLS-DA analysis (2 data sets)
    if(input$MBPLdmun == 2){
      hideTab(inputId = "MBPLtabDat", target = "MBPLdatS3")
      hideTab(inputId = "MBPLtabDat", target = "MBPLdatS4")
      hideTab(inputId = "MBPLtabDat", target = "MBPLdatS5")
      if(!is.null(metboshow$metbo_read_mbpl1) && !is.null(metboshow$metbo_read_mbpl2) ){
        if(ncol(metboshow$metbo_read_mbpl1) >2000){
          output$txtbox.mbpl <- renderPrint({
            isolate({
              cat("")
            })
          })
        }
        metboshow$nameMBPL1 <- as.character(input$nameMBPL1)
        metboshow$classM1 <- which(input$classMBPL1==colnames(metboshow$metbo_read_mbpl1))
        metboshow$firstM1 <- which(input$firstMBPL1==colnames(metboshow$metbo_read_mbpl1))
        if(input$idMBPL1 != 'none'){
          metboshow$idM1 <- which(input$idMBPL1==colnames(metboshow$metbo_read_mbpl1))
          metboshow$metbo_data_mbpl1 <- set_input_obj_decat(metboshow$metbo_read_mbpl1,metboshow$idM1,metboshow$classM1,
                                                            metboshow$firstM1)
        }else{
          metboshow$metbo_data_mbpl1 <- set_input_obj_decat(metboshow$metbo_read_mbpl1,classCol=metboshow$classM1,
                                                            xCol=metboshow$firstM1)
        }
        metboshow$nameMBPL2 <- as.character(input$nameMBPL2)
        metboshow$classM2 <- which(input$classMBPL2==colnames(metboshow$metbo_read_mbpl2))
        metboshow$firstM2 <- which(input$firstMBPL2==colnames(metboshow$metbo_read_mbpl2))
        if(input$idMBPL2 != 'none'){
          metboshow$idM2 <- which(input$idMBPL2==colnames(metboshow$metbo_read_mbpl2))
          metboshow$metbo_data_mbpl2 <- set_input_obj_decat(metboshow$metbo_read_mbpl2,metboshow$idM2,metboshow$classM2,
                                                            metboshow$firstM2)
        }else{
          metboshow$metbo_data_mbpl2 <- set_input_obj_decat(metboshow$metbo_read_mbpl2,classCol=metboshow$classM2,
                                                            xCol=metboshow$firstM2)
        }

        if(sum(is.na(metboshow$metbo_data_mbpl1$X)) > 0 || sum(is.na(metboshow$metbo_data_mbpl2$X)) > 0  ){#Data contains missing values
          output$txtbox.mbpl <- renderPrint({
            isolate({
              if(sum(is.na(metboshow$metbo_data_mbpl1$X)) > 0){
                cat("\nUploaded data set 1 with", nrow(metboshow$metbo_read_mbpl1), "rows and", ncol(metboshow$metbo_read_mbpl1), "columns.\n","class/factor column: ",metboshow$classM1,"\n")
                if(ncol(metboshow$metbo_read_mbpl1) >2000){cat("***The data set 1 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL1 == 'none'){cat("\nMissing sample ID column in data set 1; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM1,"\n")}
                cat("\nThe data set 1 contains missing values",sum(is.na(metboshow$metbo_data_mbpl1$X)), "(",round((sum(is.na(metboshow$metbo_data_mbpl1$X))/(nrow(metboshow$metbo_data_mbpl1$X)*ncol(metboshow$metbo_data_mbpl1$X)))*100,2),"%) missing values.\n")
                cat("\nPlease impute missing values in data set 1.\n")
              }else{
                cat("\nUploaded data set 1 with", nrow(metboshow$metbo_read_mbpl1), "rows and", ncol(metboshow$metbo_read_mbpl1), "columns.\n","class/factor column: ",metboshow$classM1,"\n")
                if(ncol(metboshow$metbo_read_mbpl1) >2000){cat("***The data set 1 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL1 == 'none'){cat("\nMissing sample ID column in data set 1; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM1,"\n")}
              }
              cat("\n***************************************************************")
              if(sum(is.na(metboshow$metbo_data_mbpl2$X)) > 0){
                cat("\nUploaded data set 2 with", nrow(metboshow$metbo_read_mbpl2), "rows and", ncol(metboshow$metbo_read_mbpl2), "columns.\n","class/factor column: ",metboshow$classM2,"\n")
                if(ncol(metboshow$metbo_read_mbpl2) >2000){cat("***The data set 2 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL2 == 'none'){cat("\nMissing sample ID column in data set 2; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM2,"\n")}
                cat("\nThe data set 2 contains missing values",sum(is.na(metboshow$metbo_data_mbpl2$X)), "(",round((sum(is.na(metboshow$metbo_data_mbpl2$X))/(nrow(metboshow$metbo_data_mbpl2$X)*ncol(metboshow$metbo_data_mbpl2$X)))*100,2),"%) missing values.\n")
                cat("\nPlease impute missing values in data set 2.\n")
              }else{
                cat("\nUploaded data set 2 with", nrow(metboshow$metbo_read_mbpl2), "rows and", ncol(metboshow$metbo_read_mbpl2), "columns.\n","class/factor column: ",metboshow$classM2,"\n")
                if(ncol(metboshow$metbo_read_mbpl2) >2000){cat("***The data set 1 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL2 == 'none'){cat("\nMissing sample ID column in data set 2; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM2,"\n")}
              }
              cat("\n***************************************************************")
              cat("\nNumber of data sets for MBPLS-DA analysis: 2 data sets\n\n")
              cat("Click 'Next' to impute missing value.")
              metboshow$nextMBPL_M <- TRUE

              shinyjs::show("NextMBPL")
            })
          })
        }else{
          output$txtbox.mbpl <- renderPrint({
            isolate({
              cat("\nUploaded data set 1 with", nrow(metboshow$metbo_read_mbpl1), "rows and", ncol(metboshow$metbo_read_mbpl1), "columns.\n","class/factor column: ",metboshow$classM1,"\n")
              if(ncol(metboshow$metbo_read_mbpl1) >2000){cat("***The data set 1 is too large (>2000 column), the analysis might take very long time!!!!")}
              if(input$idMBPL1 == 'none'){cat("\nMissing sample ID column in data set 1; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM1,"\n")}
              cat("\n***************************************************************")
              cat("\nUploaded data set 2 with", nrow(metboshow$metbo_read_mbpl2), "rows and", ncol(metboshow$metbo_read_mbpl2), "columns.\n","class/factor column: ",metboshow$classM2,"\n")
              if(ncol(metboshow$metbo_read_mbpl1) >2000){cat("***The data set 1 is too large (>2000 column), the analysis might take very long time!!!!")}
              if(input$idMBPL2 == 'none'){cat("\nMissing sample ID column in data set 2; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM2,"\n")}
              cat("\n***************************************************************")
              cat("\nNumber of data sets for MBPLS-DA analysis: 2 data sets\n\n")
              cat("Click 'Next' to MBPLS-DA analysis.")
              metboshow$nextMBPL_M <- FALSE
            })
            shinyjs::show("NextMBPL")
          })
        }

      }else{
        output$txtbox.mbpl <- renderPrint({
          cat("\nSome input data set are unreadable, Plese re-upload the new file.")
        })
        shinyjs::hide("NextMBPL")
        return(NULL)
      }
    }
    ###set input object for MBPLS-DA analysis (3 data sets)
    if(input$MBPLdmun == 3){
      hideTab(inputId = "MBPLtabDat", target = "MBPLdatS4")
      hideTab(inputId = "MBPLtabDat", target = "MBPLdatS5")
      if(!is.null(metboshow$metbo_read_mbpl1) && !is.null(metboshow$metbo_read_mbpl2) && !is.null(metboshow$metbo_read_mbpl3) ){
        metboshow$classM1 <- which(input$classMBPL1==colnames(metboshow$metbo_read_mbpl1))
        metboshow$firstM1 <- which(input$firstMBPL1==colnames(metboshow$metbo_read_mbpl1))
        metboshow$nameMBPL1 <- as.character(input$nameMBPL1)
        if(input$idMBPL1 != 'none'){
          metboshow$idM1 <- which(input$idMBPL1==colnames(metboshow$metbo_read_mbpl1))
          metboshow$metbo_data_mbpl1 <- set_input_obj_decat(metboshow$metbo_read_mbpl1,metboshow$idM1,metboshow$classM1,
                                                            metboshow$firstM1)
        }else{
          metboshow$metbo_data_mbpl1 <- set_input_obj_decat(metboshow$metbo_read_mbpl1,classCol=metboshow$classM1,
                                                            xCol=metboshow$firstM1)
        }

        metboshow$classM2 <- which(input$classMBPL2==colnames(metboshow$metbo_read_mbpl2))
        metboshow$firstM2 <- which(input$firstMBPL2==colnames(metboshow$metbo_read_mbpl2))
        metboshow$nameMBPL2 <- as.character(input$nameMBPL2)
        if(input$idMBPL2 != 'none'){
          metboshow$idM2 <- which(input$idMBPL2==colnames(metboshow$metbo_read_mbpl2))
          metboshow$metbo_data_mbpl2 <- set_input_obj_decat(metboshow$metbo_read_mbpl2,metboshow$idM2,metboshow$classM2,
                                                            metboshow$firstM2)
        }else{
          metboshow$metbo_data_mbpl2 <- set_input_obj_decat(metboshow$metbo_read_mbpl2,classCol=metboshow$classM2,
                                                            xCol=metboshow$firstM2)
        }

        metboshow$classM3 <- which(input$classMBPL3==colnames(metboshow$metbo_read_mbpl3))
        metboshow$firstM3 <- which(input$firstMBPL3==colnames(metboshow$metbo_read_mbpl3))
        metboshow$nameMBPL3 <- as.character(input$nameMBPL3)
        if(input$idMBPL3 != 'none'){
          metboshow$idM3 <- which(input$idMBPL3==colnames(metboshow$metbo_read_mbpl3))
          metboshow$metbo_data_mbpl3 <- set_input_obj_decat(metboshow$metbo_read_mbpl3,metboshow$idM3,metboshow$classM3,
                                                            metboshow$firstM3)
        }else{
          metboshow$metbo_data_mbpl3 <- set_input_obj_decat(metboshow$metbo_read_mbpl3,classCol=metboshow$classM3,
                                                            xCol=metboshow$firstM3)
        }

        if(sum(is.na(metboshow$metbo_data_mbpl1$X)) > 0 || sum(is.na(metboshow$metbo_data_mbpl2$X)) > 0 || sum(is.na(metboshow$metbo_data_mbpl3$X)) > 0 ){#Data contains missing values
          output$txtbox.mbpl <- renderPrint({
            isolate({
              if(sum(is.na(metboshow$metbo_data_mbpl1$X)) > 0){
                cat("\nUploaded data set 1 with", nrow(metboshow$metbo_read_mbpl1), "rows and", ncol(metboshow$metbo_read_mbpl1), "columns.\n","class/factor column: ",metboshow$classM1,"\n")
                if(ncol(metboshow$metbo_read_mbpl1) >2000){cat("***The data set 1 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL1 == 'none'){cat("\nMissing sample ID column in data set 1; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM1,"\n")}
                cat("\nThe data set 1 contains missing values",sum(is.na(metboshow$metbo_data_mbpl1$X)), "(",round((sum(is.na(metboshow$metbo_data_mbpl1$X))/(nrow(metboshow$metbo_data_mbpl1$X)*ncol(metboshow$metbo_data_mbpl1$X)))*100,2),"%) missing values.\n")
                cat("\nPlease impute missing values in data set 1.\n")
              }else{
                cat("\nUploaded data set 1 with", nrow(metboshow$metbo_read_mbpl1), "rows and", ncol(metboshow$metbo_read_mbpl1), "columns.\n","class/factor column: ",metboshow$classM1,"\n")
                if(ncol(metboshow$metbo_read_mbpl1) >2000){cat("***The data set 1 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL1 == 'none'){cat("\nMissing sample ID column in data set 1; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM1,"\n")}
              }
              cat("\n***************************************************************")
              if(sum(is.na(metboshow$metbo_data_mbpl2$X)) > 0){
                cat("\nUploaded data set 2 with", nrow(metboshow$metbo_read_mbpl2), "rows and", ncol(metboshow$metbo_read_mbpl2), "columns.\n","class/factor column: ",metboshow$classM2,"\n")
                if(ncol(metboshow$metbo_read_mbpl2) >2000){cat("***The data set 2 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL2 == 'none'){cat("\nMissing sample ID column in data set 2; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM2,"\n")}
                cat("\nThe data set 2 contains missing values",sum(is.na(metboshow$metbo_data_mbpl2$X)), "(",round((sum(is.na(metboshow$metbo_data_mbpl2$X))/(nrow(metboshow$metbo_data_mbpl2$X)*ncol(metboshow$metbo_data_mbpl2$X)))*100,2),"%) missing values.\n")
                cat("\nPlease impute missing values in data set 2.\n")
              }else{
                cat("\nUploaded data set 2 with", nrow(metboshow$metbo_read_mbpl2), "rows and", ncol(metboshow$metbo_read_mbpl2), "columns.\n","class/factor column: ",metboshow$classM2,"\n")
                if(ncol(metboshow$metbo_read_mbpl2) >2000){cat("***The data set 2 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL2 == 'none'){cat("\nMissing sample ID column in data set 2; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM2,"\n")}
              }
              cat("\n***************************************************************")
              if(sum(is.na(metboshow$metbo_data_mbpl3$X)) > 0){
                cat("\nUploaded data set 3 with", nrow(metboshow$metbo_read_mbpl3), "rows and", ncol(metboshow$metbo_read_mbpl3), "columns.\n","class/factor column: ",metboshow$classM3,"\n")
                if(ncol(metboshow$metbo_read_mbpl3) >2000){cat("***The data set 3 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL3 == 'none'){cat("\nMissing sample ID column in data set 3; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM3,"\n")}
                cat("\nThe data set 3 contains missing values",sum(is.na(metboshow$metbo_data_mbpl3$X)), "(",round((sum(is.na(metboshow$metbo_data_mbpl3$X))/(nrow(metboshow$metbo_data_mbpl3$X)*ncol(metboshow$metbo_data_mbpl3$X)))*100,2),"%) missing values.\n")
                cat("\nPlease impute missing values in data set 3.\n")
              }else{
                cat("\nUploaded data set 3 with", nrow(metboshow$metbo_read_mbpl3), "rows and", ncol(metboshow$metbo_read_mbpl3), "columns.\n","class/factor column: ",metboshow$classM3,"\n")
                if(ncol(metboshow$metbo_read_mbpl3) >2000){cat("***The data set 3 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL3 == 'none'){cat("\nMissing sample ID column in data set 3; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM3,"\n")}
              }
              cat("\n***************************************************************")
              cat("\nNumber of data sets for MBPLS-DA analysis: 3 data sets\n\n")
              cat("Click 'Next' to impute missing value.")
              metboshow$nextMBPL_M <- TRUE
              shinyjs::show("NextMBPL")
            })
          })
        }else{
          output$txtbox.mbpl <- renderPrint({
            isolate({
              cat("\nUploaded data set 1 with", nrow(metboshow$metbo_read_mbpl1), "rows and", ncol(metboshow$metbo_read_mbpl1), "columns.\n","class/factor column: ",metboshow$classM1,"\n")
              if(ncol(metboshow$metbo_read_mbpl1) >2000){cat("***The data set 1 is too large (>2000 column), the analysis might take very long time!!!!")}
              if(input$idMBPL1 == 'none'){cat("\nMissing sample ID column in data set 1; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM1,"\n")}
              cat("\n***************************************************************")
              cat("\nUploaded data set 2 with", nrow(metboshow$metbo_read_mbpl2), "rows and", ncol(metboshow$metbo_read_mbpl2), "columns.\n","class/factor column: ",metboshow$classM2,"\n")
              if(ncol(metboshow$metbo_read_mbpl2) >2000){cat("***The data set 2 is too large (>2000 column), the analysis might take very long time!!!!")}
              if(input$idMBPL2 == 'none'){cat("\nMissing sample ID column in data set 2; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM2,"\n")}
              cat("\n***************************************************************")
              cat("\nUploaded data set 3 with", nrow(metboshow$metbo_read_mbpl3), "rows and", ncol(metboshow$metbo_read_mbpl3), "columns.\n","class/factor column: ",metboshow$classM3,"\n")
              if(ncol(metboshow$metbo_read_mbpl3) >2000){cat("***The data set 3 is too large (>2000 column), the analysis might take very long time!!!!")}
              if(input$idMBPL3 == 'none'){cat("\nMissing sample ID column in data set 3; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM3,"\n")}
              cat("\n***************************************************************")
              cat("\nNumber of data sets for MBPLS-DA analysis: 3 data sets\n\n")
              cat("Click 'Next' to MBPLS-DA analysis.")
              metboshow$nextMBPL_M <- FALSE
            })
          })
          shinyjs::show("NextMBPL")
        }
      }else{
        output$txtbox.mbpl <- renderPrint({
          cat("\nSome input data set are unreadable, Plese re-upload the new file.")
        })
        shinyjs::hide("NextMBPL")
        return(NULL)
      }
    }
    ###set input object for MBPLS-DA analysis (4 data sets)
    if(input$MBPLdmun == 4){
      hideTab(inputId = "MBPLtabDat", target = "MBPLdatS5")
      if(!is.null(metboshow$metbo_read_mbpl1) && !is.null(metboshow$metbo_read_mbpl2) && !is.null(metboshow$metbo_read_mbpl3)
         && !is.null(metboshow$metbo_read_mbpl4)){
        metboshow$classM1 <- which(input$classMBPL1==colnames(metboshow$metbo_read_mbpl1))
        metboshow$firstM1 <- which(input$firstMBPL1==colnames(metboshow$metbo_read_mbpl1))
        metboshow$nameMBPL1 <- as.character(input$nameMBPL1)
        if(input$idMBPL1 != 'none'){
          metboshow$idM1 <- which(input$idMBPL1==colnames(metboshow$metbo_read_mbpl1))
          metboshow$metbo_data_mbpl1 <- set_input_obj_decat(metboshow$metbo_read_mbpl1,metboshow$idM1,metboshow$classM1,
                                                            metboshow$firstM1)
        }else{
          metboshow$metbo_data_mbpl1 <- set_input_obj_decat(metboshow$metbo_read_mbpl1,classCol=metboshow$classM1,
                                                            xCol=metboshow$firstM1)
        }

        metboshow$classM2 <- which(input$classMBPL2==colnames(metboshow$metbo_read_mbpl2))
        metboshow$firstM2 <- which(input$firstMBPL2==colnames(metboshow$metbo_read_mbpl2))
        metboshow$nameMBPL2 <- as.character(input$nameMBPL2)
        if(input$idMBPL2 != 'none'){
          metboshow$idM2 <- which(input$idMBPL2==colnames(metboshow$metbo_read_mbpl2))
          metboshow$metbo_data_mbpl2 <- set_input_obj_decat(metboshow$metbo_read_mbpl2,metboshow$idM2,metboshow$classM2,
                                                            metboshow$firstM2)
        }else{
          metboshow$metbo_data_mbpl2 <- set_input_obj_decat(metboshow$metbo_read_mbpl2,classCol=metboshow$classM2,
                                                            xCol=metboshow$firstM2)
        }

        metboshow$classM3 <- which(input$classMBPL3==colnames(metboshow$metbo_read_mbpl3))
        metboshow$firstM3 <- which(input$firstMBPL3==colnames(metboshow$metbo_read_mbpl3))
        metboshow$nameMBPL3 <- as.character(input$nameMBPL3)
        if(input$idMBPL3 != 'none'){
          metboshow$idM3 <- which(input$idMBPL3==colnames(metboshow$metbo_read_mbpl3))
          metboshow$metbo_data_mbpl3 <- set_input_obj_decat(metboshow$metbo_read_mbpl3,metboshow$idM3,metboshow$classM3,
                                                            metboshow$firstM3)
        }else{
          metboshow$metbo_data_mbpl3 <- set_input_obj_decat(metboshow$metbo_read_mbpl3,classCol=metboshow$classM3,
                                                            xCol=metboshow$firstM3)
        }

        metboshow$classM4 <- which(input$classMBPL4==colnames(metboshow$metbo_read_mbpl4))
        metboshow$firstM4 <- which(input$firstMBPL4==colnames(metboshow$metbo_read_mbpl4))
        metboshow$nameMBPL4 <- as.character(input$nameMBPL4)
        if(input$idMBPL4 != 'none'){
          metboshow$idM4 <- which(input$idMBPL4==colnames(metboshow$metbo_read_mbpl4))
          metboshow$metbo_data_mbpl4 <- set_input_obj_decat(metboshow$metbo_read_mbpl4,metboshow$idM4,metboshow$classM4,
                                                            metboshow$firstM4)
        }else{
          metboshow$metbo_data_mbpl4 <- set_input_obj_decat(metboshow$metbo_read_mbpl4,classCol=metboshow$classM4,
                                                            xCol=metboshow$firstM4)
        }

        if(sum(is.na(metboshow$metbo_data_mbpl1$X)) > 0 || sum(is.na(metboshow$metbo_data_mbpl2$X)) > 0 || sum(is.na(metboshow$metbo_data_mbpl3$X)) > 0 || sum(is.na(metboshow$metbo_data_mbpl4$X)) > 0){#Data contains missing values
          output$txtbox.mbpl <- renderPrint({
            isolate({
              if(sum(is.na(metboshow$metbo_data_mbpl1$X)) > 0){
                cat("\nUploaded data set 1 with", nrow(metboshow$metbo_read_mbpl1), "rows and", ncol(metboshow$metbo_read_mbpl1), "columns.\n","class/factor column: ",metboshow$classM1,"\n")
                if(ncol(metboshow$metbo_read_mbpl1) >2000){cat("***The data set 1 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL1 == 'none'){cat("\nMissing sample ID column in data set 1; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM1,"\n")}
                cat("\nThe data set 1 contains missing values",sum(is.na(metboshow$metbo_data_mbpl1$X)), "(",round((sum(is.na(metboshow$metbo_data_mbpl1$X))/(nrow(metboshow$metbo_data_mbpl1$X)*ncol(metboshow$metbo_data_mbpl1$X)))*100,2),"%) missing values.\n")
                cat("\nPlease impute missing values in data set 1.\n")
              }else{
                cat("\nUploaded data set 1 with", nrow(metboshow$metbo_read_mbpl1), "rows and", ncol(metboshow$metbo_read_mbpl1), "columns.\n","class/factor column: ",metboshow$classM1,"\n")
                if(ncol(metboshow$metbo_read_mbpl1) >2000){cat("***The data set 1 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL1 == 'none'){cat("\nMissing sample ID column in data set 1; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM1,"\n")}
              }
              cat("\n***************************************************************")
              if(sum(is.na(metboshow$metbo_data_mbpl2$X)) > 0){
                cat("\nUploaded data set 2 with", nrow(metboshow$metbo_read_mbpl2), "rows and", ncol(metboshow$metbo_read_mbpl2), "columns.\n","class/factor column: ",metboshow$classM2,"\n")
                if(ncol(metboshow$metbo_read_mbpl2) >2000){cat("***The data set 2 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL2 == 'none'){cat("\nMissing sample ID column in data set 2; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM2,"\n")}
                cat("\nThe data set 2 contains missing values",sum(is.na(metboshow$metbo_data_mbpl2$X)), "(",round((sum(is.na(metboshow$metbo_data_mbpl2$X))/(nrow(metboshow$metbo_data_mbpl2$X)*ncol(metboshow$metbo_data_mbpl2$X)))*100,2),"%) missing values.\n")
                cat("\nPlease impute missing values in data set 2.\n")
              }else{
                cat("\nUploaded data set 2 with", nrow(metboshow$metbo_read_mbpl2), "rows and", ncol(metboshow$metbo_read_mbpl2), "columns.\n","class/factor column: ",metboshow$classM2,"\n")
                if(ncol(metboshow$metbo_read_mbpl2) >2000){cat("***The data set 2 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL2 == 'none'){cat("\nMissing sample ID column in data set 2; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM2,"\n")}
              }
              cat("\n***************************************************************")
              if(sum(is.na(metboshow$metbo_data_mbpl3$X)) > 0){
                cat("\nUploaded data set 3 with", nrow(metboshow$metbo_read_mbpl3), "rows and", ncol(metboshow$metbo_read_mbpl3), "columns.\n","class/factor column: ",metboshow$classM3,"\n")
                if(ncol(metboshow$metbo_read_mbpl3) >2000){cat("***The data set 3 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL3 == 'none'){cat("\nMissing sample ID column in data set 3; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM3,"\n")}
                cat("\nThe data set 3 contains missing values",sum(is.na(metboshow$metbo_data_mbpl3$X)), "(",round((sum(is.na(metboshow$metbo_data_mbpl3$X))/(nrow(metboshow$metbo_data_mbpl3$X)*ncol(metboshow$metbo_data_mbpl3$X)))*100,2),"%) missing values.\n")
                cat("\nPlease impute missing values in data set 3.\n")
              }else{
                cat("\nUploaded data set 3 with", nrow(metboshow$metbo_read_mbpl3), "rows and", ncol(metboshow$metbo_read_mbpl3), "columns.\n","class/factor column: ",metboshow$classM3,"\n")
                if(ncol(metboshow$metbo_read_mbpl3) >2000){cat("***The data set 3 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL3 == 'none'){cat("\nMissing sample ID column in data set 3; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM3,"\n")}
              }
              cat("\n***************************************************************")
              if(sum(is.na(metboshow$metbo_data_mbpl4$X)) > 0){
                cat("\nUploaded data set 4 with", nrow(metboshow$metbo_read_mbpl4), "rows and", ncol(metboshow$metbo_read_mbpl4), "columns.\n","class/factor column: ",metboshow$classM4,"\n")
                if(ncol(metboshow$metbo_read_mbpl4) >2000){cat("***The data set 4 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL4 == 'none'){cat("\nMissing sample ID column in data set 4; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM4,"\n")}
                cat("\nThe data set 4 contains missing values",sum(is.na(metboshow$metbo_data_mbpl4$X)), "(",round((sum(is.na(metboshow$metbo_data_mbpl4$X))/(nrow(metboshow$metbo_data_mbpl4$X)*ncol(metboshow$metbo_data_mbpl4$X)))*100,2),"%) missing values.\n")
                cat("\nPlease impute missing values in data set 4.\n")
              }else{
                cat("\nUploaded data set 4 with", nrow(metboshow$metbo_read_mbpl4), "rows and", ncol(metboshow$metbo_read_mbpl4), "columns.\n","class/factor column: ",metboshow$classM4,"\n")
                if(ncol(metboshow$metbo_read_mbpl4) >2000){cat("***The data set 4 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL4 == 'none'){cat("\nMissing sample ID column in data set 4; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM4,"\n")}
              }
              cat("\n***************************************************************")
              cat("\nNumber of data sets for MBPLS-DA analysis: 4 data sets\n\n")
              cat("Click 'Next' to impute missing value.")
              metboshow$nextMBPL_M <- TRUE

              shinyjs::show("NextMBPL")
            })
          })
        }else{
          output$txtbox.mbpl <- renderPrint({
            isolate({
              cat("\nUploaded data set 1 with", nrow(metboshow$metbo_read_mbpl1), "rows and", ncol(metboshow$metbo_read_mbpl1), "columns.\n","class/factor column: ",metboshow$classM1,"\n")
              if(ncol(metboshow$metbo_read_mbpl1) >2000){cat("***The data set 1 is too large (>2000 column), the analysis might take very long time!!!!")}
              if(input$idMBPL1 == 'none'){cat("\nMissing sample ID column in data set 1; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM1,"\n")}
              cat("\n***************************************************************")
              cat("\nUploaded data set 2 with", nrow(metboshow$metbo_read_mbpl2), "rows and", ncol(metboshow$metbo_read_mbpl2), "columns.\n","class/factor column: ",metboshow$classM2,"\n")
              if(ncol(metboshow$metbo_read_mbpl2) >2000){cat("***The data set 2 is too large (>2000 column), the analysis might take very long time!!!!")}
              if(input$idMBPL2 == 'none'){cat("\nMissing sample ID column in data set 2; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM2,"\n")}
              cat("\n***************************************************************")
              cat("\nUploaded data set 3 with", nrow(metboshow$metbo_read_mbpl3), "rows and", ncol(metboshow$metbo_read_mbpl3), "columns.\n","class/factor column: ",metboshow$classM3,"\n")
              if(ncol(metboshow$metbo_read_mbpl3) >2000){cat("***The data set 3 is too large (>2000 column), the analysis might take very long time!!!!")}
              if(input$idMBPL3 == 'none'){cat("\nMissing sample ID column in data set 3; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM3,"\n")}
              cat("\n***************************************************************")
              cat("\nUploaded data set 4 with", nrow(metboshow$metbo_read_mbpl4), "rows and", ncol(metboshow$metbo_read_mbpl4), "columns.\n","class/factor column: ",metboshow$classM4,"\n")
              if(ncol(metboshow$metbo_read_mbpl4) >2000){cat("***The data set 4 is too large (>2000 column), the analysis might take very long time!!!!")}
              if(input$idMBPL4 == 'none'){cat("\nMissing sample ID column in data set 4; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM4,"\n")}
              cat("\n***************************************************************")
              cat("\nNumber of data sets for MBPLS-DA analysis: 4 data sets\n\n")
              cat("Click 'Next' to MBPLS-DA analysis.")
              metboshow$nextMBPL_M <- FALSE
            })
          })
          shinyjs::show("NextMBPL")

        }

      }else{
        output$txtbox.mbpl <- renderPrint({
          cat("\nSome input data set are  unreadable, Plese re-upload the new file.")
        })
        shinyjs::hide("NextMBPL")
        return(NULL)
      }
    }
    ###set input object for MBPLS-DA analysis (5 data sets)
    if(input$MBPLdmun == 5){
      shinyjs::show("MBPLsummary")
      output$txtbox.mbpl <- renderPrint({
        cat("")
      })
      shinyjs::show("NextMBPL1")
      shinyjs::show("NextMBPL2")
      shinyjs::show("NextMBPL3")
      shinyjs::show("NextMBPL4")
      shinyjs::show("NextMBPL5")
      req(metboshow$metbo_read_mbpl1);req(metboshow$metbo_read_mbpl2);req(metboshow$metbo_read_mbpl3);req(metboshow$metbo_read_mbpl4);
      req(metboshow$metbo_read_mbpl5);
      if(!is.null(metboshow$metbo_read_mbpl1) && !is.null(metboshow$metbo_read_mbpl2) && !is.null(metboshow$metbo_read_mbpl3)
         && !is.null(metboshow$metbo_read_mbpl4) && !is.null(metboshow$metbo_read_mbpl5)){
        metboshow$classM1 <- which(input$classMBPL1==colnames(metboshow$metbo_read_mbpl1))
        metboshow$firstM1 <- which(input$firstMBPL1==colnames(metboshow$metbo_read_mbpl1))
        metboshow$nameMBPL1 <- as.character(input$nameMBPL1)
        if(input$idMBPL1 != 'none'){
          metboshow$idM1 <- which(input$idMBPL1==colnames(metboshow$metbo_read_mbpl1))
          metboshow$metbo_data_mbpl1 <- set_input_obj_decat(metboshow$metbo_read_mbpl1,metboshow$idM1,metboshow$classM1,
                                                            metboshow$firstM1)
        }else{
          metboshow$metbo_data_mbpl1 <- set_input_obj_decat(metboshow$metbo_read_mbpl1,classCol=metboshow$classM1,
                                                            xCol=metboshow$firstM1)
        }

        metboshow$classM2 <- which(input$classMBPL2==colnames(metboshow$metbo_read_mbpl2))
        metboshow$firstM2 <- which(input$firstMBPL2==colnames(metboshow$metbo_read_mbpl2))
        metboshow$nameMBPL2 <- as.character(input$nameMBPL2)
        if(input$idMBPL2 != 'none'){
          metboshow$idM2 <- which(input$idMBPL2==colnames(metboshow$metbo_read_mbpl2))
          metboshow$metbo_data_mbpl2 <- set_input_obj_decat(metboshow$metbo_read_mbpl2,metboshow$idM2,metboshow$classM2,
                                                            metboshow$firstM2)
        }else{
          metboshow$metbo_data_mbpl2 <- set_input_obj_decat(metboshow$metbo_read_mbpl2,classCol=metboshow$classM2,
                                                            xCol=metboshow$firstM2)
        }

        metboshow$classM3 <- which(input$classMBPL3==colnames(metboshow$metbo_read_mbpl3))
        metboshow$firstM3 <- which(input$firstMBPL3==colnames(metboshow$metbo_read_mbpl3))
        metboshow$nameMBPL3 <- as.character(input$nameMBPL3)
        if(input$idMBPL3 != 'none'){
          metboshow$idM3 <- which(input$idMBPL3==colnames(metboshow$metbo_read_mbpl3))
          metboshow$metbo_data_mbpl3 <- set_input_obj_decat(metboshow$metbo_read_mbpl3,metboshow$idM3,metboshow$classM3,
                                                            metboshow$firstM3)
        }else{
          metboshow$metbo_data_mbpl3 <- set_input_obj_decat(metboshow$metbo_read_mbpl3,classCol=metboshow$classM3,
                                                            xCol=metboshow$firstM3)
        }

        metboshow$classM4 <- which(input$classMBPL4==colnames(metboshow$metbo_read_mbpl4))
        metboshow$firstM4 <- which(input$firstMBPL4==colnames(metboshow$metbo_read_mbpl4))
        metboshow$nameMBPL4 <- as.character(input$nameMBPL4)
        if(input$idMBPL4 != 'none'){
          metboshow$idM4 <- which(input$idMBPL4==colnames(metboshow$metbo_read_mbpl4))
          metboshow$metbo_data_mbpl4 <- set_input_obj_decat(metboshow$metbo_read_mbpl4,metboshow$idM4,metboshow$classM4,
                                                            metboshow$firstM4)
        }else{
          metboshow$metbo_data_mbpl4 <- set_input_obj_decat(metboshow$metbo_read_mbpl4,classCol=metboshow$classM4,
                                                            xCol=metboshow$firstM4)
        }

        metboshow$classM5 <- which(input$classMBPL5==colnames(metboshow$metbo_read_mbpl5))
        metboshow$firstM5 <- which(input$firstMBPL5==colnames(metboshow$metbo_read_mbpl5))
        metboshow$nameMBPL5 <- as.character(input$nameMBPL5)
        if(input$idMBPL5 != 'none'){
          metboshow$idM5 <- which(input$idMBPL5==colnames(metboshow$metbo_read_mbpl5))
          metboshow$metbo_data_mbpl5 <- set_input_obj_decat(metboshow$metbo_read_mbpl5,metboshow$idM5,metboshow$classM5,
                                                            metboshow$firstM5)
        }else{
          metboshow$metbo_data_mbpl5 <- set_input_obj_decat(metboshow$metbo_read_mbpl5,classCol=metboshow$classM5,
                                                            xCol=metboshow$firstM5)
        }

        if(sum(is.na(metboshow$metbo_data_mbpl1$X)) > 0 || sum(is.na(metboshow$metbo_data_mbpl2$X)) > 0 || sum(is.na(metboshow$metbo_data_mbpl3$X)) > 0 || sum(is.na(metboshow$metbo_data_mbpl4$X)) > 0){#Data contains missing values
          output$txtbox.mbpl <- renderPrint({
            isolate({
              if(sum(is.na(metboshow$metbo_data_mbpl1$X)) > 0){
                cat("\nUploaded data set 1 with", nrow(metboshow$metbo_read_mbpl1), "rows and", ncol(metboshow$metbo_read_mbpl1), "columns.\n","class/factor column: ",metboshow$classM1,"\n")
                if(ncol(metboshow$metbo_read_mbpl1) >2000){cat("***The data set 1 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL1 == 'none'){cat("\nMissing sample ID column in data set 1; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM1,"\n")}
                cat("\nThe data set 1 contains missing values",sum(is.na(metboshow$metbo_data_mbpl1$X)), "(",round((sum(is.na(metboshow$metbo_data_mbpl1$X))/(nrow(metboshow$metbo_data_mbpl1$X)*ncol(metboshow$metbo_data_mbpl1$X)))*100,2),"%) missing values.\n")
                cat("\nPlease impute missing values in data set 1.\n")
              }else{
                cat("\nUploaded data set 1 with", nrow(metboshow$metbo_read_mbpl1), "rows and", ncol(metboshow$metbo_read_mbpl1), "columns.\n","class/factor column: ",metboshow$classM1,"\n")
                if(ncol(metboshow$metbo_read_mbpl1) >2000){cat("***The data set 1 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL1 == 'none'){cat("\nMissing sample ID column in data set 1; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM1,"\n")}
              }
              cat("\n***************************************************************")
              if(sum(is.na(metboshow$metbo_data_mbpl2$X)) > 0){
                cat("\nUploaded data set 2 with", nrow(metboshow$metbo_read_mbpl2), "rows and", ncol(metboshow$metbo_read_mbpl2), "columns.\n","class/factor column: ",metboshow$classM2,"\n")
                if(ncol(metboshow$metbo_read_mbpl2) >2000){cat("***The data set 2 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL2 == 'none'){cat("\nMissing sample ID column in data set 2; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM2,"\n")}
                cat("\nThe data set 2 contains missing values",sum(is.na(metboshow$metbo_data_mbpl2$X)), "(",round((sum(is.na(metboshow$metbo_data_mbpl2$X))/(nrow(metboshow$metbo_data_mbpl2$X)*ncol(metboshow$metbo_data_mbpl2$X)))*100,2),"%) missing values.\n")
                cat("\nPlease impute missing values in data set 2.\n")
              }else{
                cat("\nUploaded data set 2 with", nrow(metboshow$metbo_read_mbpl2), "rows and", ncol(metboshow$metbo_read_mbpl2), "columns.\n","class/factor column: ",metboshow$classM2,"\n")
                if(ncol(metboshow$metbo_read_mbpl2) >2000){cat("***The data set 2 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL2 == 'none'){cat("\nMissing sample ID column in data set 2; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM2,"\n")}
              }
              cat("\n***************************************************************")
              if(sum(is.na(metboshow$metbo_data_mbpl3$X)) > 0){
                cat("\nUploaded data set 3 with", nrow(metboshow$metbo_read_mbpl3), "rows and", ncol(metboshow$metbo_read_mbpl3), "columns.\n","class/factor column: ",metboshow$classM3,"\n")
                if(ncol(metboshow$metbo_read_mbpl3) >2000){cat("***The data set 3 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL3 == 'none'){cat("\nMissing sample ID column in data set 3; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM3,"\n")}
                cat("\nThe data set 3 contains missing values",sum(is.na(metboshow$metbo_data_mbpl3$X)), "(",round((sum(is.na(metboshow$metbo_data_mbpl3$X))/(nrow(metboshow$metbo_data_mbpl3$X)*ncol(metboshow$metbo_data_mbpl3$X)))*100,2),"%) missing values.\n")
                cat("\nPlease impute missing values in data set 3.\n")
              }else{
                cat("\nUploaded data set 3 with", nrow(metboshow$metbo_read_mbpl3), "rows and", ncol(metboshow$metbo_read_mbpl3), "columns.\n","class/factor column: ",metboshow$classM3,"\n")
                if(ncol(metboshow$metbo_read_mbpl3) >2000){cat("***The data set 3 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL3 == 'none'){cat("\nMissing sample ID column in data set 3; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM3,"\n")}
              }
              cat("\n***************************************************************")
              if(sum(is.na(metboshow$metbo_data_mbpl4$X)) > 0){
                cat("\nUploaded data set 4 with", nrow(metboshow$metbo_read_mbpl4), "rows and", ncol(metboshow$metbo_read_mbpl4), "columns.\n","class/factor column: ",metboshow$classM4,"\n")
                if(ncol(metboshow$metbo_read_mbpl4) >2000){cat("***The data set 4 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL4 == 'none'){cat("\nMissing sample ID column in data set 4; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM4,"\n")}
                cat("\nThe data set 4 contains missing values",sum(is.na(metboshow$metbo_data_mbpl4$X)), "(",round((sum(is.na(metboshow$metbo_data_mbpl4$X))/(nrow(metboshow$metbo_data_mbpl4$X)*ncol(metboshow$metbo_data_mbpl4$X)))*100,2),"%) missing values.\n")
                cat("\nPlease impute missing values in data set 4.\n")
              }else{
                cat("\nUploaded data set 4 with", nrow(metboshow$metbo_read_mbpl4), "rows and", ncol(metboshow$metbo_read_mbpl4), "columns.\n","class/factor column: ",metboshow$classM4,"\n")
                if(ncol(metboshow$metbo_read_mbpl4) >2000){cat("***The data set 4 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL4 == 'none'){cat("\nMissing sample ID column in data set 4; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM4,"\n")}
              }
              cat("\n***************************************************************")
              if(sum(is.na(metboshow$metbo_data_mbpl5$X)) > 0){
                cat("\nUploaded data set 5 with", nrow(metboshow$metbo_read_mbpl5), "rows and", ncol(metboshow$metbo_read_mbpl5), "columns.\n","class/factor column: ",metboshow$classM5,"\n")
                if(ncol(metboshow$metbo_read_mbpl5) >2000){cat("***The data set 5 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL5 == 'none'){cat("\nMissing sample ID column in data set 5; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM5,"\n")}
                cat("\nThe data set 5 contains missing values",sum(is.na(metboshow$metbo_data_mbpl5$X)), "(",round((sum(is.na(metboshow$metbo_data_mbpl5$X))/(nrow(metboshow$metbo_data_mbpl5$X)*ncol(metboshow$metbo_data_mbpl5$X)))*100,2),"%) missing values.\n")
                cat("\nPlease impute missing values in data set 5.\n")
              }else{
                cat("\nUploaded data set 5 with", nrow(metboshow$metbo_read_mbpl5), "rows and", ncol(metboshow$metbo_read_mbpl5), "columns.\n","class/factor column: ",metboshow$classM5,"\n")
                if(ncol(metboshow$metbo_read_mbpl5) >2000){cat("***The data set 5 is too large (>2000 column), the analysis might take very long time!!!!")}
                if(input$idMBPL5 == 'none'){cat("\nMissing sample ID column in data set 5; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM5,"\n")}
              }
              cat("\n***************************************************************")
              cat("\nNumber of data sets for MBPLS-DA analysis: 4 data sets\n\n")
              cat("Click 'Next' to impute missing value.")
              metboshow$nextMBPL_M <- TRUE

              shinyjs::show("NextMBPL")
            })
          })
        }else{
          output$txtbox.mbpl <- renderPrint({
            isolate({
              cat("\nUploaded data set 1 with", nrow(metboshow$metbo_read_mbpl1), "rows and", ncol(metboshow$metbo_read_mbpl1), "columns.\n","class/factor column: ",metboshow$classM1,"\n")
              if(ncol(metboshow$metbo_read_mbpl1) >2000){cat("***The data set 1 is too large (>2000 column), the analysis might take very long time!!!!")}
              if(input$idMBPL1 == 'none'){cat("\nMissing sample ID column in data set 1; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM1,"\n")}
              cat("\n***************************************************************")
              cat("\nUploaded data set 2 with", nrow(metboshow$metbo_read_mbpl2), "rows and", ncol(metboshow$metbo_read_mbpl2), "columns.\n","class/factor column: ",metboshow$classM2,"\n")
              if(ncol(metboshow$metbo_read_mbpl2) >2000){cat("***The data set 2 is too large (>2000 column), the analysis might take very long time!!!!")}
              if(input$idMBPL2 == 'none'){cat("\nMissing sample ID column in data set 2; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM2,"\n")}
              cat("\n***************************************************************")
              cat("\nUploaded data set 3 with", nrow(metboshow$metbo_read_mbpl3), "rows and", ncol(metboshow$metbo_read_mbpl3), "columns.\n","class/factor column: ",metboshow$classM3,"\n")
              if(ncol(metboshow$metbo_read_mbpl3) >2000){cat("***The data set 3 is too large (>2000 column), the analysis might take very long time!!!!")}
              if(input$idMBPL3 == 'none'){cat("\nMissing sample ID column in data set 3; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM3,"\n")}
              cat("\n***************************************************************")
              cat("\nUploaded data set 4 with", nrow(metboshow$metbo_read_mbpl4), "rows and", ncol(metboshow$metbo_read_mbpl4), "columns.\n","class/factor column: ",metboshow$classM4,"\n")
              if(ncol(metboshow$metbo_read_mbpl4) >2000){cat("***The data set 4 is too large (>2000 column), the analysis might take very long time!!!!")}
              if(input$idMBPL4 == 'none'){cat("\nMissing sample ID column in data set 4; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM4,"\n")}
              cat("\n***************************************************************")
              cat("\nUploaded data set 5 with", nrow(metboshow$metbo_read_mbpl5), "rows and", ncol(metboshow$metbo_read_mbpl5), "columns.\n","class/factor column: ",metboshow$classM5,"\n")
              if(ncol(metboshow$metbo_read_mbpl5) >2000){cat("***The data set 5 is too large (>2000 column), the analysis might take very long time!!!!")}
              if(input$idMBPL5 == 'none'){cat("\nMissing sample ID column in data set 5; Assuming samples are independent.\n")}else{cat("ID column: ",metboshow$idM5,"\n")}
              cat("\n***************************************************************")
              cat("\nNumber of data sets for MBPLS-DA analysis: 5 data sets\n\n")
              cat("Click 'Next' to MBPLS-DA analysis.")
              metboshow$nextMBPL_M <- FALSE
            })
          })
          shinyjs::show("NextMBPL")

        }
      }else{
        output$txtbox.mbpl <- renderPrint({
          cat("\nSome input data set are  unreadable, Plese re-upload the new file.")
        })
        shinyjs::hide("NextMBPL")
        return(NULL)
      }
    }

    isolate({
      if(!is.null(metboshow$metbo_data_mbpl1)){
        showTab(inputId = "MBPLtabDat", target = "MBPLdatS1")
        output$INtableMBPL1 <- DT::renderDataTable({ as.data.frame(metboshow$metbo_data_mbpl1$inputdata) })
      }else{
        hideTab(inputId = "MBPLtabDat", target = "MBPLdatS1")
        output$INtableMBPL1 <- NULL
      }
      if(!is.null(metboshow$metbo_data_mbpl2)){
        showTab(inputId = "MBPLtabDat", target = "MBPLdatS2")
        output$INtableMBPL2 <- DT::renderDataTable({ as.data.frame(metboshow$metbo_data_mbpl2$inputdata) })
      }else{
        hideTab(inputId = "MBPLtabDat", target = "MBPLdatS2")
        output$INtableMBPL2 <- NULL
      }
      if(!is.null(metboshow$metbo_data_mbpl3)){
        showTab(inputId = "MBPLtabDat", target = "MBPLdatS3")
        output$INtableMBPL3 <- DT::renderDataTable({ as.data.frame(metboshow$metbo_data_mbpl3$inputdata) })
      }else{
        hideTab(inputId = "MBPLtabDat", target = "MBPLdatS3")
        output$INtableMBPL3 <- NULL
      }
      if(!is.null(metboshow$metbo_data_mbpl4)){
        showTab(inputId = "MBPLtabDat", target = "MBPLdatS4")
        output$INtableMBPL4 <- DT::renderDataTable({ as.data.frame(metboshow$metbo_data_mbpl4$inputdata) })
      }else{
        hideTab(inputId = "MBPLtabDat", target = "MBPLdatS4")
        output$INtableMBPL4 <- NULL
      }
      if(!is.null(metboshow$metbo_data_mbpl5)){
        showTab(inputId = "MBPLtabDat", target = "MBPLdatS5")
        output$INtableMBPL5 <- DT::renderDataTable({ as.data.frame(metboshow$metbo_data_mbpl5$inputdata) })
      }else{
        hideTab(inputId = "MBPLtabDat", target = "MBPLdatS5")
        output$INtableMBPL5 <- NULL
      }
    })
  })


  ###run impute before MBPLS-DA####
  shinyjs::hide("next_m")
  shinyjs::hide("MIsummary_m")
  shinyjs::hide("MIoutput_m")
  hideTab(inputId = "MIoutputT_m", target = "MBPLdatS1_m")
  hideTab(inputId = "MIoutputT_m", target = "MBPLdatS2_m")
  hideTab(inputId = "MIoutputT_m", target = "MBPLdatS3_m")
  hideTab(inputId = "MIoutputT_m", target = "MBPLdatS4_m")
  hideTab(inputId = "MIoutputT_m", target = "MBPLdatS5_m")

  ##option in impute missing value page#####
  observe({
    shinyjs::hide("remPercent_m")
    req(input$reall_m)
    if (input$reall_m == FALSE){
      observeEvent(input$reall_m,{shinyjs::show("remPercent_m")
        shinyjs::show("missCheck_m")} )}
    if (input$reall_m == TRUE){
      observeEvent(input$reall_m,{shinyjs::hide("remPercent_m")
        shinyjs::hide("missCheck_m")} )}
  })

  observeEvent(input$run_m, {
    shinyjs::show("MIsummary_m")
    shinyjs::show("MIoutput_m")
    shinyjs::show("next_m")
    shinyjs::show("txtbox.missing_m")
    if(!is.null(metboshow$metbo_data_mbpl1_m)){metboshow$metbo_data_mbpl1_m <- NULL}
    if(!is.null(metboshow$metbo_data_mbpl2_m)){metboshow$metbo_data_mbpl2_m <- NULL}
    if(!is.null(metboshow$metbo_data_mbpl3_m)){metboshow$metbo_data_mbpl3_m <- NULL}
    if(!is.null(metboshow$metbo_data_mbpl4_m)){metboshow$metbo_data_mbpl4_m <- NULL}
    if(!is.null(metboshow$metbo_data_mbpl5_m)){metboshow$metbo_data_mbpl5_m <- NULL}

    shinybusy::show_modal_spinner(color = "#3CB371", text = "Computing ...")
    output$txtbox.missing_m <- renderPrint({
      isolate({
        if(!is.null(metboshow$metbo_data_mbpl1)){
          if(sum(is.na(metboshow$metbo_data_mbpl1$X)) > 0){
            metboshow$metbo_data_mbpl1_m <- impute_missing_data(metboshow$metbo_data_mbpl1, method = input$missCheck_m, removeall=input$reall_m, cutoff=input$remPercent_m )
            # shinyjs::hide("txtbox.QCnorm")
            if (input$reall == TRUE){
              cat("\nRemove all variables with missing values in data set 1.\n")
            }else{
              if(input$missCheck_m=="min"){
                cat("\nImpute missing values in data set 1 with a chosen method.\n")
              }else{
                cat("\nImpute missing values in data set 1 with a chosen method.\n")
              }
            }
          }
        }

        if(!is.null(metboshow$metbo_data_mbpl2)){
          if(sum(is.na(metboshow$metbo_data_mbpl2$X)) > 0){
            metboshow$metbo_data_mbpl2_m <- impute_missing_data(metboshow$metbo_data_mbpl2, method = input$missCheck_m, removeall=input$reall_m, cutoff=input$remPercent_m )
            # shinyjs::hide("txtbox.QCnorm")
            if (input$reall == TRUE){
              cat("\nRemove all variables with missing values in data set 2.\n")
            }else{
              if(input$missCheck_m=="min"){
                cat("\nImpute missing values in data set 2 with a chosen method.\n")
              }else{
                cat("\nImpute missing values in data set 2 with a chosen method.\n")
              }
            }
          }
        }

        if(!is.null(metboshow$metbo_data_mbpl3)){
          if(sum(is.na(metboshow$metbo_data_mbpl3$X)) > 0){
            metboshow$metbo_data_mbpl3_m <- impute_missing_data(metboshow$metbo_data_mbpl3, method = input$missCheck_m, removeall=input$reall_m, cutoff=input$remPercent_m )
            # shinyjs::hide("txtbox.QCnorm")
            if (input$reall == TRUE){
              cat("\nRemove all variables with missing values in data set 3.\n")
            }else{
              if(input$missCheck_m=="min"){
                cat("\nImpute missing values in data set 3 with a chosen method.\n")
              }else{
                cat("\nImpute missing values in data set 3 with a chosen method.\n")
              }
            }
          }
        }

        if(!is.null(metboshow$metbo_data_mbpl4)){
          if(sum(is.na(metboshow$metbo_data_mbpl4$X)) > 0){
            metboshow$metbo_data_mbpl4_m <- impute_missing_data(metboshow$metbo_data_mbpl4, method = input$missCheck_m, removeall=input$reall_m, cutoff=input$remPercent_m )
            # shinyjs::hide("txtbox.QCnorm")
            if (input$reall == TRUE){
              cat("\nRemove all variables with missing values in data set 4.\n")
            }else{
              if(input$missCheck_m=="min"){
                cat("\nImpute missing values in data set 4 with a chosen method.\n")
              }else{
                cat("\nImpute missing values in data set 4 with a chosen method.\n")
              }
            }
          }
        }

        if(!is.null(metboshow$metbo_data_mbpl5)){
          if(sum(is.na(metboshow$metbo_data_mbpl5$X)) > 0){
            metboshow$metbo_data_mbpl5_m <- impute_missing_data(metboshow$metbo_data_mbpl5, method = input$missCheck_m, removeall=input$reall_m, cutoff=input$remPercent_m )
            # shinyjs::hide("txtbox.QCnorm")
            if (input$reall == TRUE){
              cat("\nRemove all variables with missing values in data set 5.\n")
            }else{
              if(input$missCheck_m=="min"){
                cat("\nImpute missing values in data set 5 with a chosen method.\n")
              }else{
                cat("\nImpute missing values in data set 5 with a chosen method.\n")
              }
            }
          }
        }
        shinybusy::remove_modal_spinner()
      })
    })

  })

  observe({
    if(!is.null(metboshow$metbo_data_mbpl1_m)){
      showTab(inputId = "MIoutputT_m", target = "MBPLdatS1_m")
      output$INtableMBPL1_m <- DT::renderDataTable({ as.data.frame(metboshow$metbo_data_mbpl1_m$X) })
    }else{
      hideTab(inputId = "MIoutputT_m", target = "MBPLdatS1_m")
      output$INtableMBPL1_m <- NULL
    }
    if(!is.null(metboshow$metbo_data_mbpl2_m)){
      showTab(inputId = "MIoutputT_m", target = "MBPLdatS2_m")
      output$INtableMBPL2_m <- DT::renderDataTable({ as.data.frame(metboshow$metbo_data_mbpl2_m$X) })
    }else{
      hideTab(inputId = "MIoutputT_m", target = "MBPLdatS2_m")
      output$INtableMBPL2_m <- NULL
    }
    if(!is.null(metboshow$metbo_data_mbpl3_m)){
      showTab(inputId = "MIoutputT_m", target = "MBPLdatS3_m")
      output$INtableMBPL3_m <- DT::renderDataTable({ as.data.frame(metboshow$metbo_data_mbpl3_m$X) })
    }else{
      hideTab(inputId = "MIoutputT_m", target = "MBPLdatS3_m")
      output$INtableMBPL3_m <- NULL
    }
    if(!is.null(metboshow$metbo_data_mbpl4_m)){
      showTab(inputId = "MIoutputT_m", target = "MBPLdatS4_m")
      output$INtableMBPL4_m <- DT::renderDataTable({ as.data.frame(metboshow$metbo_data_mbpl4_m$X) })
    }else{
      hideTab(inputId = "MIoutputT_m", target = "MBPLdatS4_m")
      output$INtableMBPL4_m <- NULL
    }
    if(!is.null(metboshow$metbo_data_mbpl5_m)){
      showTab(inputId = "MIoutputT_m", target = "MBPLdatS5_m")
      output$INtableMBPL5_m <- DT::renderDataTable({ as.data.frame(metboshow$metbo_data_mbpl5_m$X) })
    }else{
      hideTab(inputId = "MIoutputT_m", target = "MBPLdatS5_m")
      output$INtableMBPL5_m <- NULL
    }
  })

  ###RUN MBPLS-DA#####

  #treshold #npermu #nrepret
  shinyjs::hide("MBPLnpermut")
  shinyjs::hide("MBPLthreshold")
  observe({
    metboshow$MBPLtest <- as.logical(input$MBPLtest)
    if (metboshow$MBPLtest){
      shinyjs::show("MBPLnpermut")
      shinyjs::show("MBPLthreshold")
    }else{
      shinyjs::hide("MBPLnpermut")
      shinyjs::hide("MBPLthreshold")
    }
  })

  observe({
    if(is.null(metboshow$metbo_mbpl)){shinyjs::hide("MBPLplotTabAna")}else{shinyjs::show("MBPLplotTabAna")}
  })

  metbo.mbpl <- eventReactive(input$RunMBPL,{

    if(!is.null(metboshow$metbo_mbpl)){metboshow$metbo_mbpl <- NULL}

    if(!is.null(metboshow$metbo_data_mbpl1_m)){metboshow$metbo_data_mbpl1_ana <- metboshow$metbo_data_mbpl1_m}else{if(!is.null(metboshow$metbo_data_mbpl1)){metboshow$metbo_data_mbpl1_ana <- metboshow$metbo_data_mbpl1}}
    if(!is.null(metboshow$metbo_data_mbpl2_m)){metboshow$metbo_data_mbpl2_ana <- metboshow$metbo_data_mbpl2_m}else{if(!is.null(metboshow$metbo_data_mbpl2)){metboshow$metbo_data_mbpl2_ana <- metboshow$metbo_data_mbpl2}}
    if(!is.null(metboshow$metbo_data_mbpl3_m)){metboshow$metbo_data_mbpl3_ana <- metboshow$metbo_data_mbpl3_m}else{if(!is.null(metboshow$metbo_data_mbpl3)){metboshow$metbo_data_mbpl3_ana <- metboshow$metbo_data_mbpl3}}
    if(!is.null(metboshow$metbo_data_mbpl4_m)){metboshow$metbo_data_mbpl4_ana <- metboshow$metbo_data_mbpl4_m}else{if(!is.null(metboshow$metbo_data_mbpl4)){metboshow$metbo_data_mbpl4_ana <- metboshow$metbo_data_mbpl4}}
    if(!is.null(metboshow$metbo_data_mbpl5_m)){metboshow$metbo_data_mbpl5_ana <- metboshow$metbo_data_mbpl5_m}else{if(!is.null(metboshow$metbo_data_mbpl5)){metboshow$metbo_data_mbpl5_ana <- metboshow$metbo_data_mbpl5}}

    if(input$MBPLdmun == 2){metboshow$list_data_mbpl_check <- list(metboshow$metbo_data_mbpl1_ana$inputdata,metboshow$metbo_data_mbpl2_ana$inputdata)}
    if(input$MBPLdmun == 3){metboshow$list_data_mbpl_check <- list(metboshow$metbo_data_mbpl1_ana$inputdata,metboshow$metbo_data_mbpl2_ana$inputdata,metboshow$metbo_data_mbpl3_ana$inputdata)}
    if(input$MBPLdmun == 4){metboshow$list_data_mbpl_check <- list(metboshow$metbo_data_mbpl1_ana$inputdata,metboshow$metbo_data_mbpl2_ana$inputdata,metboshow$metbo_data_mbpl3_ana$inputdata,metboshow$metbo_data_mbpl4_ana$inputdata)}
    if(input$MBPLdmun == 5){metboshow$list_data_mbpl_check <- list(metboshow$metbo_data_mbpl1_ana$inputdata,metboshow$metbo_data_mbpl2_ana$inputdata,metboshow$metbo_data_mbpl3_ana$inputdata,metboshow$metbo_data_mbpl4_ana$inputdata,metboshow$metbo_data_mbpl5_ana$inputdata)}

    for(i in 1:length(metboshow$list_data_mbpl_check)){
      for (j in 1:length(metboshow$list_data_mbpl_check)){
        if(i != j){
          #cat(i,"\t",j,"\n")
          #cat(identical(metboshow$list_data_mbpl_check[[i]],metboshow$list_data_mbpl_check[[j]]),"\n")
          if(identical(metboshow$list_data_mbpl_check[[i]],metboshow$list_data_mbpl_check[[j]]) == TRUE){
            cat("\nSome data sets are identical please check your input data.\n")
            cat("\nData was not analyzed.\n")
            shinyjs::hide("MBPLplotTabAna")
            return(NULL)
          }
        }
      }
    }

    if(input$MBPLdmun == 2){
      metboshow$list_data_mbpl <- list(metboshow$metbo_data_mbpl1_ana$X,metboshow$metbo_data_mbpl2_ana$X)
      names(metboshow$list_data_mbpl) <- c(metboshow$nameMBPL1,metboshow$nameMBPL2)
    }
    if(input$MBPLdmun == 3){
      metboshow$list_data_mbpl <- list(metboshow$metbo_data_mbpl1_ana$X,metboshow$metbo_data_mbpl2_ana$X,metboshow$metbo_data_mbpl3_ana$X)
      names(metboshow$list_data_mbpl) <- c(metboshow$nameMBPL1,metboshow$nameMBPL2,metboshow$nameMBPL3)
    }
    if(input$MBPLdmun == 4){
      metboshow$list_data_mbpl <- list(metboshow$metbo_data_mbpl1_ana$X,metboshow$metbo_data_mbpl2_ana$X,metboshow$metbo_data_mbpl3_ana$X,metboshow$metbo_data_mbpl4_ana$X)
      names(metboshow$list_data_mbpl) <- c(metboshow$nameMBPL1,metboshow$nameMBPL2,metboshow$nameMBPL3,metboshow$nameMBPL4)
    }
    if(input$MBPLdmun == 5){
      metboshow$list_data_mbpl <- list(metboshow$metbo_data_mbpl1_ana$X,metboshow$metbo_data_mbpl2_ana$X,metboshow$metbo_data_mbpl3_ana$X,metboshow$metbo_data_mbpl4_ana$X,metboshow$metbo_data_mbpl5_ana$X)
      names(metboshow$list_data_mbpl) <- c(metboshow$nameMBPL1,metboshow$nameMBPL2,metboshow$nameMBPL3,metboshow$nameMBPL4,metboshow$nameMBPL5)
    }

    if(metboshow$MBPLtest){
      shinybusy::show_modal_spinner(color = "#3CB371", text = "Computing and performing testing, this process might take long time ...")
    }else{
      shinybusy::show_modal_spinner(color = "#3CB371", text = "Computing ...")
    }
    #capture >2000 column
    #cat(input$varRatio,"\t",input$partitionSize,"\t")
    if(!is.null(metboshow$list_data_mbpl)){
      if(input$MBPLnrepet < 2 ){ #base function use 1.5
        metboshow$MBPLnrepet <- 10
        cat("Number of repret paramerter was set lower than 2","\tThe default paramerter is used: 'Number of repret = 10'")
      }else{
        metboshow$MBPLnrepet <- as.numeric(input$MBPLnrepet)
      }

      if(input$MBPLthreshold > 1 ){
        metboshow$MBPLthreshold <- 0.5
        cat("Threshold paramerter was set higher than 1.","\tThe default paramerter is used: 'threshold = 0.5'")
      }else{
        metboshow$MBPLthreshold <- as.numeric(input$MBPLthreshold)
      }

      metboshow$metbo_mbpl = tryCatch({
        mbplsda_analyze(data.frame(metboshow$metbo_data_mbpl1_ana$Y), metboshow$list_data_mbpl,nf=as.numeric(input$MBPLnf),optdim=as.numeric(input$MBPLoptdim),nrepet=metboshow$MBPLnrepet, npermut= as.numeric(input$MBPLnpermut), nboot= as.numeric(input$MBPLnboot), testmodel = metboshow$MBPLtest, threshold=metboshow$MBPLthreshold)
      },error=function(e){
        shinybusy::remove_modal_spinner()
        message(e)
        return()
      })
      cat("\nDone!!!\n")
    }else{
      cat("\nERROR! Data was not analyzed.\n")
      return()
      shinybusy::remove_modal_spinner()
    }
    shinybusy::remove_modal_spinner()
  })

  observe({
    if(!is.null(metboshow$metbo_mbpl)){
      updateSelectInput(session, "MBPLselectPC",
                        choices=c(1,1:length(metboshow$metbo_mbpl$result$res_boot$faX)))
    }
  })

  observe({
    req(metboshow$metbo_mbpl)
    if (!is.null(metboshow$metbo_mbpl)){
      isolate({
        output$plotScreeMBPL <- renderPlotly({ggplotly(mbplsda_screeplot(metboshow$metbo_mbpl$result$base_model$eig, plot_title="Scree plot"))})
        if(metboshow$MBPLtest){
          output$plotTestdimMBPL <- renderPlotly({ggplotly(mbplsda_plottestdim(metboshow$metbo_mbpl$result$res_optimal) )})
          output$plotPermutMBPL <- renderPlotly({ggplotly(mbplsda_plotpermut(metboshow$metbo_mbpl$result$res_permut, plot_title="Permutation testing plot"))})
        }else{
          output$plotTestdimMBPL <- renderPlotly({plotly_empty()})
          output$plotPermutMBPL <- renderPlotly({plotly_empty()})
        }
        output$plotBipMBPL <- renderPlotly({ggplotly(mbplsda_plotboot_bipc(metboshow$metbo_mbpl$result$res_boot, plot_title="Block importance plot"))})
        output$plotVipMBPL <- renderPlotly({ggplotly(mbplsda_plotboot_vipc(metboshow$metbo_mbpl$result$res_boot,  propbestvar = input$MBPLselectvipprop/100, plot_title="Variable importance plot"))})
        output$plotLoadMBPL <- renderPlotly({ggplotly(mbplsda_plotboot_loading(metboshow$metbo_mbpl$result$res_boot, pc= as.numeric(input$MBPLselectPC), propbestvar = input$MBPLselectloadprop/100, plot_title="Loading plot by PC"))})
      })
    }else{
      output$plotScreeMBPL <- renderPlotly({plotly_empty()})
      output$plotTestdimMBPL <- renderPlotly({plotly_empty()})
      output$plotPermutMBPL <- renderPlotly({plotly_empty()})
      output$plotBipMBPL <- renderPlotly({plotly_empty()})
      output$plotVipMBPL <- renderPlotly({plotly_empty()})
      output$plotLoadMBPL <- renderPlotly({plotly_empty()})
    }
  })

  output$txtbox.mbplAna <- renderPrint({
    if(is.null(metbo.mbpl())){
      cat("")
    }
    if(is.null(metboshow$metbo_mbpl)){
      cat("")
    }
  })

  #####DE page##########

  ####run overerrep analysis#####
  observeEvent(input$revOVER, {
    if(!is.null(input$infileOver)){
      if(file.exists(input$infileOver$datapath)){ file.remove(input$infileOver$datapath)}
    }
    metboshow$infileOver <- NULL
    metboshow$EntityIDs <- NULL
    metboshow$metbo_overr <- NULL
    shinyjs::enable("EntityIDs")
    reset('EntityIDs')
    reset('INfileOVER')
    reset('OVERnode')
    reset('OVERset')
    reset('OVERsize')
    shinyjs::hide('OVERsummary')
    shinyjs::hide('OVERoutput')
  })
  shinyjs::hide('OVERoutput')
  shinyjs::hide('OVERsummary')

  observe({
    req(input$INfileOVER)
    shinyjs::show('OVERsummary')
    isolate({
      if(file.exists(input$INfileOVER$datapath)){
        if(tools::file_ext(input$INfileOVER$datapath) != "csv"){
          output$txtbox.Overr <- renderPrint({
            cat("ERROR: File type is not supported, please use 'csv' file as input\n")
          })
          metboshow$infileOver <- NULL
          if(!is.null(input$infileOver)){
            if(file.exists(input$infileOver$datapath)){ file.remove(input$infileOver$datapath)}
          }
          reset('EntityIDs')
          reset('INfileOVER')
          reset('OVERnode')
          reset('OVERset')
          reset('OVERsize')
          shinyjs::enable("EntityIDs")
          return(metboshow$metbo_overr <- NULL)
        }
        metboshow$infileOver <- read.csv(input$INfileOVER$datapath, header = F)
        updateTextAreaInput(session,"EntityIDs",value="")
        shinyjs::disable("EntityIDs")
      }
    })
  })

  observe({
    if(!is.null(metboshow$infileOver)){
      output$txtbox.Overr <- renderPrint({
        cat("Uploaded data with", nrow(metboshow$infileOver), "variables.\n")
      })
      updateTextAreaInput(session,"EntityIDs",value="")
      shinyjs::disable("EntityIDs")
    }else{
      shinyjs::enable("EntityIDs")
    }
    if(input$EntityIDs!=""){
      metboshow$EntityIDs <- strsplit(as.character(input$EntityIDs), split = "\n")
      output$txtbox.Overr <- renderPrint({
        cat("Uploaded data with", length(metboshow$EntityIDs[[1]]), "variables.\n")
      })
      shinyjs::hide("INfileOVER")
      reset('INfileOVER')
      if(!is.null(input$infileOver)){
        if(file.exists(input$infileOver$datapath)){ file.remove(input$infileOver$datapath)}
      }
    }else{
      shinyjs::show("INfileOVER")
    }
    if(is.null(metboshow$infileOver) && input$EntityIDs==""){
      shinyjs::disable("runOVER")
    }else{
      shinyjs::enable("runOVER")
    }
  })

  observeEvent(input$runOVER, {
    if(!is.null(metboshow$infileOver)){
      shinybusy::show_modal_spinner(color = "#3CB371", text = "Computing ...")
      metboshow$metbo_overr <- overrep_analyze(metboshow$infileOver, nodetype= input$OVERnode, settype=input$OVERset, organism="hsa", size= input$OVERsize)
      output$txtbox.Overr <- renderPrint({
        cat(str(metboshow$infileOver),"Returning ",nrow(metboshow$metbo_overr$enrichment)," overrepresented sets ...\n")
      })
    }
    else if (input$EntityIDs!=""){
      shinybusy::show_modal_spinner(color = "#3CB371", text = "Computing ...")
      # metboshow$EntityIDs <- strsplit(as.character(input$EntityIDs), split = "\n")
      metboshow$metbo_overr <- overrep_analyze(metboshow$EntityIDs, nodetype= input$OVERnode, settype=input$OVERset, organism="hsa", size= input$OVERsize)
      output$txtbox.Overr <- renderPrint({
        cat(str(metboshow$EntityIDs),"Returning ",nrow(metboshow$metbo_overr$enrichment)," overrepresented sets ...\n")
      })
    }else(
      return(NULL)
    )
    shinybusy::remove_modal_spinner()
    # if(is.null(input$INfileOVER) & input$EntityIDs==""){
    #   output$txtbox.Overr <- renderPrint({
    #     cat('No Input Found!!')
    #   })
    # }
    # else{return(NULL)}
  })

  observe({
    if(!is.null(metboshow$metbo_overr)){
      shinyjs::show('OVERsummary')
      shinyjs::show('OVERoutput')
      output$INtableOVER <- DT::renderDataTable({data.frame(metboshow$metbo_overr$enrichment)})
    }else{
      shinyjs::hide('OVERoutput')
      output$INtableOVER <- NULL
    }
  })

  ####run Enrichment analysis#####
  observeEvent(input$revENR, {
    if(!is.null(input$infileOver)){
      if(file.exists(input$infileOver$datapath)){ file.remove(input$infileOver$datapath)}
    }
    metboshow$metbo_ENRfile <- NULL
    metboshow$metbo_enr <- NULL
    reset('INfileENR')
    reset('ENRfccol')
    reset('ENRpcol')
    reset('ENRset')
    reset('ENRsize')
    reset('ENRmed')
    reset('ENRnode')
    shinyjs::hide('ENRsummary')
    shinyjs::hide('ENRoutput')
  })

  shinyjs::hide('ENRsummary')
  shinyjs::hide('ENRoutput')

  observe({
    req(input$INfileENR)
    shinyjs::show('ENRsummary')
    isolate({
      if(file.exists(input$INfileENR$datapath)){
        if(tools::file_ext(input$INfileENR$datapath) != "csv"){
          output$txtbox.ENR <- renderPrint({
            cat("ERROR: File type is not supported, please use 'csv' file as input\n")
          })
          metboshow$metbo_ENRfile <- NULL
          reset('INfileENR')
          reset('ENRfccol')
          reset('ENRpcol')
          return(metboshow$metbo_enr <- NULL)
        }
        metboshow$metbo_ENRfile <- data.frame(read.csv(input$INfileENR$datapath, header = T))
        if(!is.null(metboshow$metbo_ENRfile)){
          if(ncol(metboshow$metbo_ENRfile) < 2){
            output$txtbox.ENR <- renderPrint({
              cat("ERROR: Input file must contain at least 2 columns\n")
            })
            metboshow$metbo_ENRfile <- NULL
            reset('INfileENR')
            reset('ENRfccol')
            reset('ENRpcol')
            return(metboshow$metbo_enr <- NULL)
          }
          if(ncol(metboshow$metbo_ENRfile) == 2){
            output$txtbox.ENR <- renderPrint({
              cat("Uploaded data with", nrow(metboshow$metbo_ENRfile), "rows and", ncol(metboshow$metbo_ENRfile), "columns.\n")
            })
            updateSelectInput(session, "ENRpcol",
                              choices=c(colnames(metboshow$metbo_ENRfile[2]),colnames(metboshow$metbo_ENRfile)))
            updateSelectInput(session, "ENRfccol", choices=c("None" = "FALSE"))
          }
          if(ncol(metboshow$metbo_ENRfile) > 2){
            output$txtbox.ENR <- renderPrint({
              cat("Uploaded data with", nrow(metboshow$metbo_ENRfile), "rows and", ncol(metboshow$metbo_ENRfile), "columns.\n")
            })
            updateSelectInput(session, "ENRpcol",
                              choices=c(colnames(metboshow$metbo_ENRfile[2]),colnames(metboshow$metbo_ENRfile)))
            updateSelectInput(session, "ENRfccol",
                                choices=c("None" = "FALSE",colnames(metboshow$metbo_ENRfile)))
          }
        }
      }
    })
  })

  observe({
    if(!is.null(metboshow$metbo_ENRfile)){
      shinyjs::enable("runENR")
    }
    if(is.null(metboshow$metbo_ENRfile)){
      shinyjs::disable("runENR")
    }
  })

  observeEvent(input$runENR,{
    shinyjs::show('ENRsummary')
    if(!is.null(metboshow$metbo_enr)){metboshow$metbo_enr <- NULL}
    shinybusy::show_modal_spinner(color = "#3CB371", text = "Computing ...")
    if(!is.null(metboshow$metbo_ENRfile)){
      metboshow$ENRpcol <- which(input$ENRpcol==colnames(metboshow$metbo_ENRfile))
      if (input$ENRfccol == "FALSE"){
        metboshow$ENRfccol <- NULL
      }else{
        metboshow$ENRfccol <- which(input$ENRfccol==colnames(metboshow$metbo_ENRfile))
      }
      metboshow$metbo_enr <- enrichment_analyze(metboshow$metbo_ENRfile, pcol = metboshow$ENRpcol, fccol = metboshow$ENRfccol,
                                                method = input$ENRmed,nodetype= input$ENRnode, settype=input$ENRset, organism="hsa", size= input$ENRsize)
      output$txtbox.ENR <- renderPrint({
        cat(str(metboshow$metbo_ENRfile),"Returning ",nrow(metboshow$metbo_enr$enrichment)," enrichment sets ...\n")
      })
    }else{
      return(NULL)
    }
    shinybusy::remove_modal_spinner()
  })

  observe({
    req(metboshow$metbo_enr)
    if(!is.null(metboshow$metbo_enr)){
      shinyjs::show('ENRsummary')
      shinyjs::show('ENRoutput')
      output$INtableENR <- DT::renderDataTable({data.frame(metboshow$metbo_enr$enrichment)})
    }else{
      shinyjs::hide('ENRoutput')
      output$INfileENR <- NULL
    }
  })

  ####run combine ORA#####
  #error reset is not clear value in run
  # shinyjs::disable("runORA")
  # observe({
  #   if(!is.null(input$INfileORA)){shinyjs::enable("runORA")}else{shinyjs::disable("runORA")}
  # })

  observeEvent(input$revORA, {
    if(!is.null(input$INfileORA)){
      if(file.exists(input$INfileORA$datapath)){ file.remove(input$INfileORA$datapath)}
    }
    metboshow$metbo_cORA_data <- NULL
    metboshow$metbo_cORA <- NULL
    reset('INfileORA')
    reset('ORAsize')
    shinyjs::hide('ORAsummary')
    shinyjs::hide('ORAoutput')
  })

  shinyjs::hide('ORAsummary')
  shinyjs::hide('ORAoutput')

  observe({
    req(input$INfileORA)
    shinyjs::show('ORAsummary')
    isolate({
      if(file.exists(input$INfileORA$datapath)){
        if(tools::file_ext(input$INfileORA$datapath) != "csv"){
          output$txtbox.ORA <- renderPrint({
            cat("ERROR: File type is not supported, please use 'csv' file as input\n")
          })
          if(!is.null(input$INfileORA)){
            if(file.exists(input$INfileORA$datapath)){ file.remove(input$INfileORA$datapath)}
          }
          metboshow$metbo_cORA_data <- NULL
          reset('INfileORA')
          reset('ORAsize')
          return(metboshow$metbo_cORA <- NULL)
        }
        metboshow$metbo_cORA_data <- read.csv(input$INfileORA$datapath, header = TRUE)
        if(ncol(metboshow$metbo_cORA_data)>=2){
          output$txtbox.ORA <- renderPrint({
            cat("Uploaded data with", nrow(metboshow$metbo_cORA_data), "rows and", ncol(metboshow$metbo_cORA_data), "columns.\n")
          })
        }else{
          output$txtbox.ORA <- renderPrint({
            cat("ERROR: Input file must contain at least 2 columns\n")
          })
          if(!is.null(input$INfileORA)){
            if(file.exists(input$INfileORA$datapath)){ file.remove(input$INfileORA$datapath)}
          }
          metboshow$metbo_cORA_data <- NULL
          reset('INfileORA')
          reset('ORAsize')
          return(metboshow$metbo_cORA <- NULL)
        }
      }
    })
  })

  observe({
    if(!is.null(metboshow$metbo_cORA_data)){
      shinyjs::enable("runORA")
    }
    if(is.null(metboshow$metbo_cORA_data)){
      shinyjs::disable("runORA")
    }
  })

  observeEvent(input$runORA,{
    if(!is.null(metboshow$metbo_cORA_data)){
      shinybusy::show_modal_spinner(color = "#3CB371", text = "Computing ...")
      metboshow$metbo_cORA <- tryCatch({#catch error
        comb_overrep_analyze(data.frame(metboshow$metbo_cORA_data[,1:2]),size = input$ORAsize)
      },
      error=function(e){
        cat(e$message)
        cat(str(metboshow$metbo_cORA_data),"Data was not analyzed. please check data format.\n")
        NULL
      })
      if(!is.null(metboshow$metbo_cORA)){
        output$txtbox.ORA <- renderPrint({
          cat(str(metboshow$metbo_cORA_data),"Returning ",nrow(metboshow$metbo_cORA$enrichment)," overrepresented sets ...\n")
        })
      }else{
        output$txtbox.ORA <- renderPrint({
          cat(str(metboshow$metbo_cORA_data),"Data was not analyzed. please check data format.\n")
        })
      }
      shinybusy::remove_modal_spinner()
    }
  })

  observe({
    if(!is.null(metboshow$metbo_cORA)){
      shinyjs::show('ORAsummary')
      shinyjs::show('ORAoutput')
      output$INtableORA <- DT::renderDataTable({data.frame(metboshow$metbo_cORA$enrichment)})
    }else{
      shinyjs::hide('ORAoutput')
      output$INtableORA <- NULL
    }
  })

  ######export results######

  observe({
    if(is.null(metboshow$metbo_missing)){shinyjs::hide('export_mis')}else{shinyjs::show('export_mis')};
    if(is.null(metboshow$metbo_QCnorm)){shinyjs::hide('export_normIS')}else{shinyjs::show('export_normIS')};
    if(is.null(metboshow$keepValueMP)){shinyjs::hide('export_normDAT')}else{shinyjs::show('export_normDAT')};
    if(is.null(metboshow$metbo_uni)){shinyjs::hide('export_uni')}else{shinyjs::show('export_uni')};
    if(is.null(metboshow$metbo_multi)){shinyjs::hide('export_multi')}else{shinyjs::show('export_multi')};
    if(is.null(metboshow$metbo_corr)){shinyjs::hide('export_corr')}else{shinyjs::show('export_corr')};
    if(is.null(metboshow$metbo_lme)){shinyjs::hide('export_lme')}else{shinyjs::show('export_lme')};
    if(is.null(metboshow$metbo_ml)){shinyjs::hide('export_muvr')}else{shinyjs::show('export_muvr')};
    if(is.null(metboshow$metbo_overr)){shinyjs::hide('export_opa')}else{shinyjs::show('export_opa')};
    if(is.null(metboshow$metbo_enr)){shinyjs::hide('export_enr')}else{shinyjs::show('export_enr')};
    if(is.null(metboshow$metbo_cORA)){shinyjs::hide('export_cora')}else{shinyjs::show('export_cora')};
    if(is.null(metboshow$metbo_mbpl)){shinyjs::hide('export_mbpls')}else{shinyjs::show('export_mbpls')};
  })

  ######export results functions######

  ex_mis <- function(){
    if(!is.null(metboshow$metbo_missing)){
      generate_report(metboshow$metbo_missing,reportfile="impute_data_report")
      write.csv(cbind(metboshow$metbo_missing$inputdata[,1:metboshow$firstV-1],metboshow$metbo_missing$X),"imputed_output_table.csv")
      files <- c("impute_data_report.pdf","imputed_output_table.csv")
      return(files)
    }else(return(NULL))
  }

  ex_ISnorm <- function(){
    if(!is.null(metboshow$metbo_QCnorm)){
      generate_report(metboshow$metbo_QCnorm,reportfile="normalize_byqc_report")
      write.csv(cbind(metboshow$metbo_QCnorm$inputdata[,1:metboshow$firstV-1],metboshow$metbo_QCnorm$X),"ISQC_normalization_output_table.csv")
      dt_plot = list();
      dt_plot[['pca1']] = pcaplot_overview(metboshow$keepValueN, scale=FALSE, plot_title="Before normalization")
      dt_plot[['pca2']] = pcaplot_overview(metboshow$metbo_QCnorm,scale=FALSE, plot_title="After normalization")
      dt_plot[['rla1']] = rlaplot_overview(metboshow$keepValueN, limitx = TRUE, dolog = TRUE, plot_title="Before normalization")
      dt_plot[['rla2']] = rlaplot_overview(metboshow$metbo_QCnorm, limitx = TRUE, dolog = TRUE, plot_title="After normalization")
      # dt_plot[['box1']] = boxplot_overview(metboshow$keepValueN, plot_title="Before normalization")
      # dt_plot[['box2']] = boxplot_overview(metboshow$metbo_QCnorm, plot_title="After normalization")
      dt_plot[['dens1']] = densityplot_overview(metboshow$keepValueN$X, plot_title="Before normalization (variables)")
      dt_plot[['dens2']] = densityplot_overview(metboshow$metbo_QCnorm$X, plot_title="After normalization (variables)")
      dt_plot[['dens_samp1']] = densityplot_overview(metboshow$keepValueN$X, plotvar=FALSE, plot_title="Before normalization (samples)")
      dt_plot[['dens_samp2']] = densityplot_overview(metboshow$metbo_QCnorm$X, plotvar=FALSE, plot_title="After normalization (samples)")
      ggsave("ISQC_normalization_output_figure.pdf",
             marrangeGrob(grobs = dt_plot, nrow=1, ncol=2, as.table=TRUE), dpi = 600, units = "in", width = 18, height = 10, device = "pdf")
      files <- c("normalize_byqc_report.pdf","ISQC_normalization_output_table.csv","ISQC_normalization_output_figure.pdf")
      return(files)
    }else{
      return(NULL)
    }
  }

  ex_DATnorm <- function(){
    if(!is.null(metboshow$keepValueMP)){
      dt_plot = list();
      if(!is.null(metboshow$metbo_norm)){
        generate_report(metboshow$metbo_norm,reportfile="normalize_bydata_report")
      }
      if(!is.null(metboshow$metbo_tran)){
        generate_report(metboshow$metbo_tran,reportfile="transform_data_report")
      }
      write.csv(cbind(metboshow$keepValueMP$inputdata[,1:metboshow$firstV-1],metboshow$keepValueMP$X),"Data_normalization_output_table.csv")
      dt_plot[['pca1']] = pcaplot_overview(metboshow$keepValueDP,scale=FALSE, plot_title="Before normalization")
      if(!is.null(metboshow$metbo_scal)){
        generate_report(metboshow$metbo_scal,reportfile="scale_data_report")
        dt_plot[['pca2']] = pcaplot_overview(metboshow$keepValueMP,scale=FALSE, plot_title="After normalization")
      }
      if(is.null(metboshow$metbo_scal)){
        dt_plot[['pca2']] = pcaplot_overview(metboshow$keepValueMP,scale=FALSE, plot_title="After normalization")
      }
      dt_plot[['rla1']] = rlaplot_overview(metboshow$keepValueDP, limitx = TRUE, dolog = TRUE, plot_title="Before normalization")
      if(!is.null(metboshow$metbo_tran)|!is.null(metboshow$metbo_scal)){
        dt_plot[['rla2']] = rlaplot_overview(metboshow$keepValueMP, limitx = TRUE, dolog = FALSE, plot_title="After normalization")
      }else{
        dt_plot[['rla2']] = rlaplot_overview(metboshow$keepValueMP, limitx = TRUE, dolog = TRUE, plot_title="After normalization")#dolog for normalization
      }
      # dt_plot[['box1']] = boxplot_overview(metboshow$keepValueDP, plot_title="Before normalization")
      # dt_plot[['box2']] = boxplot_overview(metboshow$keepValueMP, plot_title="After normalization")
      dt_plot[['dens1']] = densityplot_overview(metboshow$keepValueDP$X, plot_title="Before normalization (variables)")
      dt_plot[['dens2']] = densityplot_overview(metboshow$keepValueMP$X, plot_title="After normalization (variables)")
      dt_plot[['dens_samp1']] = densityplot_overview(metboshow$keepValueDP$X, plotvar=FALSE, plot_title="Before normalization (samples)")
      dt_plot[['dens_samp2']] = densityplot_overview(metboshow$keepValueMP$X, plotvar=FALSE, plot_title="After normalization (samples)")
      ggsave("Data_normalization_output_figure.pdf",
             marrangeGrob(grobs = dt_plot, nrow=1, ncol=2, as.table=TRUE), dpi = 600, units = "in", width = 18, height = 10, device = "pdf")
      files <- c("Data_normalization_output_table.csv","Data_normalization_output_figure.pdf")
      if(!is.null(metboshow$metbo_norm)){
        files <- c(files,"normalize_bydata_report.pdf")
      }
      if(!is.null(metboshow$metbo_tran)){
        files <- c(files,"transform_data_report.pdf")
      }
      if(!is.null(metboshow$metbo_scal)){
        files <- c(files,"scale_data_report.pdf")
      }
      return(files)
    }else{
      return(NULL)
    }
  }

  ######export mis norm######
  output$export_mis <- downloadHandler(
    filename = function(){
      paste("Metabox_report","tar",sep=".")
    },
    content = function(file){
      shinybusy::show_modal_spinner(color = "#3CB371", text = "Exporting results ...")
      files_miss <- ex_mis()
      shinybusy::remove_modal_spinner()
      tar(tarfile= file, files = files_miss)
    },
    contentType = "application/tar"
  )

  ######export IS norm######
  output$export_normIS <- downloadHandler(
    filename = function(){
      paste("Metabox_report","tar",sep=".")
    },
    content = function(file){
      shinybusy::show_modal_spinner(color = "#3CB371", text = "Exporting results ...")
      files_miss <- ex_mis()
      files_IS <- ex_ISnorm()
      shinybusy::remove_modal_spinner()
      tar(tarfile= file, files = c(files_miss,files_IS))
    },
    contentType = "application/tar"
  )

  ######export DATA norm######
  output$export_normDAT <- downloadHandler(
    filename = function(){
      paste("Metabox_report","tar",sep=".")
    },
    content = function(file){
      shinybusy::show_modal_spinner(color = "#3CB371", text = "Exporting results ...")
      files_miss <- ex_mis()
      files_IS <- ex_ISnorm()
      files_DATnorm <- ex_DATnorm()
      shinybusy::remove_modal_spinner()
      tar(tarfile= file, files = c(files_miss,files_IS,files_DATnorm))
    },
    contentType = "application/tar"
  )

  ######export uni analyze######
  ex_uni <- function(){
    if(!is.null(metboshow$metbo_uni)){
      generate_report(metboshow$keepValueM,datsummary1=metboshow$metbo_uni,reportfile="univariate_report")
      f1inx = grep("adjusted_pvalue", colnames(metboshow$univ_merge_table))[1]
      pvtb = data.frame(metboshow$univ_merge_table[,f1inx], row.names = metboshow$univ_merge_table$variable)
      pvtb = pvtb[order(pvtb),, drop=FALSE]
      filtpv = head(pvtb,n=100)
      ggsave("Univariate_analysis_output_figure.pdf", pvalplot_overview(filtpv, plot_title = "Statistical significance plot"), dpi = 600, units = "in", device = "pdf", width = 12, height = 6)
      if(metboshow$metbo_uni$details$posthocTest == "No test"){
        write.csv(metboshow$univ_merge_table,"Univariate_analysis_output_table.csv",row.names = T)
      }
      else{
        write.csv(metboshow$poshoc_merge_table,"Univariate_analysis_output_table.csv",row.names = T)
      }
      files <- c("Univariate_analysis_output_figure.pdf","univariate_report.pdf","Univariate_analysis_output_table.csv")
      return(files)
    }else{
      return(NULL)
    }
  }

  output$export_uni <- downloadHandler(
    filename = function(){
      paste("Metabox_report","tar",sep=".")
    },
    content = function(file){
      shinybusy::show_modal_spinner(color = "#3CB371", text = "Exporting results ...")

      if(metboshow$full == TRUE){
        files_miss <- ex_mis()
        files_IS <- ex_ISnorm()
        files_DATnorm <- ex_DATnorm()
      }else{
        files_miss <- NULL
        files_IS <- NULL
        files_DATnorm <- NULL
      }
      if(metboshow$S_missing_M1 == TRUE){
        files_miss <- ex_mis()
      }else{
        files_miss <- NULL
      }
      files_uni <- ex_uni()
      files_multi <- ex_multi()
      files_corr <- ex_corr()
      files_lme <- ex_lme()
      shinybusy::remove_modal_spinner()
      tar(tarfile= file, files = c(files_uni,files_multi,files_corr,files_lme,files_miss,files_IS,files_DATnorm))
    },
    contentType = "application/tar"
  )

  ######export multi analyze######
  ex_multi <- function(){
    if(!is.null(metboshow$metbo_multi)){
      generate_report(metboshow$keepValueM,datsummary1=metboshow$metbo_multi,reportfile="multivariate_report")
      write.csv(metboshow$mul_table,"Multivariate_analysis_output_table.csv")
      dt_plot = list();
      if((input$multiM == "pca" || input$multiM == "pls") && ncol(metboshow$metbo_multi$score_val) > 1 && input$MULScorePC2 != "NA"){#pca or pls
        dt_plot[['scpl']] = multiv_scoreplot(metboshow$keepValueM,cbind(PCx=metboshow$metbo_multi$score_val[,input$MULScorePC1],PCY=metboshow$metbo_multi$score_val[,input$MULScorePC2]),
                                             pcx=metboshow$metbo_multi$model_summary$R2X[metboshow$MULScoreCol1],pcy=metboshow$metbo_multi$model_summary$R2X[metboshow$MULScoreCol2], plot_title="Score plot")
        dt_plot[['loadingpl']] = multiv_loadingplot(cbind(PCX=metboshow$metbo_multi$loading_val[,input$MULScorePC1],PCY=metboshow$metbo_multi$loading_val[,input$MULScorePC2]), plot_title="Loading plot")
        if(input$multiM == "pls"){#vip pls
          filtvp = sort(metboshow$metbo_multi$vip_val, decreasing = TRUE)
          dt_plot[['vippl']] = multiv_vipplot(head(filtvp,n=100))
        }
        filtld = metboshow$metbo_multi$loading_val[order(abs(metboshow$metbo_multi$loading_val[,input$MULselectPC]),decreasing = TRUE),]
        dt_plot[['loadingpc']] = multiv_loadingplot_bypc(head(filtld,n=100), pc = input$MULselectPC, plot_title="Loadings by PC")
      }else if(input$multiM == "opls"){#opls
        if(!is.na(metboshow$metbo_multi$score_val)[,1] && !is.na(metboshow$metbo_multi$oscore_val)[,1]){
          dt_plot[['scpl']] = multiv_scoreplot(metboshow$keepValueM,score_data=metboshow$metbo_multi$score_val,
                                               pcx=metboshow$metbo_multi$model_summary$R2X[1],pcy=metboshow$metbo_multi$model_summary$R2X[2],
                                               oscore_data=metboshow$metbo_multi$oscore_val, plot_title="Score plot")
          dt_plot[['loadingpl']] = multiv_loadingplot(metboshow$metbo_multi$loading_val,oloading_data=metboshow$metbo_multi$oloading_val,plot_title="Loading plot")
          filtvp = sort(metboshow$metbo_multi$vip_val, decreasing = TRUE)
          dt_plot[['vippl']] = multiv_vipplot(head(filtvp,n=100))
          filtld = metboshow$metbo_multi$loading_val[order(abs(metboshow$metbo_multi$loading_val[,1]),decreasing = TRUE),,drop=FALSE]
          dt_plot[['loadingpc']] = multiv_loadingplot_bypc(head(filtld,n=100), pc = 1, plot_title="Loadings of PC1")
        }
      }else{#no model
        return(NULL)
      }
      ggsave("Multivariate_analysis_output_figure.pdf",
             marrangeGrob(grobs = dt_plot, nrow=1, ncol=1, as.table=TRUE), dpi = 600, units = "in", width = 11, height = 7, device = "pdf")
      files <- c("multivariate_report.pdf","Multivariate_analysis_output_table.csv","Multivariate_analysis_output_figure.pdf")
      return(files)
    }else{
      return(NULL)}
  }

  output$export_multi <- downloadHandler(
    filename = function(){
      paste("Metabox_report","tar",sep=".")
    },
    content = function(file){
      shinybusy::show_modal_spinner(color = "#3CB371", text = "Exporting results ...")
      if(metboshow$full == TRUE){
        files_miss <- ex_mis()
        files_IS <- ex_ISnorm()
        files_DATnorm <- ex_DATnorm()
      }else{
        files_miss <- NULL
        files_IS <- NULL
        files_DATnorm <- NULL
      }

      if(metboshow$S_missing_M1 == TRUE){
        files_miss <- ex_mis()
      }else{
        files_miss <- NULL
      }
      files_uni <- ex_uni()
      files_multi <- ex_multi()
      files_corr <- ex_corr()
      files_lme <- ex_lme()
      shinybusy::remove_modal_spinner()
      tar(tarfile= file, files = c(files_uni,files_multi,files_corr,files_lme,files_miss,files_IS,files_DATnorm))
    },
    contentType = "application/tar"
  )

  ######export corr analyze######
  ex_corr <- function(){
    if(!is.null(metboshow$metbo_corr)){
      generate_report(metboshow$keepValueM,datsummary1=metboshow$metbo_corr,reportfile="correlation_report")
      ggsave("corrplot_heatmap.pdf",corrplot_heatmap(metboshow$metbo_corr$corr_data, plot_title="Correlation heatmap"), dpi = 600, units = "in", device = "pdf", width = 10, height = 10)
      write.csv(data.frame(metboshow$metbo_corr$corr_data),"Result_table_correlation.csv")
      files <- c("correlation_report.pdf","corrplot_heatmap.pdf","Result_table_correlation.csv")
      return(files)
    }else{
      return(NULL)}
  }

  output$export_corr <- downloadHandler(
    filename = function(){
      paste("Metabox_report","tar",sep=".")
    },
    content = function(file){
      shinybusy::show_modal_spinner(color = "#3CB371", text = "Exporting results ...")
      if(metboshow$full == TRUE){
        files_miss <- ex_mis()
        files_IS <- ex_ISnorm()
        files_DATnorm <- ex_DATnorm()
      }else{
        files_miss <- NULL
        files_IS <- NULL
        files_DATnorm <- NULL
      }

      if(metboshow$S_missing_M1 == TRUE){
        files_miss <- ex_mis()
      }else{
        files_miss <- NULL
      }

      files_uni <- ex_uni()
      files_multi <- ex_multi()
      files_corr <- ex_corr()
      files_lme <- ex_lme()
      shinybusy::remove_modal_spinner()
      tar(tarfile= file, files = c(files_uni,files_multi,files_corr,files_lme,files_miss,files_IS,files_DATnorm))
    },
    contentType = "application/tar"
  )

  ######export lme analyze######
  ex_lme <- function(){
    if(!is.null(metboshow$metbo_lme)){
      generate_report(metboshow$keepValueM,datsummary1=metboshow$metbo_lme,reportfile="lme_report")
      write.csv(data.frame(metboshow$metbo_lme$testTable),"Result_table_lme.csv")
      files <- c("lme_report.pdf","Result_table_lme.csv")
      return(files)
    }else{
      return(NULL)}
  }

  output$export_lme <- downloadHandler(
    filename = function(){
      paste("Metabox_report","tar",sep=".")
    },
    content = function(file){
      shinybusy::show_modal_spinner(color = "#3CB371", text = "Exporting results ...")
      if(metboshow$full == TRUE){
        files_miss <- ex_mis()
        files_IS <- ex_ISnorm()
        files_DATnorm <- ex_DATnorm()
      }else{
        files_miss <- NULL
        files_IS <- NULL
        files_DATnorm <- NULL
      }

      if(metboshow$S_missing_M1 == TRUE){
        files_miss <- ex_mis()
      }else{
        files_miss <- NULL
      }

      files_uni <- ex_uni()
      files_multi <- ex_multi()
      files_corr <- ex_corr()
      files_lme <- ex_lme()
      shinybusy::remove_modal_spinner()
      tar(tarfile= file, files = c(files_uni,files_multi,files_corr,files_lme,files_miss,files_IS,files_DATnorm))
    },
    contentType = "application/tar"
  )

  ######export muvr analyze######
  ex_muvr <- function(){
    if(!is.null(metboshow$metbo_ml)){
      generate_report(metboshow$keepValueM,datsummary1=metboshow$metbo_ml,reportfile="biomarker_report")
      write.csv(MUVR_getvip(metboshow$metbo_ml,model="min"),"Result_table_biomarker_min_model.csv")
      write.csv(MUVR_getvip(metboshow$metbo_ml,model="mid"),"Result_table_biomarker_mid_model.csv")
      write.csv(MUVR_getvip(metboshow$metbo_ml,model="max"),"Result_table_biomarker_max_model.csv")
      dt_plot = list();
      dt_plot[['plotval']] = MUVR_plotval(metboshow$metbo_ml,plot_title="Performnce evaluation plot")
      dt_plot[['plotvip_min']] = MUVR_plotvip(metboshow$metbo_ml,model="min",plot_title="VIP plot of min-optimal model")
      dt_plot[['plotvip_mid']] = MUVR_plotvip(metboshow$metbo_ml,model="mid",plot_title="VIP plot of mid-optimal model")
      dt_plot[['plotvip_max']] = MUVR_plotvip(metboshow$metbo_ml,model="max",plot_title="VIP plot of max-optimal model")
      ggsave("Biomarker_analysis_output_figure.pdf",
             marrangeGrob(grobs = dt_plot, nrow=1, ncol=1, as.table=TRUE), dpi = 600, units = "in", width = 11, height = 7, device = "pdf")
      files <- c("biomarker_report.pdf","Result_table_biomarker_min_model.csv","Result_table_biomarker_mid_model.csv","Result_table_biomarker_max_model.csv","Biomarker_analysis_output_figure.pdf")
      return(files)
    }else{
      return(NULL)}
  }
  output$export_muvr <- downloadHandler(
    filename = function(){
      paste("Metabox_report","tar",sep=".")
    },
    content = function(file){
      shinybusy::show_modal_spinner(color = "#3CB371", text = "Exporting results ...")
      if(metboshow$full == TRUE){
        files_miss <- ex_mis()
        files_IS <- ex_ISnorm()
        files_DATnorm <- ex_DATnorm()
      }else{
        files_miss <- NULL
        files_IS <- NULL
        files_DATnorm <- NULL
      }

      # if(metboshow$B_missing == TRUE){
      #   files_miss <- ex_mis()
      # }else{
      #   files_miss <- NULL
      # }

      files_muvr <- ex_muvr()
      shinybusy::remove_modal_spinner()
      tar(tarfile= file, files = c(files_muvr,files_miss,files_IS,files_DATnorm))
    },
    contentType = "application/tar"
  )

  #######export MBPLS-DA#####
  ex_mbpls <- function(){
    if(!is.null(metboshow$metbo_mbpl)){
      generate_report(metboshow$metbo_mbpl,datsummary1=metboshow$metbo_mbpl,reportfile="integrative_report")
      dt_plot = list();
      dt_plot[['plotscree']] = mbplsda_screeplot(metboshow$metbo_mbpl$result$base_model$eig,plot_title="Scree plot")
      if(metboshow$MBPLtest){
        dt_plot[['plottest']] = mbplsda_plottestdim(metboshow$metbo_mbpl$result$res_optimal)
        dt_plot[['plotpermu']] = mbplsda_plotpermut(metboshow$metbo_mbpl$result$res_permut,plot_title="Permutation testing plot")
      }
      dt_plot[['plotbip']] = mbplsda_plotboot_bipc(metboshow$metbo_mbpl$result$res_boot,plot_title="Block importance plot")
      dt_plot[['plotvip']] = mbplsda_plotboot_vipc(metboshow$metbo_mbpl$result$res_boot,propbestvar=input$MBPLselectvipprop/100, plot_title=paste0("VIP plot of top ",input$MBPLselectvipprop," percent variables"))
      dt_plot[['plotloadi1']] = mbplsda_plotboot_loading(metboshow$metbo_mbpl$result$res_boot, pc=as.numeric(input$MBPLselectPC), propbestvar=input$MBPLselectloadprop/100,
                                                         plot_title=paste0("Loading plot of top ",input$MBPLselectloadprop, " percent variables in PC", as.numeric(input$MBPLselectPC)))
      write.csv(metboshow$metbo_mbpl$result$res_boot$vipc,"Result_table_integrative_vip.csv")
      write.csv(metboshow$metbo_mbpl$result$res_boot$bipc,"Result_table_integrative_bip.csv")
      file_loading <- c()
      for(axis in 1:length(metboshow$metbo_mbpl$result$res_boot$faX)){
        data_load <- as.data.frame(metboshow$metbo_mbpl$result$res_boot$faX[axis])
        write.csv(data_load,paste0("Result_table_integrative_loading_",names(metboshow$metbo_mbpl$result$res_boot$faX[axis]),".csv"))
        file_loading <- c(file_loading,paste0("Result_table_integrative_loading_",names(metboshow$metbo_mbpl$result$res_boot$faX[axis]),".csv"))
      }
      ggsave("Integrative_analysis_output_figure.pdf",
             marrangeGrob(grobs = dt_plot, nrow=1, ncol=1, as.table=TRUE), dpi = 450, units = "in", width = 11, height = 7, device = "pdf")
      files <- c("integrative_report.pdf","Integrative_analysis_output_figure.pdf","Result_table_integrative_vip.csv","Result_table_integrative_bip.csv",file_loading)
      return(files)
    }else{
      return(NULL)}
  }

  output$export_mbpls <- downloadHandler(
    filename = function(){
      paste("Metabox_report","tar",sep=".")
    },
    content = function(file){
      shinybusy::show_modal_spinner(color = "#3CB371", text = "Exporting results ...")

      files_mbpls <- ex_mbpls()

      shinybusy::remove_modal_spinner()
      tar(tarfile= file, files = files_mbpls)
    },
    contentType = "application/tar"
  )

  ######export ora######
  ex_opa <- function(){
    if(!is.null(metboshow$metbo_overr)){
      generate_report(datsummary1=metboshow$metbo_overr,reportfile="overrep_report")
      write.csv(data.frame(metboshow$metbo_overr$enrichment),"Result_table_overrepresentation.csv")
      files <- c("overrep_report.pdf","Result_table_overrepresentation.csv")
      return(files)
    }else{
      return(NULL)}
  }

  output$export_opa <- downloadHandler(
    filename = function(){
      paste("Metabox_report","tar",sep=".")
    },
    content = function(file){
      shinybusy::show_modal_spinner(color = "#3CB371", text = "Exporting results ...")

      files_cora <- ex_cora()
      files_enr <- ex_enr()
      files_opa <- ex_opa()
      shinybusy::remove_modal_spinner()
      tar(tarfile= file, files = c(files_opa,files_enr,files_cora))
    },
    contentType = "application/tar"
  )


  ######export enr######
  ex_enr <- function(){
    if(!is.null(metboshow$metbo_enr)){
      generate_report(datsummary1=metboshow$metbo_enr,reportfile="enrichment_report")
      write.csv(data.frame(metboshow$metbo_enr$enrichment),"Result_table_enrichment.csv")
      files <- c("enrichment_report.pdf","Result_table_enrichment.csv")
      return(files)
    }else{
      return(NULL)}
  }

  output$export_enr <- downloadHandler(
    filename = function(){
      paste("Metabox_report","tar",sep=".")
    },
    content = function(file){
      shinybusy::show_modal_spinner(color = "#3CB371", text = "Exporting results ...")

      files_cora <- ex_cora()
      files_enr <- ex_enr()
      files_opa <- ex_opa()
      shinybusy::remove_modal_spinner()
      tar(tarfile= file, files = c(files_opa,files_enr,files_cora))
    },
    contentType = "application/tar"
  )


  ######export cora######
  ex_cora <- function(){
    if(!is.null(metboshow$metbo_cORA)){
      generate_report(datsummary1=metboshow$metbo_cORA,reportfile="comb_overrep_report")
      write.csv(data.frame(metboshow$metbo_cORA$enrichment),"Result_table_integratedORA.csv")
      files <- c("comb_overrep_report.pdf","Result_table_integratedORA.csv")
      return(files)
    }else{
      return(NULL)}
  }

  output$export_cora <- downloadHandler(
    filename = function(){
      paste("Metabox_report","tar",sep=".")
    },
    content = function(file){
      shinybusy::show_modal_spinner(color = "#3CB371", text = "Exporting results ...")
      files_cora <- ex_cora()
      files_enr <- ex_enr()
      files_opa <- ex_opa()
      shinybusy::remove_modal_spinner()
      tar(tarfile= file, files = c(files_opa,files_enr,files_cora))
    },
    contentType = "application/tar"
  )


  ####Clear tmp files####
  onStop(function() {
    #metboshow <- NULL
    unlink(tempdir, recursive = TRUE)
    #stopApp()
  })
}
