ui <- fluidPage(theme = bs_theme(version = 3, bootswatch = "lumen", base_font = font_google("Nunito"),
                                 code_font = font_google("Space Mono")),
                # tags$head(
                #   tags$style(
                #     HTML(".shiny-notification {position:fixed;top: calc(90%);left: calc(1%);}")
                #   )
                # ),
                shinyjs::useShinyjs(),
                shinyjs::extendShinyjs(text = app_jscode, functions = c('disableTab','enableTab')),
                extendShinyjs(text = reset_all_evenData, functions = c("resetSelected")),
                navbarPage(windowTitle="Metabox 2.0",img(src="metabox2.png",width=105), id="tabs",

                           ###################################################################
                           ############################HOME PAGE##############################
                           ###################################################################
                           tabPanel("Home",value = "home",
                                    tags$div(HTML("<h2><center>WELCOME TO <b>METABOX</b> 2.0 </center></h2>"), class = "header"),
                                    tags$div(HTML("<h5><center>A toolbox for thorough metabolomic data analysis, integration and interpretation</center></h5>"),class = "header"),
                                    br(),br(),
                                    HTML('<center><img src="first_pipeline.png" width="1150" height = "450" ></center>'),
                                    br(), br(), hr(),
                                    fluidRow(
                                      column(1),
                                      column(2,
                                             actionButton("start", HTML("Start<br>Data Processing<br>(1)"),
                                                          style="margin:0; padding:0; color: #464546; background-color: #87C0CA; width: 150px; height: 65px; font-size: 12px; font-family:  Nunito; font-weight: bold; text-align: center;"
                                             )),
                                      column(2,
                                             actionButton("start2", HTML("Start<br>Statistical Analysis<br>(2)"),
                                                          style="margin:0; padding:0; color: #464546; background-color: #96B8CF; width: 150px; height: 65px; font-size: 12px; font-family:  Nunito; font-weight: bold; text-align: center;"
                                             )),
                                      column(2,
                                             actionButton("start3", HTML("Start<br>Biomaker Analysis<br>(3)"),
                                                          style="margin:0; padding:0; color: #464546; background-color: #8DA7CF; width: 150px; height: 65px; font-size: 12px; font-family:  Nunito; font-weight: bold; text-align: center;"
                                             )),

                                      column(2,
                                             actionButton("start5", HTML("Start<br>Integrative Analysis<br>(4)"),
                                                          style="margin:0; padding:0; color: #464546; background-color: #A6C47D; width: 150px; height: 65px; font-size: 12px; font-family:  Nunito; font-weight: bold; text-align: center;"
                                             )),
                                      column(2,
                                             actionButton("start4", HTML("Start<br>Data Interpretation<br>(5)"),
                                                          style="margin:0; padding:0; color: #464546; background-color: #DDE892; width: 150px; height: 65px; font-size: 12px; font-family:  Nunito; font-weight: bold; text-align: center;"
                                             ))
                                    ),br()
                           ),
                           navbarMenu("More",
                                      ###################################################################
                                      ############################Data PAGE############################
                                      ###################################################################
                                      tabPanel("Example Data",value = "exampleDT",
                                               source("data_page.R", local=TRUE)[1]
                                      ),

                                      ###################################################################
                                      ############################Help PAGE############################
                                      ###################################################################
                                      tabPanel("About",value = "helpPage",
                                               source("about_page.R", local=TRUE)[1]
                                      )),

                           ###################################################################
                           ############################Upload PAGE##############################
                           ###################################################################
                           tabPanel("", value = "upload",
                                    fluidRow(width=12,
                                             sidebarLayout(
                                               column(width=3,
                                                      sidebarPanel(width=12,style = "overflow-y:scroll; max-height: 900px; position:relative;",
                                                                   h4("Start upload your data"),
                                                                   h6("*'Sample_ID' and 'ClassCol can't be used as column names."),
                                                                   br(),
                                                                   fileInput("INfile",h5("Upload a csv file"),accept=".csv"),
                                                                   br(),
                                                                   selectInput("idcol","Sample ID column",choices=NULL, selected=NULL),
                                                                   br(),
                                                                   selectInput("classS","Class/response column (Y)",choices=NULL, selected=NULL),
                                                                   br(),
                                                                   selectInput("firstV","1st variable/feature column (X)",choices=NULL, selected=NULL),
                                                                   br(),
                                                                   actionButton("next1", "Next", align = "right")
                                                      )),
                                               column(width=9,
                                                      mainPanel(width=12,
                                                                fluidRow(
                                                                  column(12,id = "UPsummary",
                                                                         tabsetPanel(id = "UPsummaryT",
                                                                                      tabPanel("Summary", value = "UPsummaryTab",
                                                                                               div(verbatimTextOutput("txtbox"),style = "overflow-y:scroll;  max-height: 350px;")
                                                                                      )),hr(style="border-top: 1px dashed #B5C7DA;")
                                                                  ),
                                                                  column(12,id = "UPdata",
                                                                         tabsetPanel(id = "UPdataT",
                                                                                      tabPanel("Data table",value = "UPdataTab",br(),
                                                                                               div(DT::dataTableOutput('INtable'), style = "overflow-y:scroll; overflow-x:scroll;max-height: 600px; position:relative; max-height: 600px;"),br()
                                                                                      ))
                                                                  )
                                                                ))
                                               ))
                                    )),

                           ###################################################################
                           ############################DP PAGE##############################
                           ###################################################################
                           tabPanel("",value = "DPtabs",
                                    tabsetPanel( id = "DPtabset",
                                                 #####Impute missing value#########################################
                                                 tabPanel("Missing value imputation", value = "missing",
                                                          fluidRow(width=12,
                                                                   sidebarLayout(
                                                                     column(width=3,
                                                                            sidebarPanel(width=12,style = "overflow-y:scroll; max-height: 900px; position:relative;",
                                                                                         h4("Missing value imputation"),
                                                                                         radioButtons("reall",
                                                                                                      h5("Remove all variables with missing values?"),
                                                                                                      choices = list("Remove all" = 'TRUE',
                                                                                                                     "Not remove all" = 'FALSE'
                                                                                                      ),
                                                                                                      selected = 'FALSE'),
                                                                                         br(),
                                                                                         sliderInput("remPercent", h5("Remove variables with missing values > (%)"),
                                                                                                     min = 1, max = 99, value = 20),
                                                                                         br(),
                                                                                         radioButtons("missCheck", h5("Choose imputation method"),
                                                                                                      choices = list("Replace by zero" = "zero",
                                                                                                                     "Half-minimum" = "halfmin",
                                                                                                                     "Minimum" = "min",
                                                                                                                     "Mean" = "mean",
                                                                                                                     "Median" = "median",
                                                                                                                     "KNN" = "knn",
                                                                                                                     "BPCA" = "bpca",
                                                                                                                     "PPCA" = "ppca",
                                                                                                                     "SVD" = "svd",
                                                                                                                     "RF" = "rf"
                                                                                                      ),selected = 'mean'),
                                                                                         br(),
                                                                                         actionButton("run1", "Run", class = "btn btn-primary"),
                                                                                         actionButton("back1", "Back"),
                                                                                         actionButton("next2", "Next"),
                                                                                         actionButton("skip1", "Skip This Step"),br(),br(),
                                                                                         downloadButton("export_mis", "Export results", align = "center", class = "btn btn-success")
                                                                            )),
                                                                     column(width=9,
                                                                            mainPanel(width=12,br(),
                                                                                      fluidRow(
                                                                                        column(12,id = "MIsummary",
                                                                                               tabsetPanel(id = "MIsummaryT",
                                                                                                           tabPanel("Summary", value = "MIsummaryTab",
                                                                                                                    div(verbatimTextOutput("txtbox.missing"),style = "overflow-y:scroll; max-height: 350px;")
                                                                                                           )
                                                                                               ),hr(style="border-top: 1px dashed #B5C7DA;")),
                                                                                        column(12, id = "MIoutput",
                                                                                               tabsetPanel(id = "MIoutputT",
                                                                                                            tabPanel("Output table", value = "MIoutputTab",br(),
                                                                                                                     div(DT::dataTableOutput('INtable3'), style = "overflow-y:scroll; overflow-x:scroll;max-height: 600px; position:relative; max-height: 600px;"),br()
                                                                                                            )
                                                                                               )

                                                                                        ))
                                                                            )
                                                                     )
                                                                   ))
                                                 ),
                                                 #####IS-QC base norm##############################################
                                                 tabPanel("IS- or QC-based normalization", value = "QCnorm",
                                                          fluidRow(width=12,
                                                                   sidebarLayout(
                                                                     column(3,
                                                                            sidebarPanel(width=12,style = "overflow-y:scroll; max-height: 900px; position:relative; height: 600px;",
                                                                                         h4("IS- or QC-based normalization"),
                                                                                         h6("*serrf and loess methods take long execution time."),
                                                                                         radioButtons("QCnormSelect", h5("Choose normalization method"),
                                                                                                      choices = list("none" = "none",
                                                                                                                     "ccmn" = "ccmn",
                                                                                                                     "nomis" = "nomis",
                                                                                                                     "ruvrand" = "ruvrand",
                                                                                                                     "serrf" = "serrf",
                                                                                                                     "loess" = "loess"
                                                                                                      ),selected = "none"),
                                                                                         selectInput("exdesign","Class/factor column",choices=NULL, selected=NULL),
                                                                                         selectInput("istd","Internal standard column(s)",choices=NULL, selected=NULL,multiple = TRUE),
                                                                                         selectInput("sampleType","sampleType column (require at least 3 QCs)",choices=NULL, selected=NULL),
                                                                                         selectInput("injectionOrder","injectionOrder column",choices=NULL, selected=NULL),
                                                                                         selectInput("batch","batch column",choices=NULL, selected=NULL),
                                                                                         actionButton("run2", "Run", class = "btn btn-primary"),
                                                                                         actionButton("back2", "Back"),
                                                                                         actionButton("next3", "Next"),
                                                                                         actionButton("skip2", "Skip This Step"),br(),br(),
                                                                                         downloadButton("export_normIS", "Export results", align = "center", class = "btn btn-success")
                                                                            )),
                                                                     column(9,
                                                                            mainPanel(width=12,br(),
                                                                                      fluidRow(
                                                                                        column(12, id = "QCsummary",
                                                                                               tabsetPanel(id = "QCsummaryT",
                                                                                                            tabPanel("Summary", value = "QCsummaryTab",
                                                                                                                     div(verbatimTextOutput("txtbox.QCnorm"),style = " overflow-y:scroll; max-height: 350px;")
                                                                                                            )
                                                                                               ),hr(style="border-top: 1px dashed #B5C7DA;")),
                                                                                        column(12, id = "QCplotTab",
                                                                                               tabsetPanel(id = "QCtabset",
                                                                                                           tabPanel("Overview", value = "QCoverview",br(),
                                                                                                                    column(12,
                                                                                                                           h5("PCA score plot of the 1st two PCs provides the overview of the data."),hr(style="border-top: 1px dashed #D3D3D3;"),
                                                                                                                           column(6,
                                                                                                                                  plotlyOutput(outputId = "plotQCoverOri",height = "480px")
                                                                                                                           ),
                                                                                                                           column(6,
                                                                                                                                  div(verbatimTextOutput("txtbox.QCoverOri"),style = " overflow-y:scroll; max-height: 350px;")
                                                                                                                           )
                                                                                                                    ),
                                                                                                                    column(12,br(),br(),
                                                                                                                           column(6,
                                                                                                                                  plotlyOutput(outputId = "plotQCover",height = "480px"),br()
                                                                                                                           ),
                                                                                                                           column(6,
                                                                                                                                  h5("*Drag on the plot to select sample(s), double-click to unselect."),
                                                                                                                                  h5("*To exclude sample(s)/outlier(s) in further analysis, click REMOVE SAMPLE."),
                                                                                                                                  div(verbatimTextOutput("txtbox.QCover"),style = "overflow-y:scroll; max-height: 350px;"),
                                                                                                                                  actionButton("rev_is", "Remove sample",style="float:right;"),br()
                                                                                                                           )
                                                                                                                    ),
                                                                                                           ),
                                                                                                           # tabPanel("Boxplot", value = "QCboxplot",
                                                                                                           #          column(12,
                                                                                                           #                 plotlyOutput(outputId = "plotQCboxOri",height = 900)
                                                                                                           #          ),
                                                                                                           #          column(12,
                                                                                                           #                 plotlyOutput(outputId = "plotQCbox",height = 900)
                                                                                                           #          )),
                                                                                                           tabPanel("RLA plot", value = "QCrlaplot",br(),
                                                                                                                    column(12,
                                                                                                                           h5("RLA plot provides the overview of the replicates within each group."),h5("*Display only the first 100 samples."),
                                                                                                                           hr(style="border-top: 1px dashed #D3D3D3;"),
                                                                                                                           plotOutput(outputId = "plotQCboxOriS",height = 800)
                                                                                                                    ),
                                                                                                                    column(12,
                                                                                                                           plotOutput(outputId = "plotQCboxS",height = 800)
                                                                                                                    )
                                                                                                           ),
                                                                                                           tabPanel("Density", value = "QCdensity",br(),
                                                                                                                    column(12,
                                                                                                                           h5("Density plot provides the distribution of variables and samples."),hr(style="border-top: 1px dashed #D3D3D3;")),
                                                                                                                    column(12,
                                                                                                                           column(6,
                                                                                                                                  plotOutput(outputId = "plotQCdenOri")
                                                                                                                           ),
                                                                                                                           column(6,
                                                                                                                                  plotOutput(outputId = "plotQCden")
                                                                                                                           )),
                                                                                                                    column(12,br(),
                                                                                                                           column(6,
                                                                                                                                  plotOutput(outputId = "plotQCdenOriS")
                                                                                                                           ),
                                                                                                                           column(6,
                                                                                                                                  plotOutput(outputId = "plotQCdenS")
                                                                                                                           ))
                                                                                                           ),
                                                                                                           tabPanel("Output table", value = "QCoutputTab",br(),
                                                                                                                    actionButton("deleteC_qc", "Delete column"),br(),
                                                                                                                    h5("*To exclude variable(s) in further analysis, select column(s) and click DELETE COLUMN."),hr(style="border-top: 1px dashed #D3D3D3;"),
                                                                                                                    div(DT::dataTableOutput('INtable4'), style = " overflow-y:scroll; overflow-x:scroll;max-height: 600px; position:relative; max-height: 600px;"),br()
                                                                                                           )
                                                                                               )
                                                                                        ))



                                                                            )
                                                                     ))
                                                          )),
                                                 #####Data-driven Normalization######################################
                                                 tabPanel("Data-driven normalization", value = "norm",
                                                          fluidRow(width=12,
                                                                   sidebarLayout(
                                                                     column(3,
                                                                            sidebarPanel(width=12,style = "overflow-y:scroll; max-height: 900px; position:relative;",
                                                                                         h4("Sample-based normalization"),
                                                                                         radioButtons("DATAnormSelect", h5("Choose normalization method"),
                                                                                                      choices = list("none" = "none",
                                                                                                                     "contrast" = "contrast",
                                                                                                                     "cubic" = "cubic",
                                                                                                                     "eigenms" = "eigenms",
                                                                                                                     "linear" = "linear",
                                                                                                                     "liwong" = "liwong",
                                                                                                                     "loess" = "loess",
                                                                                                                     "mean" = "mean",
                                                                                                                     "median" = "median",
                                                                                                                     "pqn" = "pqn",
                                                                                                                     "quantile" = "quantile",
                                                                                                                     "sum" = "sum",
                                                                                                                     "vsn" = "vsn"
                                                                                                      ),selected = "none"),
                                                                                         selectInput("ref","Reference row",choices=NULL, selected=NULL,multiple = TRUE),
                                                                                         br(),
                                                                                         h4("Transformation"),
                                                                                         radioButtons("tranM", h5("Choose transformation method"),
                                                                                                      choices = list("none" = "none",
                                                                                                                     "cube root" = "cube",
                                                                                                                     "log10" = "log10",
                                                                                                                     "log2" = "log2",
                                                                                                                     "glog10" = "glog10",
                                                                                                                     "glog2" = "glog2",
                                                                                                                     "square root" = "sqrt"

                                                                                                      ),selected = "none"),
                                                                                         br(),
                                                                                         h4("Scaling"),
                                                                                         radioButtons("scalM", h5("Choose scaling method"),
                                                                                                      choices = list("none" = "none",
                                                                                                                     "auto" = "auto",
                                                                                                                     "level" = "level",
                                                                                                                     "pareto" = "pareto",
                                                                                                                     "power" = "power",
                                                                                                                     "range" = "range",
                                                                                                                     "vast" = "vast"
                                                                                                      ),selected = "none"),
                                                                                         br(),
                                                                                         actionButton("run3", "Run", class = "btn btn-primary"),
                                                                                         actionButton("back3", "Back"),
                                                                                         br(),br(),
                                                                                         downloadButton("export_normDAT", "Export results", align = "center", class = "btn btn-success"),
                                                                                         br(),hr(style="border: 1px dashed #ABB2B9;"),
                                                                                         actionButton("next4", "Statistical Analysis", style="color: #464546; background-color: #96B8CF; font-weight: bold; height: 45px; "),
                                                                                         actionButton("next5", "Biomarker Analysis", style="color: #464546; background-color: #8DA7CF; font-weight: bold; height: 45px; ")
                                                                            )),
                                                                     column(9,
                                                                            mainPanel(width=12,br(),
                                                                                      fluidRow(
                                                                                        column(12, id = "Nsummary",
                                                                                               tabsetPanel(id = "NsummaryT",
                                                                                                            tabPanel("Summary", value = "NsummaryTab",
                                                                                                                     div(verbatimTextOutput("txtbox.DATnorm"),style = " overflow-y:scroll; max-height: 350px;")
                                                                                                            )
                                                                                               ),hr(style="border-top: 1px dashed #B5C7DA;")),
                                                                                        column(12, id = "NplotTab",
                                                                                               tabsetPanel(id = "Ntabset",
                                                                                                            tabPanel("Overview", value = "Noverview", br(),
                                                                                                                     column(12,
                                                                                                                            h5("PCA score plot of the first two PCs provides the overview of the data."),hr(style="border-top: 1px dashed #D3D3D3;"),
                                                                                                                            column(6,
                                                                                                                                   plotlyOutput(outputId = "plotNoverOri",height = "480px", width = "100%")
                                                                                                                            ),
                                                                                                                            column(6,
                                                                                                                                   div(verbatimTextOutput("txtbox.NoverOri"),style = "overflow-y:scroll; max-height: 350px;")
                                                                                                                            )
                                                                                                                     ),
                                                                                                                     column(12,br(),br(),
                                                                                                                            column(6,
                                                                                                                                   plotlyOutput(outputId = "plotNover",height = "480px", width = "100%"),br()
                                                                                                                            ),
                                                                                                                            column(6,
                                                                                                                                   h5("*Drag on the plot to select sample(s), double-click to unselect."),
                                                                                                                                   h5("*To exclude sample(s)/outlier(s) in further analysis, click REMOVE SAMPLE."),
                                                                                                                                   div(verbatimTextOutput("txtbox.Nover"),style = "overflow-y:scroll; max-height: 350px;"),
                                                                                                                                   actionButton("rev_N", "Remove sample",style="float:right;"),br()
                                                                                                                            )
                                                                                                                     )
                                                                                                            ),
                                                                                                            # tabPanel("Boxplot", value = "Nboxplot",
                                                                                                            #          column(6,
                                                                                                            #                 htmltools::div(style = "display:inline-block", plotlyOutput(outputId = "plotNboxOri", width = 250, height = 400))
                                                                                                            #          ),
                                                                                                            #          column(6,style = "max-height:600px; overflow-y:scroll; ",
                                                                                                            #                 plotlyOutput(outputId = "plotNbox",inline = T)
                                                                                                            #          ),br(),
                                                                                                            #          column(6,style = "max-height:600px; overflow-y:scroll; ",
                                                                                                            #                 plotlyOutput(outputId = "plotNboxOriS",inline = T)
                                                                                                            #          ),
                                                                                                            #          column(6,style = "max-height:600px; overflow-y:scroll; ",
                                                                                                            #                 plotlyOutput(outputId = "plotNboxS",inline = T)
                                                                                                            #          )
                                                                                                            # ),
                                                                                                            # tabPanel("Boxplot", value = "Nboxplot",
                                                                                                            #          column(12,
                                                                                                            #                 plotlyOutput(outputId = "plotNboxOri",height = 900)
                                                                                                            #          ),
                                                                                                            #          column(12,
                                                                                                            #                 plotlyOutput(outputId = "plotNbox",height = 900)
                                                                                                            #          )),
                                                                                                            tabPanel("RLA plot", value = "Nrlaplot",br(),
                                                                                                                     column(12,
                                                                                                                            h5("RLA plot provides the overview of the replicates within each group."),h5("*Display only the first 100 samples."),
                                                                                                                            hr(style="border-top: 1px dashed #D3D3D3;"),
                                                                                                                            plotOutput(outputId = "plotNboxOriS",height = 800)
                                                                                                                     ),
                                                                                                                     column(12,
                                                                                                                            plotOutput(outputId = "plotNboxS",height = 800)
                                                                                                                     )
                                                                                                            ),
                                                                                                            tabPanel("Density", value = "Ndensity",br(),
                                                                                                                     column(12,
                                                                                                                            h5("Density plot provides the distribution of variables and samples."),hr(style="border-top: 1px dashed #D3D3D3;")),
                                                                                                                     column(12,
                                                                                                                            column(6,
                                                                                                                                   plotOutput(outputId = "plotNdenOri")
                                                                                                                            ),
                                                                                                                            column(6,
                                                                                                                                   plotOutput(outputId = "plotNden")
                                                                                                                            )),
                                                                                                                     column(12,br(),
                                                                                                                            column(6,
                                                                                                                                   plotOutput(outputId = "plotNdenOriS")
                                                                                                                            ),
                                                                                                                            column(6,
                                                                                                                                   plotOutput(outputId = "plotNdenS")
                                                                                                                            ))
                                                                                                            ),
                                                                                                            tabPanel("Output table", value = "NoutputTab",br(),
                                                                                                                     actionButton("deleteC_norm", "Delete column"),br(),
                                                                                                                     h5("*To exclude variable(s) in further analysis, select column(s) and click DELETE COLUMN."),hr(style="border-top: 1px dashed #D3D3D3;"),
                                                                                                                     div(DT::dataTableOutput('INtable5'), style = "overflow-x:scroll;max-height: 600px; position:relative; max-height: 600px;"),br()
                                                                                                            )
                                                                                               )
                                                                                        ))
                                                                            )
                                                                     ))
                                                          ))
                                    )),

                           ###################################################################
                           ############################stat PAGE#############################
                           ###################################################################

                           tabPanel("", value = "statTabs",
                                    tabsetPanel(id = "statTabset",

                                                #####Univariate analysis#############################
                                                tabPanel("Univariate analysis", value = "uni",
                                                         fluidRow(width=12,
                                                                  sidebarLayout(
                                                                    column(3,
                                                                           sidebarPanel(width=12,style = "overflow-y:scroll; max-height: 900px; position:relative; ",
                                                                                        h4("Univariate analysis"),
                                                                                        radioButtons("var.equal",
                                                                                                     h5("Equal variance?"),
                                                                                                     choices = list("yes" = 'TRUE',
                                                                                                                    "no" = 'FALSE'
                                                                                                     ),
                                                                                                     selected = 'TRUE'),
                                                                                        radioButtons("ispara",
                                                                                                     h5("Parametric test?"),
                                                                                                     choices = list("yes" = 'TRUE',
                                                                                                                    "no" = 'FALSE'
                                                                                                     ),
                                                                                                     selected = 'TRUE'),
                                                                                        selectInput("factor2Col","2nd factor column (for 2-way ANOVA)",choices=NULL, selected=NULL),
                                                                                        radioButtons("doposthoc",
                                                                                                     h5("Post hoc test? (for 1-way and 2-way ANOVA)"),
                                                                                                     choices = list("yes" = 'TRUE',
                                                                                                                    "no" = 'FALSE'
                                                                                                     ),
                                                                                                     selected = 'FALSE'),
                                                                                        h6("*Post hoc test might take a long execution time."),
                                                                                        br(),
                                                                                        actionButton("run6", "Run", class = "btn btn-primary"),
                                                                                        actionButton("back6", "Back"),
                                                                                        br(),br(),
                                                                                        downloadButton("export_uni", "Export results", align = "center", class = "btn btn-success"),
                                                                                        br(),hr(style="border: 1px dashed #ABB2B9;"),
                                                                                        actionButton("nextS_B1", "Biomarker Analysis", style="color: #464546; background-color: #8DA7CF; font-weight: bold; height: 45px; "),
                                                                                        actionButton("nextS_M1", "Impute missing value",style="color: #464546; background-color: #FFDCD1; font-weight: bold; height: 45px; ")
                                                                           )),
                                                                    column(9,
                                                                           mainPanel(width=12,br(),
                                                                                     fluidRow(
                                                                                       column(12, id = "UNIsummary",
                                                                                              tabsetPanel(id = "UNIsummaryT",
                                                                                                           tabPanel("Summary", value = "UNIsummaryTab",
                                                                                                                    div(verbatimTextOutput("txtbox.uni"),style = "overflow-y:scroll;  max-height: 350px;")
                                                                                                           )
                                                                                              ),hr(style="border-top: 1px dashed #B5C7DA;")),
                                                                                       column(12, id = "UNIplotTab",
                                                                                              tabsetPanel(id = "UNItabset",
                                                                                                           tabPanel("Overview", value = "UNIoverview",br(),
                                                                                                                    column(12,
                                                                                                                           h5("Statistical significance plot provides the -log10(adjusted p-value) of each variable.
                                                                                                                               Dashed line represents statistical significance cutoff (adjusted p-value < 0.05)."),
                                                                                                                           h5("*Display only the top 100 variables, sorte by adjusted p-values. Click on a dot to toggle its boxplot."),
                                                                                                                           hr(style="border-top: 1px dashed #D3D3D3;"),
                                                                                                                           plotlyOutput(outputId = "plotPvalUNI",height = "480px")
                                                                                                                    ),
                                                                                                                    column(12,
                                                                                                                           plotlyOutput(outputId = "plotBoxUNI",height = "480px")
                                                                                                                    )
                                                                                                           ),
                                                                                                           tabPanel("Output table", value = "UNIoutputTab",
                                                                                                                    h5("*For 2-level factor/class, fold change (fc) = level1/level2, for more than 2 levels, fc = each_level/all_mean. Statistical significance for post hoc test is adjusted p-value < 0.05."),
                                                                                                                    div(DT::dataTableOutput('INtableUNI'), style = "overflow-y:scroll; overflow-x:scroll;max-height: 600px; position:relative; max-height: 600px;"),
                                                                                                                    br()
                                                                                                           )

                                                                                              )
                                                                                       )
                                                                                     ))
                                                                    ))
                                                         )),

                                                #####Multivariate analysis#############################
                                                tabPanel("Multivariate analysis", value = "multi",
                                                         fluidRow(width=12,
                                                                  sidebarLayout(
                                                                    column(width=3,
                                                                           sidebarPanel(width=12,style = "overflow-y:scroll; max-height: 900px; position:relative;",
                                                                                        h4("Multivariate analysis"),
                                                                                        radioButtons("multiM", h5("Choose analysis method"),
                                                                                                     choices = list(
                                                                                                       "PCA" = "pca",
                                                                                                       "PLS" = "pls",
                                                                                                       "OPLS" = "opls"
                                                                                                     ),selected = "pca"),
                                                                                        radioButtons("Mscale", h5("Choose scaling method, if data wasn't processed"),
                                                                                                     choices = list(
                                                                                                       "None" = "center",
                                                                                                       "Pareto" = "pareto",
                                                                                                       "Auto" = "standard"
                                                                                                     ),selected = "center"),
                                                                                        actionButton("run7", "Run", class = "btn btn-primary"),
                                                                                        actionButton("back7", "Back"),
                                                                                        br(),br(),
                                                                                        downloadButton("export_multi", "Export results", align = "center", class = "btn btn-success"),
                                                                                        br(),hr(style="border: 1px dashed #ABB2B9;"),
                                                                                        actionButton("nextS_B2", "Biomarker Analysis", style="color: #464546; background-color: #8DA7CF; font-weight: bold; height: 45px; "),
                                                                                        actionButton("nextS_M2", "Impute missing value",style="color: #464546; background-color: #FFDCD1; font-weight: bold; height: 45px; ")
                                                                           )),
                                                                    column(width=9,
                                                                           mainPanel(width=12,br(),
                                                                                     fluidRow(
                                                                                       column(12, id = "MULsummary",
                                                                                              tabsetPanel(id = "MULsummaryT",
                                                                                                          tabPanel("Summary",value = "MULsummaryTab",
                                                                                                                   div(verbatimTextOutput("txtbox.multi"),style = "overflow-y:scroll; max-height: 350px;")                                                                                                                     )
                                                                                              ),hr(style="border-top: 1px dashed #B5C7DA;")),
                                                                                       column(12, id = "MULplotTab",
                                                                                              tabsetPanel(id = "MULtabset",
                                                                                                           tabPanel("Overview", value = "MULoverview",br(),
                                                                                                                    column(12,
                                                                                                                           h5("Score and loading plots of 2 principal components (PC)."),
                                                                                                                           hr(style="border-top: 1px dashed #D3D3D3;"),
                                                                                                                           column(3,selectInput("MULScorePC1","Select PC to display on x-axis (PCX)",choices=NULL, selected=NULL)),
                                                                                                                           column(3,selectInput("MULScorePC2","Select PC to display on y-axis (PCY)",choices=NULL, selected=NULL)),br()
                                                                                                                    ),
                                                                                                                    column(6,
                                                                                                                           plotlyOutput(outputId = "plotScoreMUL",height = "480px"),br()
                                                                                                                    ),
                                                                                                                    column(6,
                                                                                                                           plotlyOutput(outputId = "plotLoadMUL",height = "480px"),br()
                                                                                                                    )
                                                                                                           ),
                                                                                                           tabPanel("Loading", value = "MULloading",br(),
                                                                                                                    h5("Loading plot of a principal component (PC)."),
                                                                                                                    h5("*Display only the top 100 variables, sorted by loading values."),
                                                                                                                    hr(style="border-top: 1px dashed #D3D3D3;"),
                                                                                                                    selectInput("MULselectPC","Select PC to display",choices=NULL, selected=NULL),
                                                                                                                    column(12,
                                                                                                                           plotlyOutput(outputId = "plotPCloadMUL",height = "480px")
                                                                                                                    )
                                                                                                           ),
                                                                                                           tabPanel("VIP", value = "MULvip",br(),
                                                                                                                    column(12,
                                                                                                                           h5("Variable importance in projection (VIP) plot."),
                                                                                                                           h5("*Display only the top 100 variables, sorted by VIP values."),
                                                                                                                           hr(style="border-top: 1px dashed #D3D3D3;"),
                                                                                                                           plotlyOutput(outputId = "plotvipMUL",height = "480px")
                                                                                                                    )
                                                                                                           ),
                                                                                                           tabPanel("Output table", value = "MULoutputTab",br(),
                                                                                                                    h5("Output table contains loading values of variables in each PC and VIP values for PLS and OPLS."),
                                                                                                                    hr(style="border-top: 1px dashed #D3D3D3;"),
                                                                                                                    div(DT::dataTableOutput('INtableMUL'), style = "overflow-y:scroll; overflow-x:scroll;max-height: 600px; position:relative; max-height: 600px;"),br()
                                                                                                           )

                                                                                              )
                                                                                       )

                                                                                     )
                                                                           )
                                                                    ))
                                                         )),

                                                #####Correlation analysis#############################
                                                tabPanel("Correlation analysis", value = "corr",
                                                         fluidRow(width=12,
                                                                  sidebarLayout(
                                                                    column(width=3,
                                                                           sidebarPanel(width=12,style = "overflow-y:scroll; max-height: 900px; position:relative;",
                                                                                        h4("Correlation analysis"),
                                                                                        radioButtons("corrM", h5("Choose correlation method"),
                                                                                                     choices = list(
                                                                                                       "pearson" = "pearson",
                                                                                                       "spearman" = "spearman"
                                                                                                     ),selected = "pearson"),
                                                                                        br(),
                                                                                        # fileInput("INfileCorr",h5("Upload your csv file")),
                                                                                        # br(),
                                                                                        actionButton("run8", "Run", class = "btn btn-primary"),
                                                                                        actionButton("back8", "Back"),
                                                                                        br(),br(),
                                                                                        downloadButton("export_corr", "Export results", align = "center", class = "btn btn-success"),
                                                                                        br(),hr(style="border: 1px dashed #ABB2B9;"),
                                                                                        actionButton("nextS_B3", "Biomarker Analysis", style="color: #464546; background-color: #8DA7CF; font-weight: bold; height: 45px; "),
                                                                                        actionButton("nextS_M3", "Impute missing value",style="color: #464546; background-color: #FFDCD1; font-weight: bold; height: 45px; ")
                                                                           )),
                                                                    column(9,
                                                                           mainPanel(width=12,br(),
                                                                                     fluidRow(
                                                                                       column(12, id = "CORsummary",
                                                                                              tabsetPanel(id = "CORsummaryT",
                                                                                                          tabPanel("Summary",value = "CORsummaryTab",
                                                                                                                   div(verbatimTextOutput("txtbox.corr"),style = "overflow-y:scroll; max-height: 350px;")
                                                                                                          )
                                                                                              ),hr(style="border-top: 1px dashed #B5C7DA;")),
                                                                                       column(12, id = "CORplotTab",
                                                                                              tabsetPanel(id = "CORtabset",
                                                                                                           tabPanel("Overview", value = "CORoverview",br(),
                                                                                                                    column(12,
                                                                                                                           h5("*Drag on the plot to zoom, double-click to reset."),
                                                                                                                           h5("*Plotting will take awhile for a lot of variables."),br(),
                                                                                                                           plotlyOutput(outputId = "plotHeatCOR", height = "680px")
                                                                                                                    )

                                                                                                           ),
                                                                                                           tabPanel("Output table", value = "CORoutputTab",br(),
                                                                                                                    h5("Output table contains correlation coefficient, p-value and adjusted p-value between two variables."),
                                                                                                                    hr(style="border-top: 1px dashed #D3D3D3;"),
                                                                                                                    div(DT::dataTableOutput('INtableCOR'), style = "overflow-y:scroll; overflow-x:scroll;max-height: 600px; position:relative; max-height: 600px;"),br()
                                                                                                           )
                                                                                              )
                                                                                       )

                                                                                     )
                                                                           )
                                                                    ))
                                                         )),

                                                #####Linear mixed-effects modeling#############################
                                                tabPanel("Linear mixed-effects modeling", value = "linear",
                                                         fluidRow(width=12,
                                                                  sidebarLayout(
                                                                    column(width=3,
                                                                           sidebarPanel(width=12,style = "overflow-y:scroll; max-height: 900px; position:relative; ",
                                                                                        h4("Linear mixed-effects modeling"),
                                                                                        tags$div(HTML("<p>*Implemented from lmm2met package,
                                                                                        see <a href='https://doi.org/10.1016/j.csbj.2019.04.009', target='_blank'>reference</a></p>")),
                                                                                        selectInput("LMEfix",h5("Choose Fixed-effect column(s)"),choices=NULL, selected=NULL,multiple = TRUE, width = 300),
                                                                                        br(),
                                                                                        textInput("LMErandom",h5("Fill in Random-effect formular e.g. (1|SampID)"), width = 300),
                                                                                        br(),
                                                                                        actionButton("run9", "Run", class = "btn btn-primary"),
                                                                                        actionButton("back9", "Back"),
                                                                                        br(),br(),
                                                                                        downloadButton("export_lme", "Export results", align = "center", class = "btn btn-success"),
                                                                                        br(),hr(style="border: 1px dashed #ABB2B9;"),
                                                                                        actionButton("nextS_B4", "Biomarker Analysis", style="color: #464546; background-color: #8DA7CF; font-weight: bold; height: 45px; "),
                                                                                        actionButton("nextS_M4", "Impute missing value",style="color: #464546; background-color: #FFDCD1; font-weight: bold; height: 45px; ")
                                                                           )),
                                                                    column(9,
                                                                           mainPanel(width=12,br(),
                                                                                     fluidRow(
                                                                                       column(12, id = "LMEsummary",
                                                                                              tabsetPanel(id = "LMEsummaryT",
                                                                                                          tabPanel("Summary",value = "LMEsummaryTab",
                                                                                                                   div(verbatimTextOutput("txtbox.lme"),style = "overflow-y:scroll; max-height: 350px;")
                                                                                                          )
                                                                                              ),hr(style="border-top: 1px dashed #B5C7DA;")),
                                                                                       column(12, id = "LMEplotTab",
                                                                                              tabsetPanel(id = "LMEtabset",
                                                                                                           tabPanel("Output table", value = "LMEoutputTab",br(),
                                                                                                                    h5("Output table contains Intercept, coefficients and chi-square test results for significant fixed effects."),
                                                                                                                    hr(style="border-top: 1px dashed #D3D3D3;"),
                                                                                                                    div(DT::dataTableOutput('INtableLME'), style = "overflow-y:scroll; overflow-x:scroll;max-height: 600px; position:relative; max-height: 600px;"),br()
                                                                                                           )
                                                                                              )
                                                                                       )
                                                                                     )

                                                                           ))
                                                                  )
                                                         ))
                                    )),


                           ###################################################################
                           ############################Biomarker PAGE#########################
                           ###################################################################
                           tabPanel("", value = "bio",
                                    fluidRow(width=12,
                                             sidebarLayout(
                                               column(width=3,
                                                      sidebarPanel(width=12,style = "overflow-y:scroll; max-height: 900px; position:relative; height: 750px;",
                                                                   h4("Biomarker analysis"),
                                                                   tags$div(HTML("<p>*Implemented from MUVR package, see <a href='https://doi.org/10.1093/bioinformatics/bty710', target='_blank'>reference</a></p>")),
                                                                   h6("*It might take long execution time depended on data set size and the number of repetitions."),
                                                                   radioButtons("muvr", h5("Choose ML analysis method"),
                                                                                choices = list("PLS" = "PLS",
                                                                                               "RF" = "RF"
                                                                                ),selected = "PLS"),
                                                                   numericInput("nRep", "Number of repetitions", 3, min = 1),
                                                                   numericInput("varRatio", "Ratio of variables included in inner loop iterations (max = 1)", 0.75, min = 0, max = 1, step=0.1),
                                                                   numericInput("partitionSize", "Ratio of training set (max = 1)", 0.7, min = 0, max = 1, step=0.1),
                                                                   radioButtons("da", h5("Choose analysis task"),
                                                                                choices = list("Classification" = "tr",
                                                                                               "Regression" = "fa"
                                                                                ),selected = "tr"),
                                                                   radioButtons("fitn", h5("Choose performance evaluation method"),
                                                                                choices = list("AUROC" = "AUROC",
                                                                                               "MISS" = "MISS"
                                                                                ),selected = "AUROC"),
                                                                   actionButton("run5", "Run", class = "btn btn-primary"),
                                                                   actionButton("back5", "Back"),
                                                                   br(),br(),
                                                                   downloadButton("export_muvr", "Export results", align = "center", class = "btn btn-success"),
                                                                   br(),hr(style="border: 1px dashed #ABB2B9;"),
                                                                   actionButton("nextB_S", "Statistical Analysis", style="color: #464546; background-color: #96B8CF; font-weight: bold; height: 45px; "),
                                                                   actionButton("nextB_M", "Impute missing value",style="color: #464546; background-color: #FFDCD1; font-weight: bold; height: 45px; ")

                                                      )),
                                               column(9,
                                                      mainPanel(width=12,br(),
                                                                fluidRow(
                                                                  column(12, id = "MUVsummary",
                                                                         tabsetPanel(id = "MUVsummaryT",
                                                                                      tabPanel("Summary",value = "MUVsummaryTab",
                                                                                               div(verbatimTextOutput("txtbox.muv"),style = "overflow-y:scroll;  max-height: 350px;")
                                                                                      )
                                                                         ),hr(style="border-top: 1px dashed #B5C7DA;")),
                                                                  column(12, id = "MUVplotTab",
                                                                         tabsetPanel(id = "MUVtabset",
                                                                                      tabPanel("VIP plot", value = "MUVplotVIP",br(),
                                                                                               column(12,
                                                                                                      plotlyOutput(outputId = "plotVIP", width = "100%", height = "650px")
                                                                                               )
                                                                                      ),
                                                                                      tabPanel("Model evaluation", value = "MUVplotVAL1",br(),
                                                                                               column(12,
                                                                                                      plotlyOutput(outputId = "plotVAL1", width = "100%", height = "650px")
                                                                                               )
                                                                                      ),
                                                                                      tabPanel("Output table", value = "MUVoutputTab",br(),
                                                                                               h5("Output table contains order and average VIP rank of variables in minimal-optimal model."),
                                                                                               hr(style="border-top: 1px dashed #D3D3D3;"),
                                                                                               div(DT::dataTableOutput('INtableVIP'), style = "overflow-y:scroll; overflow-x:scroll;max-height: 600px; position:relative; max-height: 600px;"),br()
                                                                                      )
                                                                         )
                                                                  )

                                                                )
                                                      )
                                               ))
                                    )),

                           ###################################################################
                           ############################MBPL PAGE##############################
                           ###################################################################
                           tabPanel("", value = "mbpl",
                                    tabsetPanel(id = "MBPLTabset",
                                                 ###upload data###
                                                 tabPanel("Upload data", value = "mbplUP",
                                                          fluidRow(width=12,
                                                                   sidebarLayout(
                                                                     column(3,
                                                                            ##text input data set name
                                                                            sidebarPanel(width=12,style = "overflow-y:scroll; max-height: 900px; position:relative; height: 750px;",
                                                                                         h4("Start upload your data"),
                                                                                         br(),
                                                                                         selectInput("MBPLdmun","Number of data sets",choices=c(2,3,4,5), selected=2,multiple = F, width = 300),
                                                                                         h6("*The data must be in a CSV file, limit to 5 data sets."),
                                                                                         actionButton("NextMBPL1", "Data set 1",style="color: #464546; background-color: #c5d7d9;"),
                                                                                         fileInput("INfileMBPL1",h5("Upload a csv file of data set 1")),
                                                                                         textInput("nameMBPL1",label=h5("Data set name")),
                                                                                         selectInput("idMBPL1","Sample ID column",choices=NULL, selected=NULL),
                                                                                         selectInput("classMBPL1","Class/response column (Y)",choices=NULL, selected=NULL),
                                                                                         selectInput("firstMBPL1","1st variable/feature column (X)",choices=NULL, selected=NULL),
                                                                                         br(),
                                                                                         actionButton("NextMBPL2", "Data set 2",style="color: #464546; background-color: #c5d7d9;"),
                                                                                         fileInput("INfileMBPL2",h5("Upload a csv file of data set 2")),
                                                                                         textInput("nameMBPL2",label=h5("Data set name")),
                                                                                         selectInput("idMBPL2","Sample ID column",choices=NULL, selected=NULL),
                                                                                         selectInput("classMBPL2","Class/response column (Y)",choices=NULL, selected=NULL),
                                                                                         selectInput("firstMBPL2","1st variable/feature column (X)",choices=NULL, selected=NULL),
                                                                                         br(),
                                                                                         actionButton("NextMBPL3", "Data set 3",style="color: #464546; background-color: #c5d7d9;"),
                                                                                         fileInput("INfileMBPL3",h5("PUpload a csv file of data set 3")),
                                                                                         textInput("nameMBPL3",label=h5("Data set name")),
                                                                                         selectInput("idMBPL3","Sample ID column",choices=NULL, selected=NULL),
                                                                                         selectInput("classMBPL3","Class/response column (Y)",choices=NULL, selected=NULL),
                                                                                         selectInput("firstMBPL3","1st variable/feature column (X)",choices=NULL, selected=NULL),
                                                                                         br(),
                                                                                         actionButton("NextMBPL4", "Data set 4",style="color: #464546; background-color: #c5d7d9;"),
                                                                                         fileInput("INfileMBPL4",h5("Upload a csv file of data set 4")),
                                                                                         textInput("nameMBPL4",label=h5("Data set name")),
                                                                                         selectInput("idMBPL4","Sample ID column",choices=NULL, selected=NULL),
                                                                                         selectInput("classMBPL4","Class/response column (Y)",choices=NULL, selected=NULL),
                                                                                         selectInput("firstMBPL4","1st variable/feature column (X)",choices=NULL, selected=NULL),
                                                                                         br(),
                                                                                         actionButton("NextMBPL5", "Data set 5",style="color: #464546; background-color: #c5d7d9;"),
                                                                                         fileInput("INfileMBPL5",h5("Upload a csv file of data set 5")),
                                                                                         textInput("nameMBPL5",label=h5("Data set name")),
                                                                                         selectInput("idMBPL5","Sample ID column",choices=NULL, selected=NULL),
                                                                                         selectInput("classMBPL5","Class/response column (Y)",choices=NULL, selected=NULL),
                                                                                         selectInput("firstMBPL5","1st variable/feature column (X)",choices=NULL, selected=NULL),
                                                                                         br(),br(),
                                                                                         actionButton("UploadMBPL", "Upload", class = "btn btn-primary"),
                                                                                         actionButton("BackMBPL", "Back"),
                                                                                         actionButton("NextMBPL", "Next"),
                                                                                         br(),br(),
                                                                            )),
                                                                     column(9,
                                                                            mainPanel(width=12,br(),
                                                                                      fluidRow(
                                                                                        column(12, id = "MBPLsummary",
                                                                                               tabsetPanel(id = "MBPLsummaryT",
                                                                                                            tabPanel("Summary",value = "MBPLsummaryTab",
                                                                                                                     div(verbatimTextOutput("txtbox.mbpl"),style = "overflow-y:scroll; max-height: 350px;")
                                                                                                            )
                                                                                               ),hr(style="border-top: 1px dashed #B5C7DA;")),
                                                                                        column(12, id = "MBPLplotTab",
                                                                                               tabsetPanel(id = "MBPLtabDat",
                                                                                                            tabPanel("Data set 1", value = "MBPLdatS1",br(),
                                                                                                                     div(DT::dataTableOutput('INtableMBPL1'), style = "overflow-y:scroll; overflow-x:scroll;max-height: 600px; position:relative; max-height: 600px;"),br()
                                                                                                            ),
                                                                                                            tabPanel("Data set 2", value = "MBPLdatS2",br(),
                                                                                                                     div(DT::dataTableOutput('INtableMBPL2'), style = "overflow-y:scroll; overflow-x:scroll;max-height: 600px; position:relative; max-height: 600px;"),br()
                                                                                                            ),
                                                                                                            tabPanel("Data set 3", value = "MBPLdatS3",br(),
                                                                                                                     div(DT::dataTableOutput('INtableMBPL3'), style = "overflow-y:scroll; overflow-x:scroll;max-height: 600px; position:relative; max-height: 600px;"),br()
                                                                                                            ),
                                                                                                            tabPanel("Data set 4", value = "MBPLdatS4",br(),
                                                                                                                     div(DT::dataTableOutput('INtableMBPL4'), style = "overflow-y:scroll; overflow-x:scroll;max-height: 600px; position:relative; max-height: 600px;"),br()
                                                                                                            ),
                                                                                                            tabPanel("Data set 5", value = "MBPLdatS5",br(),
                                                                                                                     div(DT::dataTableOutput('INtableMBPL5'), style = "overflow-y:scroll; overflow-x:scroll;max-height: 600px; position:relative; max-height: 600px;"),br()
                                                                                                            )
                                                                                               )
                                                                                        )
                                                                                      )

                                                                            )
                                                                     )
                                                                   ))
                                                 ),
                                                 ###impute missing value before MBPLS-DA###
                                                 tabPanel("Missing value imputation", value = "missing_m",
                                                          fluidRow(width=12,
                                                                   sidebarLayout(
                                                                     column(width=3,
                                                                            sidebarPanel(width=12,style = "overflow-y:scroll; max-height: 900px; position:relative;",
                                                                                         h4("Missing value imputation"),
                                                                                         radioButtons("reall_m",
                                                                                                      h5("Remove all variables with missing values?"),
                                                                                                      choices = list("Remove all" = 'TRUE',
                                                                                                                     "Not Remove all" = 'FALSE'
                                                                                                      ),
                                                                                                      selected = 'FALSE'),
                                                                                         br(),
                                                                                         sliderInput("remPercent_m", h5("Remove variables with missing values > (%)"),
                                                                                                     min = 1, max = 99, value = 20),
                                                                                         br(),
                                                                                         radioButtons("missCheck_m", h5("Choose imputation method"),
                                                                                                      choices = list("Replace by minimum value" = "min",
                                                                                                                     "Mean" = "mean",
                                                                                                                     "Median" = "median",
                                                                                                                     "KNN" = "knn",
                                                                                                                     "BPCA" = "bpca",
                                                                                                                     "PPCA" = "ppca",
                                                                                                                     "SVD" = "svd",
                                                                                                                     "RF" = "rf"
                                                                                                      ),selected = 'median'),
                                                                                         br(),
                                                                                         actionButton("run_m", "Run", class = "btn btn-primary"),
                                                                                         actionButton("back_m", "Back"),
                                                                                         actionButton("next_m", "Next")

                                                                            )),
                                                                     column(width=9,
                                                                            mainPanel(width=12,br(),
                                                                                      fluidRow(
                                                                                        column(12,id = "MIsummary_m",
                                                                                               tabsetPanel(id = "MIsummaryT_m",
                                                                                                            tabPanel("Summary",value = "MIsummaryTab_m",
                                                                                                                     div(verbatimTextOutput("txtbox.missing_m"),style = "overflow-y:scroll; max-height: 350px;")
                                                                                                            )
                                                                                               ),hr(style="border-top: 1px dashed #B5C7DA;")),
                                                                                        column(12, id = "MIoutput_m",
                                                                                               tabsetPanel(id = "MIoutputT_m",
                                                                                                            tabPanel("Data set 1", value = "MBPLdatS1_m",br(),
                                                                                                                     div(DT::dataTableOutput('INtableMBPL1_m'), style = "overflow-y:scroll; overflow-x:scroll;max-height: 600px; position:relative; max-height: 600px;"),br()
                                                                                                            ),
                                                                                                            tabPanel("Data set 2", value = "MBPLdatS2_m",br(),
                                                                                                                     div(DT::dataTableOutput('INtableMBPL2_m'), style = "overflow-y:scroll; overflow-x:scroll;max-height: 600px; position:relative; max-height: 600px;"),br()
                                                                                                            ),
                                                                                                            tabPanel("Data set 3", value = "MBPLdatS3_m",br(),
                                                                                                                     div(DT::dataTableOutput('INtableMBPL3_m'), style = "overflow-y:scroll; overflow-x:scroll;max-height: 600px; position:relative; max-height: 600px;"),br()
                                                                                                            ),
                                                                                                            tabPanel("Data set 4", value = "MBPLdatS4_m",br(),
                                                                                                                     div(DT::dataTableOutput('INtableMBPL4_m'), style = "overflow-y:scroll; overflow-x:scroll;max-height: 600px; position:relative; max-height: 600px;"),br()
                                                                                                            ),
                                                                                                            tabPanel("Data set 5", value = "MBPLdatS5_m",br(),
                                                                                                                     div(DT::dataTableOutput('INtableMBPL5_m'), style = "overflow-y:scroll; overflow-x:scroll;max-height: 600px; position:relative; max-height: 600px;"),br()
                                                                                                            )
                                                                                               )

                                                                                        ))
                                                                            )
                                                                     )
                                                                   ))
                                                 ),
                                                 ###MBPLS-DA Analysis###
                                                 tabPanel("MBPLS-DA analysis", value = "mbplAna",
                                                          fluidRow(width=12,
                                                                   sidebarLayout(
                                                                     column(width=3,
                                                                            sidebarPanel(width=12,style = "overflow-y:scroll; max-height: 900px; position:relative; height: 750px;",
                                                                                         h4("MBPLS-DA analysis"),
                                                                                         tags$div(HTML("<p>*Implemented from packMBPLSDA package, see <a href='https://doi.org/10.1007/s11306-019-1598-y', target='_blank'>reference</a></p>")),
                                                                                         h6("*It might take long execution time depended on data set size, the number of repetitions and model testing process."),br(),
                                                                                         numericInput("MBPLnf", "Maximum number of scanned principle components", 10, min = 1, max=10),
                                                                                         numericInput("MBPLoptdim", "Optimal number of principle components", 2, min = 1, max=10),
                                                                                         numericInput("MBPLnboot", "Number of bootstrap replications", 10, min = 1),
                                                                                         numericInput("MBPLnrepet", "Number of cross-validation repetitions", 3, min = 1),
                                                                                         radioButtons("MBPLtest",
                                                                                                      h5("Perform model components testing and permutation?"),
                                                                                                      choices = list("Yes" = 'TRUE',
                                                                                                                     "No" = 'FALSE'
                                                                                                      ),
                                                                                                      selected = 'FALSE'),
                                                                                         numericInput("MBPLnpermut", "Number of permutations", 10, min = 1),
                                                                                         numericInput("MBPLthreshold", "Prediction threshold (max = 1)", 0.5, min = 0, max=1, step=0.1),
                                                                                         actionButton("RunMBPL", "Run", class = "btn btn-primary"),
                                                                                         actionButton("BackMBPL2", "Back"),br(),br(),
                                                                                         downloadButton("export_mbpls", "Export results", align = "center", class = "btn btn-success")
                                                                            )),
                                                                     column(9,
                                                                            mainPanel(width=12,br(),
                                                                                      fluidRow(
                                                                                        column(12, id = "MBPLsummaryAna",
                                                                                               tabsetPanel(id = "MBPLsummaryTAna",
                                                                                                            tabPanel("Summary",value = "MBPLsummaryTabAna",
                                                                                                                     div(verbatimTextOutput("txtbox.mbplAna"),style = "overflow-y:scroll; max-height: 350px;")
                                                                                                            )
                                                                                               ),hr(style="border-top: 1px dashed #B5C7DA;")),
                                                                                        column(12, id = "MBPLplotTabAna",
                                                                                               tabsetPanel(id = "MBPLtabsetAna",
                                                                                                            tabPanel("BIP plot", value = "MBPLbipDat",br(),
                                                                                                                     column(12,
                                                                                                                      plotlyOutput(outputId = "plotBipMBPL", width = "100%", height = "480px")
                                                                                                                     )
                                                                                                            ),
                                                                                                            tabPanel("VIP plot", value = "MBPLvip",br(),
                                                                                                                     column(12,
                                                                                                                      numericInput("MBPLselectvipprop", "Select top n% of variables to display (max = 100%)", 10, min = 5, max=100, step=5),
                                                                                                                      plotlyOutput(outputId = "plotVipMBPL", width = "100%", height = "480px")
                                                                                                                     )
                                                                                                            ),tabPanel("Loading plot", value = "MBPLload",br(),
                                                                                                                     column(12,
                                                                                                                      column(6,selectInput("MBPLselectPC","Select PC to display",choices=NULL, selected=NULL)),
                                                                                                                      column(6,numericInput("MBPLselectloadprop", "Select top n% of variables to display (max = 100%)", 10, min = 5, max=100, step=5)),
                                                                                                                      column(12,plotlyOutput(outputId = "plotLoadMBPL", width = "100%", height = "480px"))
                                                                                                                     )
                                                                                                            ),tabPanel("Model validation", value = "MBPLmodelVa",br(),
                                                                                                                     column(12,
                                                                                                                      plotlyOutput(outputId = "plotScreeMBPL", width = "90%", height = "90%"),br(),
                                                                                                                      plotlyOutput(outputId = "plotTestdimMBPL", width = "90%", height = "90%"),br(),
                                                                                                                      plotlyOutput(outputId = "plotPermutMBPL", width = "90%", height = "90%"),br()
                                                                                                                     )
                                                                                                            ),

                                                                                               )
                                                                                        )

                                                                                      )
                                                                            )
                                                                     ))
                                                          ))

                                    )),

                           ###################################################################
                           ############################Enrich PAGE############################
                           ###################################################################
                           tabPanel("",value = "DIP",
                                    tabsetPanel( id = "statTabset",
                                                 ##Over-representation analysis page######
                                                 tabPanel("Overrepresentation", value = "Overr",
                                                          fluidRow(width=12,
                                                                   sidebarLayout(
                                                                     column(width=3,
                                                                            sidebarPanel(
                                                                              width=12,style = "overflow-y:scroll; max-height: 850px; position:relative;",
                                                                              h4("Overrepresentation analysis (ORA)"),br(),
                                                                              textAreaInput("EntityIDs", "Enter list of variables:",height = 250),
                                                                              h6("*Use KEGG or HMDB for compounds, UniProt for proteins or Ensembl for genes."),
                                                                              h6("*Separate each ID by a new line."),
                                                                              fileInput("INfileOVER",h5("or Upload a csv file (one column table without header)")),
                                                                              radioButtons("OVERnode", h5("Choose variable type"),
                                                                                           choices = list(
                                                                                             "Compound" = "compound",
                                                                                             "Protein" = "protein",
                                                                                             "Gene" = "gene"
                                                                                           ),selected = "compound"),
                                                                              radioButtons("OVERset", h5("Choose annotation type"),
                                                                                           choices = list(
                                                                                             "KEGG pathway" = "pathway",
                                                                                             "HMDB chemical class" = "chemicalclass"
                                                                                           ),selected = "pathway"),
                                                                              h6("*KEGG pathway is for KEGG compounds, UniProt proteins or Ensembl genes."),
                                                                              h6("*HMDB chemical class is for HMDB compounds only."),br(),
                                                                              numericInput("OVERsize", "Minimum number of members", 3, min = 1),
                                                                              br(),
                                                                              actionButton("runOVER", "Run", class = "btn btn-primary"),
                                                                              actionButton("backOVER", "Back"),
                                                                              actionButton("revOVER", "Reset"),
                                                                              br(),br(),
                                                                              downloadButton("export_opa", "Export results", align = "center", class = "btn btn-success")
                                                                            )
                                                                     ),
                                                                     column(width=9,
                                                                            mainPanel(width=12,br(),
                                                                                      fluidRow(
                                                                                        column(12, id = "OVERsummary",
                                                                                               tabsetPanel(id = "OVERsummaryT",
                                                                                                            tabPanel("Summary", value = "OVERsummaryTab",
                                                                                                                     div(verbatimTextOutput("txtbox.Overr"),style = "overflow-y:scroll;  max-height: 350px;")
                                                                                                            )
                                                                                               ),hr(style="border-top: 1px dashed #B5C7DA;")),
                                                                                        column(12, id = "OVERoutput",
                                                                                               tabsetPanel( id = "OVERoutputT",
                                                                                                            tabPanel("Output table", value = "OVERoutputTab",br(),
                                                                                                                     div(DT::dataTableOutput('INtableOVER'), style = "overflow-y:scroll; overflow-x:scroll;max-height: 900px; position:relative;"),br()
                                                                                                            )
                                                                                               )
                                                                                        )

                                                                                      )
                                                                            ))
                                                                   ))

                                                 ),

                                                 ##### Enrichment analysis ######
                                                 tabPanel("Enrichment", value = "Enrich",
                                                          fluidRow(width=12,
                                                                   sidebarLayout(
                                                                     column(width=3,
                                                                            sidebarPanel(
                                                                              width=12,style = "overflow-y:scroll; max-height: 900px; position:relative;",
                                                                              h4("Enrichment analysis (ORA)"),
                                                                              tags$div(HTML("<p>*Implemented from piano package, see <a href='https://doi.org/10.1093/nar/gkt111', target='_blank'>reference</a></p>")),
                                                                              h6("*Use KEGG or HMDB for compounds, UniProt for proteins or Ensembl for genes."),
                                                                              h6("*Variables are in the first column."),br(),
                                                                              fileInput("INfileENR",h5("Upload a csv file")),
                                                                              selectInput("ENRpcol","P-value column",choices=NULL, selected=NULL),
                                                                              selectInput("ENRfccol","Fold-change column (optional)",choices=NULL, selected=NULL),
                                                                              radioButtons("ENRmed", h5("Choose analysis method"),
                                                                                           choices = list(
                                                                                             "Reporter" = "reporter",
                                                                                             "Fisher" = "fisher",
                                                                                             "Median" = "median",
                                                                                             "Mean" = "mean",
                                                                                             "Stouffer" = "stouffer"
                                                                                           ),selected = "reporter"),
                                                                              radioButtons("ENRnode", h5("Choose variable type"),
                                                                                           choices = list(
                                                                                             "Compounds" = "compound",
                                                                                             "Proteins" = "protein",
                                                                                             "Genes" = "gene"
                                                                                           ),selected = "compound"),
                                                                              radioButtons("ENRset", h5("Choose annotation type"),
                                                                                           choices = list(
                                                                                             "KEGG pathway" = "pathway",
                                                                                             "HMDB chemical class" = "chemicalclass"
                                                                                           ),selected = "pathway"),
                                                                              h6("*KEGG pathway is for KEGG compounds, UniProt proteins or Ensembl genes."),
                                                                              h6("*HMDB chemical class is for HMDB compounds only."),br(),
                                                                              numericInput("ENRsize", "Minimum number of members", 3, min = 1),br(),
                                                                              actionButton("runENR", "Run", class = "btn btn-primary"),
                                                                              actionButton("backENR", "Back"),
                                                                              actionButton("revENR", "Reset"),
                                                                              br(),br(),
                                                                              downloadButton("export_enr", "Export results", align = "center", class = "btn btn-success")
                                                                            )),
                                                                     column(width=9,
                                                                            mainPanel(width=12,br(),
                                                                                      fluidRow(
                                                                                        column(12, id = "ENRsummary",
                                                                                               tabsetPanel(id = "ENRsummaryT",
                                                                                                            tabPanel("Summary", value = "ENRsummaryTab",
                                                                                                                     div(verbatimTextOutput("txtbox.ENR"),style = "overflow-y:scroll;  max-height: 350px;")
                                                                                                            )
                                                                                               ),hr(style="border-top: 1px dashed #B5C7DA;")),
                                                                                        column(12, id = "ENRoutput",
                                                                                               tabsetPanel(id = "ENRoutputT",
                                                                                                            tabPanel("Output table", value = "ENRoutputTab",br(),
                                                                                                                     div(DT::dataTableOutput('INtableENR'), style = "overflow-y:scroll; overflow-x:scroll;max-height: 900px; position:relative;"),br()
                                                                                                            )
                                                                                               )
                                                                                        )

                                                                                      )
                                                                            )))
                                                          )
                                                 ),

                                                 ##### Combine ORA ######
                                                 tabPanel("Integrated ORA", value = "cORA",
                                                          fluidRow(width=12,
                                                                   sidebarLayout(
                                                                     column(width=3,
                                                                            sidebarPanel(width=12,style = "overflow-y:scroll; max-height: 600px; position:relative; ",
                                                                                         h4("Integrated pathway overrepresentation analysis"),
                                                                                         h6("*Use KEGG for compounds, UniProt for proteins or Ensembl for genes."),
                                                                                         h6("*The 1st column contains variables and the 2nd column contains variable types."),br(),
                                                                                         fileInput("INfileORA",h5("Upload a csv file")),
                                                                                         numericInput("ORAsize", "Minimum number of members", 3, min = 1),br(),
                                                                                         actionButton("runORA", "Run", class = "btn btn-primary"),
                                                                                         actionButton("backORA", "Back"),
                                                                                         actionButton("revORA", "Reset"),
                                                                                         br(),br(),
                                                                                         downloadButton("export_cora", "Export results", align = "center", class = "btn btn-success")
                                                                            )),
                                                                     column(width=9,
                                                                            mainPanel(width=12,br(),
                                                                                      fluidRow(
                                                                                        column(12, id = "ORAsummary",
                                                                                               tabsetPanel(id = "ORAsummaryT",
                                                                                                            tabPanel("Summary", value = "ORAsummaryTab",
                                                                                                                     div(verbatimTextOutput("txtbox.ORA"),style = "overflow-y:scroll;  max-height: 350px;")                                                                                                            )

                                                                                               ),hr(style="border-top: 1px dashed #B5C7DA;")),
                                                                                        column(12, id = "ORAoutput",
                                                                                             tabsetPanel(id = "ORAoutputT",
                                                                                                          tabPanel("Output table", value = "ORAoutputTab",br(),
                                                                                                                   div(DT::dataTableOutput('INtableORA'), style = "overflow-y:scroll; overflow-x:scroll;max-height: 900px; position:relative;"),br()
                                                                                                          )
                                                                                             )
                                                                                        ))
                                                                            )
                                                                     )))
                                                 )
                                    )


                           )
                )## navbarPage ##

)
