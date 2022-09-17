fluidRow(
  column(12,
         #includeMarkdown("../../README.md")
         h2("About"),hr(style="border-top: 1px dashed #B5C7DA;"),
         fluidRow(
           h4("Metabox 2.0"),
           column(12,
                  tags$p("It is an R package for comprehensive metabolomics data analysis.","It is an R package for comprehensive metabolomics data analysis.","It is an R package for comprehensive metabolomics data analysis.")
           )),hr(style="border-top: 1px dashed #B5C7DA;"),
         fluidRow(
           h5("News"),
           column(12,
                  tags$p(a(href="https://github.com/kwanjeeraw/metabox2","Metabox 2.0"),
                         "- An updated version for comprehensive metabolomics data analysis.",
                         a(href="http://commons.wikimedia.org/wiki/User:Sfoskett","[Ref]")),
                  tags$p(a(href="https://github.com/kwanjeeraw/metabox","Metabox 1.0"),
                         "- A toolbox for metabolomic data analysis, interpretation and integrative exploration was released in 2016.",
                         a(href="https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0171046","[Ref]"))
           ))
  ))
