fluidRow(
  column(12,
         fluidRow(column(5,h2("About")),
                  column(7,align="right",
                         tags$div(tags$a(href="https://www.si.mahidol.ac.th/",tags$img(src="mulogo.png")),
                                  tags$a(href="https://www.si.mahidol.ac.th/th/research-academics/research/simpc/",tags$img(src="simpclogo.png")),
                                  tags$a(href="http://metsysbio.com/",tags$img(src="msblogo.png")))
                  )),hr(style="border-top: 1px dashed #B5C7DA;"),
         fluidRow(
           h4("Metabox 2.0"),
           column(12,
                  tags$p("Metabox 2.0: A toolbox for thorough metabolomic data analysis, integration and interpretation.
                  Metabox 2.0 is an updated version of the R package ",
                  a(href="https://github.com/kwanjeeraw/metabox","Metabox"),", released in 2016.
                  The tool includes several methods for data processing, statistical analysis, biomarker analysis, integrative analysis and data interpretation.
                  This GUI supports a wide range of users, from bench biologists to experienced bioinformaticians.
                  It comes with an intuitive web interface for simple data analysis.
                  We recommend the R command line version for custom pipelines and other exclusive projects."),
                  #tags$p("This project is supported by Program Management Unit for Human Resources & Institutional Development, Research and Innovation (PMU-B).")
           )),hr(style="border-top: 1px dashed #B5C7DA;"),
         fluidRow(
           h4("Updates"),
           column(12,
                  HTML('<h5>version 2.11 (DEC 2024)</h5>
                  <ul><li>New combine_statplot function for command version</li>
                  <li>Fix bugs for impute missing values in GUI version</li></ul>
                  <h5>version 2.10 (NOV 2024)</h5>
                  <ul><li>New multiv_viploadingplot function</li>
                  <li>New scoreplot with shape for command version</li>
                  <li>New impute missing values by LOD and PCA overview</li>
                  <li>Update annotation data</li>
                  <li>Remove test_multinormality and test_multiequalvar functions</li>
                  <li>Change to glog2 for rlaplot</li>
                  <li>Change default no. of predI to 5</li>
                  <li>Change default ptsize to 3</li>
                  <li>Change pdf size of a report</li>
                  <li>Change transparency of score plot</li>
                  <li>Update and bugs fix for GUI version</li>
                  <li>Fix bugs for enrichment and ORA analysis</li></ul>
                  <p><a href="https://github.com/kwanjeeraw/metabox2">See all updates</a>.</p>'),
           )),hr(style="border-top: 1px dashed #B5C7DA;"),
         fluidRow(
           h4("References"),
           column(12,
                  tags$p(a(href="https://github.com/kwanjeeraw/metabox2","Metabox 2.0"),
                         "- A toolbox for thorough metabolomic data analysis, integration and interpretation, see",
                         a(href="https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giae005/7629842","[Ref]")),
                  tags$p(a(href="https://github.com/kwanjeeraw/metabox","Metabox 1.0"),
                         "- A toolbox for metabolomic data analysis, interpretation and integrative exploration was released in 2016, see",
                         a(href="https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0171046","[Ref]"))
           ))
  ))
