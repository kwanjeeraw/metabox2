fluidRow(
  column(12,
         h2("About"),hr(style="border-top: 1px dashed #B5C7DA;"),
         fluidRow(
           h4("Metabox 2.0"),
           column(12,
                  tags$p("Metabox 2.0: A toolbox for thorough metabolomic data analysis, integration and interpretation.
                         It is a substantial update of the metabox R package, released in 2016.
                         The tool includes several methods for data processing, statistical analysis, biomarker analysis, integrative analysis and data interpretation.
                         This GUI supports a wide range of users, from bench biologists to experienced bioinformaticians.
                         It comes with an intuitive web interface for simple data analysis.
                         We recommend the R command line version for custom pipelines and other exclusive projects."),
                  tags$p("This project is supported by Program Management Unit for Human Resources & Institutional Development, Research and Innovation (PMU-B).")
           )),hr(style="border-top: 1px dashed #B5C7DA;"),
         fluidRow(
           h4("Updates"),
           column(12,
                  HTML('<h5>version 2.4 (MAR 2023)</h5>
                  <ul><li>Add example data sets for GUI version</li></ul>
                  <h5>version 2.3 (FEB 2023)</h5>
                  <ul><li>Summarize coefficient of variation (cv) and normality</li>
                  <li>Fix default scaling of PCA plot</li>
                  <li>Fix default color</li></ul>
                  <h5>version 2.2 (OCT 2022)</h5>
                  <ul><li>Fix bug when running MUVR</li></ul>
                  <h5>version 2.1 (SEP 2022)</h5>
                  <ul><li>Add imputation methods: zero, half-min</li>
                  <li>Allow scaling and block weighting by the block inertia</li>
                  <li>Allow missing values in normalization, transformation and scaling</li></ul>
                  <h5>version 2.0 (JUL 2022)</h5>
                  <ul><li>Initial release</li></ul>'),
           )),hr(style="border-top: 1px dashed #B5C7DA;"),
         fluidRow(
           h5("References"),
           column(12,
                  tags$p(a(href="https://github.com/kwanjeeraw/metabox2","Metabox 2.0"),
                         "- A toolbox for thorough metabolomic data analysis, integration and interpretation, see",
                         a(href="#","[Ref]")),
                  tags$p(a(href="https://github.com/kwanjeeraw/metabox","Metabox 1.0"),
                         "- A toolbox for metabolomic data analysis, interpretation and integrative exploration was released in 2016, see",
                         a(href="https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0171046","[Ref]"))
           )),hr(),
         fluidRow(
           column(12,align="center",
                  tags$div(tags$a(href="https://www.si.mahidol.ac.th/",tags$img(src="mulogo.png")),
                  tags$a(href="https://www.si.mahidol.ac.th/th/research-academics/research/simpc/",tags$img(src="simpclogo.png")),
                  tags$a(href="http://metsysbio.com/",tags$img(src="msblogo.png")))
          ))
  ))
