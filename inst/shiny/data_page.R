fluidRow(
  column(12,
         h2("Example data sets"),hr(style="border-top: 1px dashed #B5C7DA;"),
         fluidRow(
         h4("Data processing and analysis pipeline"),
         column(12,
                tags$p(a(href="GCGC_DM_Samples.csv","GCGC_DM"),
                       "- plasma samples from diabetes mellitus cohort: 15 Pooled-QC samples and 2 ISs, 60 samples x 89 metabolites, see",
                       a(href="https://doi.org/10.1016/j.jbc.2022.102445","[Ref]")),
                tags$p(a(href="GC_Milk_Samples.csv","GC_Milk"),
                       "- retail milk samples: 10 Pooled-QC samples and 1 IS, 87 samples x 16 fatty acids, see",
                       a(href="https://doi.org/10.1016/j.jfca.2020.103785","[Ref]")),
                tags$p(a(href="LC_LN_Samples.csv","LC_LN"),
                       "- urine samples from lupus nephritis cohort: 1 IS, 116 samples x 8 metabolites, see",
                       a(href="https://doi.org/10.1016/j.isci.2021.103355","[Ref]"))
         )),
         fluidRow(
         h4("Data integration pipeline"),
         column(12,
                tags$p(a(href="GC_Fat_Tissue.csv","GC_Fat_Tissue"),
                       "and",
                       a(href="LC_Fat_Tissue.csv","LC_Fat_Tissue"),
                       "- adipose tissue samples from colorectal carcinoma cohort, metabolites are measured by 2 platforms, see",
                       a(href="https://doi.org/10.3945/ajcn.114.103804","[Ref]"))
         )),
         fluidRow(
         h4("Data interpretation pipeline"),
         column(12,
                tags$p(a(href="metabolite_list.csv","metabolite_list"),
                       "- table of metabolites, p-values, and fold-changes for ORA or enrichment analysis"),
                tags$p(a(href="metabolite_protein_list.csv","metabolite_protein_list"),
                       "- list of metabolites and proteins for integrated pathway analysis")
         ))
  ))
