# Metabox 2.0
Metabox 2.0: A toolbox for thorough metabolomics data processing, analysis, integration and interpretation. The tool includes several methods for data processing, statistical analysis, biomarker analysis and data interpretation.

![demo](metabox2_img.png)

Installation
============
* Install metabox2 and required packages using the following commands
```
##Install metabox2 ##
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("kwanjeeraw/metabox2")
library(metabox2)

##(Optional) Use graphical user interface ##
install.packages("shiny") #Skip this step, if shiny package is alreay installed
launch_gui()

##Run data transformation ##
input_dat = read_input_data('filename')
input_obj = set_input_obj(input_dat,idCol=1,classCol=2,xCol=3)
output_transform = transform_input_data(input_obj, method="log10")

##Run multivariate analysis ##
input_dat = read_input_data('filename')
input_obj = set_input_obj(input_dat,idCol=1,classCol=2,xCol=3)
output_multiv = multiv_analyze(input_obj, method = "pca", scale = "standard")
```

Updates
=========
#### version 1.0 (JUL 2022)
* Initial version

References
=========
- Wanichthanarak K, Fan S, Grapov D, Barupal DK, Fiehn O (2017) Metabox: A Toolbox for Metabolomic Data Analysis, Interpretation and Integrative Exploration. PLOS ONE 12(1): e0171046. [![DOI](https://doi.org/10.1371/journal.pone.0171046)

License
=========
[GNU General Public License (v3)](https://github.com/kwanjeeraw/metabox2/blob/master/LICENSE)
