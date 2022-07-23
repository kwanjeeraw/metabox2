# Metabox 2.0
Metabox 2.0: A toolbox for thorough metabolomics data processing, analysis, integration and interpretatin.
The tool includes several methods for data processing, statistical analysis, biomarker analysis and data interpretation.
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

##(optional) Using graphical user interface, require shiny package ##
install.packages("shiny") #Skip this step, if shiny is alreay installed
launch_gui()

##Example data transformation ##
input_dat = read_input_data('filename')
input_obj = set_input_obj(input_dat,idCol=1,classCol=2,xCol=3)
output_obj = transform_input_data(input_obj, method="log10")
```

Updates
=========
#### version 1.0 (JUL 2022)
* Initial version

License
=========
[GNU General Public License (v3)](https://github.com/kwanjeeraw/grinnGUI/blob/master/LICENSE)
