# Metabox 2.0
Metabox 2.0: A toolbox for thorough metabolomic data analysis, integration and interpretation. Metabox 2.0 is an updated version of the [R package Metabox](https://github.com/kwanjeeraw/metabox), released in 2016. The tool includes several methods for data processing, statistical analysis, biomarker analysis, integrative analysis and data interpretation. Metabox 2.0 supports a wide range of users, from bench biologists to experienced bioinformaticians. It comes with an intuitive web interface for simple data analysis. We recommend the R command line version for custom pipelines and other exclusive projects.

![demo](https://github.com/kwanjeeraw/metabox2/blob/main/inst/shiny/www/first_pipeline.png)

How to use
==========
Here are some alternative ways to use Metabox 2.0:
* [Install as a standard R package](#install-as-a-standard-r-package)
* [Use GUI version from Docker](#build-a-docker-image-and-deploy)
* [Use an online version](#use-an-online-version)

### Install as a standard R package
* Download and install [R software](https://www.r-project.org/)
* For R 4.4.1, Install metabox2 and required packages using the following commands:
```
##Install R dependencies ##
if (!require("BiocManager"))
    install.packages("BiocManager")
if (!require("remotes"))
  install.packages("remotes")
BiocManager::install("affy", update=FALSE, version="3.20")
BiocManager::install("pcaMethods", update=FALSE, version="3.20")
BiocManager::install("preprocessCore", update=FALSE, version="3.20")
BiocManager::install('impute', update=FALSE, version='3.20')
BiocManager::install("vsn", update=FALSE, version="3.20")
BiocManager::install("ropls", update=FALSE, version="3.20")
remotes::install_version('igraph',version='2.0.3',repos='https://cran.rstudio.org/')
BiocManager::install("piano", update=FALSE, version="3.20")
remotes::install_gitlab("CarlBrunius/MUVR")
install.packages("https://cran.r-project.org/src/contrib/Archive/MetNorm/MetNorm_0.1.tar.gz", repo=NULL, method = "libcurl")

##Install metabox2 ##
remotes::install_github("kwanjeeraw/metabox2", dependencies = "Imports", force = FALSE, upgrade ="never")

##Run metabox2
library(metabox2)

##Run data transformation ##
input_dat = read_input_data('filename')
input_obj = set_input_obj(input_dat,idCol=1,classCol=2,xCol=3)
output_transform = transform_input_data(input_obj, method="log10")

##Run multivariate analysis ##
input_dat = read_input_data('filename')
input_obj = set_input_obj(input_dat,idCol=1,classCol=2,xCol=3)
output_multiv = multiv_analyze(input_obj, method = "pca", scale = "standard")

##(Optional) Use graphical user interface ##
launch_gui()
```
### Use GUI version from Docker
* Use metabox2 GUI version from Docker:
1) Download and install [Docker](https://www.docker.com/)
2) Download a Dockerfile from [here](https://github.com/kwanjeeraw/metabox2/blob/main/Dockerfile)
3) Build and deploy metaboxweb Docker using the following commands:
```
cd [go to Dockerfile location]
docker build --no-cache=true --platform linux/x86_64 -t metaboxweb .
docker run --name mbdocker -p 8081:3838 metaboxweb
```
### Use an online version
* Use an online version from our servers:

* [server1](http://metabox.metsysbio.com:3838/metaboxweb/)
* [server2](http://metabox.metsysbio.com:3838/metaboxweb2/)
* [server3](http://metabox.metsysbio.com:3838/metaboxweb3/)
* [server4](http://metabox.metsysbio.com:3838/metaboxweb4/)
* [server5](http://metabox.metsysbio.com:3838/metaboxweb5/)

*Note: An online version is in high demand. Users might experience slow page loading. Currently, we are expanding our server and creating a portable docker image.*

Updates
=======
#### version 2.11 (DEC 2024)
* New combine_statplot function for command version
* Fix bugs for impute missing values in GUI version
#### version 2.10 (NOV 2024)
* New multiv_viploadingplot function
* New scoreplot with shape for command version
* New impute missing values by LOD and PCA overview
* Update annotation data
* Remove test_multinormality and test_multiequalvar functions
* Change to glog2 for rlaplot
* Change default no. of predI to 5
* Change default ptsize to 3
* Change pdf size of a report
* Change transparency of score plot
* Update and bugs fix for GUI version
* Fix bugs for enrichment and ORA analysis
#### version 2.9 (SEP 2024)
* Exclude liwong and LOESS normalization
* New messages
* install.packages('MetNorm') for ruv
* New check_serrf input data
#### version 2.8 (AUG 2024)
* Change number of QC for normalization by QC-based methods
#### version 2.7 (APR 2024)
* Change number of cross-validation segments for multivariate analysis
* Change figure title
#### version 2.6 (JUNE 2023)
* Change report location
* Set default package color
* Update MUVR to current version
* Fix bug when running univariate analysis
#### version 2.5 (MAY 2023)
* Update pathway data for enrichment analysis
#### version 2.4 (MAR 2023)
* Add example data sets for GUI version
#### version 2.3 (FEB 2023)
* Summarize coefficient of variation (cv) and normality
* Fix default scaling of PCA plot
* Fix default color
#### version 2.2 (OCT 2022)
* Fix bug when running MUVR
#### version 2.1 (SEP 2022)
* More imputation methods: zero, half-min
* Allow scaling and block weighting by the block inertia
* Allow missing values in normalization, transformation and scaling
#### version 2.0 (JUL 2022)
* Initial version

References
=========
- Wanichthanarak K, In-On A, Fan S, Fiehn O, Wangwiwatsin A, Khoomrung S (2024) Data processing solutions to render metabolomics more quantitative: case studies in food and clinical metabolomics using Metabox 2.0. GigaScience 13(2024): [giae005](https://10.1093/gigascience/giae005)
- Wanichthanarak K, Fan S, Grapov D, Barupal DK, Fiehn O (2017) Metabox: A Toolbox for Metabolomic Data Analysis, Interpretation and Integrative Exploration. PLOS ONE 12(1): [e0171046](https://doi.org/10.1371/journal.pone.0171046)

License
=========
[GNU General Public License (v3)](https://github.com/kwanjeeraw/metabox2/blob/master/LICENSE)
