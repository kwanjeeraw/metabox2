# Base R Shiny image
FROM rocker/shiny:4.4.1

# Install git
RUN apt-get update && apt-get install -y git

# Make a directory in the container
RUN mkdir /home/metabox2

# Install metabox2 and R dependencies
RUN install2.r --error --skipinstalled \
    remotes BiocManager \
    && rm -rf /tmp/downloaded_packages

RUN R -e "BiocManager::install('affy', update=FALSE, version='3.20')"
RUN R -e "BiocManager::install('pcaMethods', update=FALSE, version='3.20')"
RUN R -e "BiocManager::install('preprocessCore', update=FALSE, version='3.20')"
RUN R -e "BiocManager::install('impute', update=FALSE, version='3.20')"
RUN R -e "BiocManager::install('vsn', update=FALSE, version='3.20')"
RUN R -e "BiocManager::install('ropls', update=FALSE, version='3.20')"
RUN R -e "remotes::install_version('igraph',version='2.0.3',repos='https://cran.rstudio.org/')"
RUN R -e "BiocManager::install('piano', update=FALSE, version='3.20')"
RUN R -e "remotes::install_gitlab('CarlBrunius/MUVR', dependencies = 'Imports', version = '0.0.976', force = FALSE, upgrade = 'never')"
RUN R -e "remotes::install_github('kwanjeeraw/metabox2', dependencies = 'Imports', force = FALSE, upgrade ='never')"

# COPY MetNorm_0.1.tar.gz /home/metaboxr/MetNorm_0.1.tar.gz
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/MetNorm/MetNorm_0.1.tar.gz', dependencies = 'Imports', repos = NULL, type = 'source')"

# COPY DiscriMiner_0.1-29.tar.gz /home/metaboxr/DiscriMiner_0.1-29.tar.gz
#RUN R -e "install.packages('http://cran.nexr.com/src/contrib/DiscriMiner_0.1-29.tar.gz', dependencies = 'Imports', repos = NULL, type = 'source')"

# Copy the Shiny app code
#COPY metaboxweb/ /srv/shiny-server/metaboxweb/
RUN git clone https://github.com/kwanjeeraw/metaboxweb.git --depth 1 --branch=main /srv/shiny-server/metaboxweb/

RUN sudo chown -R shiny:shiny /srv/shiny-server/metaboxweb/

EXPOSE 3838

# Run the R Shiny app
CMD ["/usr/bin/shiny-server"]
