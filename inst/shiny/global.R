app_jscode <-
  "shinyjs.disableTab = function(name) {
    var tab = $('.nav li a[data-value=' + name + ']');
    tab.bind('click.tab', function(e) {
      e.preventDefault();
      return false;
    });
    tab.addClass('disabled');
  }
  shinyjs.enableTab = function(name) {
    var tab = $('.nav li a[data-value=' + name + ']');
    tab.unbind('click.tab');
    tab.removeClass('disabled');
  }"

reset_all_evenData <-
    "shinyjs.resetSelected = function() { Shiny.setInputValue('plotly_selected-pca_isO', 'null', {priority: 'event'});
                              Shiny.setInputValue('plotly_selected-pca_is', 'null', {priority: 'event'});
                              Shiny.setInputValue('plotly_selected-pca_NO', 'null', {priority: 'event'});
                              Shiny.setInputValue('plotly_selected-pca_N', 'null', {priority: 'event'});
                              Shiny.setInputValue('plotly_click-A', 'null', {priority: 'event'});
                                                                          }"

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")}

install_pkgs('shiny')
install_pkgs('shinydashboard')
install_pkgs('shinyjs')
install_pkgs('dplyr')
install_pkgs('MUVR')
install_pkgs('V8')
install_pkgs('DT')
install_pkgs('ggplot2')
install_pkgs('plotly')
install_pkgs('ggthemes')
install_pkgs('shinythemes')
install_pkgs('bslib')
install_pkgs('stringr')
install_pkgs('shinybusy')
install_pkgs('shinyWidgets')
install_pkgs('ade4')
install_pkgs('packMBPLSDA')
install_pkgs('gridExtra')

library(shiny)
library(shinydashboard)
library(shinyjs)
library(dplyr)
library(MUVR)
library(V8)
library(DT)
library(ggplot2)
library(RColorBrewer)
library(plotly)
library(ggthemes)
library(shinythemes)
library(bslib)
library(stringr)
library(shinybusy)
library(shinyWidgets)
library(ade4)
library(packMBPLSDA)
library(gridExtra)
