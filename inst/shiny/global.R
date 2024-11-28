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
                              Shiny.setInputValue('plotly_selected-pca_mi', 'null', {priority: 'event'});
                              Shiny.setInputValue('plotly_selected-pca_miO', 'null', {priority: 'event'});
                              Shiny.setInputValue('plotly_selected-pca_isO', 'null', {priority: 'event'});
                              Shiny.setInputValue('plotly_selected-pca_is', 'null', {priority: 'event'});
                              Shiny.setInputValue('plotly_selected-pca_NO', 'null', {priority: 'event'});
                              Shiny.setInputValue('plotly_selected-pca_N', 'null', {priority: 'event'});
                              Shiny.setInputValue('plotly_click-A', 'null', {priority: 'event'});
                                                                          }"

library(shiny)
library(shinydashboard)
library(shinyjs)
library(dplyr)
#library(MUVR)
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
#library(ade4)
#library(packMBPLSDA)
library(gridExtra)
