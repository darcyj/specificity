# get packages required by specificity.shiny
library(shiny)
library(ggplot2)
library(DT)
library(colourpicker)

# source specificity.shiny portable functions
source('funs.r')

# load data for the app
load('appdata.rdata')
# run visualization
plot_specs_shiny(sl, fd, fd_id_col)
