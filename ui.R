#######################
# general information #
#######################

# file:       ui.R
# author(s):  Marcel Schilling <marcel.schilling@mdc-berlin.de>
# created:    2017-02-21
# purpose:    define front end for tomo-seq shiny app


######################################
# change log (reverse chronological) #
######################################

# 2017-02-21: initial version (app title only)


#############
# libraries #
#############

# get pipe operators
require(magrittr)


##############
# parameters #
##############

# load parameter definitions
source("params.R")


########################
# shiny user interface #
########################

# generate single page user interface with title panel

# generate title panel
titlePanel(

  # set app title
  title=params$app.title

  # end title panel definition
  ) %>%

# embed title panel in page
fluidPage %>%

# initialize user interface
shinyUI
