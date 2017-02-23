#######################
# general information #
#######################

# file:         ui.R
# author(s):    Marcel Schilling <marcel.schilling@mdc-berlin.de>
# created:      2017-02-21
# last update:  2017-02-23
# purpose:      define front end for tomo-seq shiny app


######################################
# change log (reverse chronological) #
######################################

# 2017-02-23: replaced gene names output panel by profile plot output panel
# 2017-02-21: added gene names input/output panels
#             initial version (app title only)


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
fluidPage(

  # add gene name input panel
  textInput(

    # name gene names input
    inputId="gene.names"

    # label gene names input panel
    ,label=params$gene.names.input.label %>%

      # make label 3rd level header
      h3

    # set default input gene names
    ,value=params$gene.names.input.default

    # set gene names input placeholder
    ,placeholder=params$gene.names.input.placeholder

    # end gene names input panel definition
    ) %>%

  # embed gene name input panel in sidebar
  sidebarPanel %>%

  # initialize side bar layout
  sidebarLayout(

    # generate gene names output panel
    plotOutput(

      # name gene names output
      outputId="profile.plot"
      ) %>%

    # embed gene names output panel in main panel
    mainPanel

    # end layout definition
    )

  # end page definition
  ) %>%

# initialize user interface
shinyUI
