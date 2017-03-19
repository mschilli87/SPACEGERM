# tomo-seq shiny app user interface script
# Copyright (C) 2017  Marcel Schilling
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#######################
# general information #
#######################

# file:         ui.R
# author(s):    Marcel Schilling <marcel.schilling@mdc-berlin.de>
# created:      2017-02-21
# last update:  2017-03-19
# license:      GNU Affero General Public License Version 3 (GNU AGPL v3)
# purpose:      define front end for tomo-seq shiny app


######################################
# change log (reverse chronological) #
######################################

# 2017-03-19: added plot columns count input panel
#             fixed copy-and-paste error in comment
# 2017-02-24: added license comment
#             added dynamically generated sample shifts input panel
# 2017-02-23: added plot options input panel
#             added sample names input panel
#             replaced gene names output panel by profile plot output panel
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


########
# data #
########

# load input data
source("data.R")


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
  sidebarPanel(

    # add sample names input panel
    checkboxGroupInput(

      # name sample names input
      inputId="sample.names"

      # label gene names input panel
      ,label=params$sample.names.input.label %>%

        # make label 3rd level header
        h3

      # set choices for sample names input panel
      ,choices=input.data$sample.names

      # set default selection for sample names input panel
      ,selected=params$sample.names.input.default

      # end sample names input panel definition
      )

    # add plot options input panel
    ,checkboxGroupInput(

      # name plot options input
      inputId="plot.options"

      # label plot options input panel
      ,label=params$plot.options.input.label %>%

        # make label 3rd level header
        h3

      # set choices for plot options input panel
      ,choices=params$plot.options

      # set default selection for plot options input panel
      ,selected=params$plot.options.input.default

      # end plot options input panel definition
      )

    # add plot columns count input panel
    ,numericInput(

      # name plot columns count input
      inputId="ncols.plot"

      # label plot columns count input panel
      ,label=params$ncols.plot.input.label %>%

        # make label 3rd level header
        h3

      # set minimal value for plot columns count input panel
      ,min=params$ncols.plot.input.min

      # set default value for plot columns count input panel
      ,value=params$ncols.plot.input.default

      # end plot columns count input panel definition
      )

    # label sample shifts input panel
    ,params$sample.shifts.input.label %>%

        # make label 3rd level header
        h3

    # add dynamically generated sample shifts panel
    ,uiOutput(

      # name sample shifts input panel output
      outputId="shifts.input"

      # end sample shifts input panel definition
      )

    ) %>%

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
