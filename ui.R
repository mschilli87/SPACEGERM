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
# last update:  2017-04-18
# license:      GNU Affero General Public License Version 3 (GNU AGPL v3)
# purpose:      define front end for tomo-seq shiny app


######################################
# change log (reverse chronological) #
######################################

# 2017-04-18: added dynamically generated gene type input panel
# 2017-04-11: added gene list file import panel
# 2017-04-10: added gene table XLSX export button
# 2017-04-06: added gene table output
# 2017-04-05: fixed code indentation
#             added heatmap options input panel
#             added gene cluster count input panel
# 2017-03-29: added heatmap tab panel (incl. sample description & genotype input & heatmap output
#             panels)
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

# get plotlyOutput
require(plotly)


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
    ) %>%

  # embed gene profiles page into tab panel
  tabPanel(

    # label gene profiles tab panel
    title=params$gene.profiles.tab.title

    # end gene profiles tab panel definition
    ) %>%

  # embed gene profiles tab panel in tabset panel
  tabsetPanel(

    # add heatmap tab panel
    tabPanel(

      # label heatmap tab panel
      title=params$heatmap.tab.title

      # add sidebar layout to heatmap tab panel
      ,sidebarLayout(

        # define sidebar panel for heatmap tab panel
        sidebarPanel(

          # add sample description input panel
          selectInput(

            # name sample description input
            inputId="sample.description"

            # label sample description input panel
            ,label=params$sample.description.input.label %>%

              # make label 3rd level header
              h3

            # set choices for sample description input panel
            ,choices=input.data$sample.descriptions

            # set default selection for sample description input panel
            ,selected=params$sample.description.input.default

            # end sample description input panel definition
            )

          # add dynamically generated genotype input panel
          ,uiOutput(

            # name genotype input panel output
            outputId="genotype.input"

            # end genotype input panel definition
            )

          # add dynamically generated gene type input panel
          ,uiOutput(

            # name gene type input panel output
            outputId="gene.type.input"

            # end gene type input panel definition
            )

          # add gene cluster count input panel
          ,numericInput(

            # name gene cluster count input
            inputId="nclust.genes"

            # label gene cluster count input panel
            ,label=params$nclust.genes.input.label %>%

              # make label 3rd level header
              h3

            # set minimal value for gene cluster count input panel
            ,min=params$nclust.genes.input.min

            # set default value for gene cluster count input panel
            ,value=params$nclust.genes.input.default

            # end plot columns count input panel definition
            )

          # add heatmap options input panel
          ,checkboxGroupInput(

            # name heatmap options input
            inputId="heatmap.options"

            # label heatmap options input panel
            ,label=params$heatmap.options.input.label %>%

              # make label 3rd level header
              h3

            # set choices for heatmap options input panel
            ,choices=params$heatmap.options

            # set default selection for heatmap options input panel
            ,selected=params$heatmap.options.input.default

            # end heatmap options input panel definition
            )

          # add gene list file import panel
          ,fileInput(

            # name gene list file input
            inputId="gene.list.file"

            # label gene list file import panel
            ,label=params$gene.list.import.label %>%

              # make label 3rd level header
              h3

            # set accepted MIME types for gene list file import
            ,accept=params$gene.list.file.import.mime.accept

            # label gene list file import button
            ,buttonLabel=params$gene.list.file.import.button.label

            # set gene list file import placeholder
            ,placeholder=params$gene.list.file.import.placeholder

            # end gene list file import
            )

          # end sidebar panel definition for heatmap tab panel
          )

        # define main panel for heatmap tab panel
        ,mainPanel(

          # generate heatmap output panel
          plotlyOutput(

            # name heatmaps output
            outputId="heatmap"

            # end heatmap output panel definition
            )

          # generate gene table output panel
          ,dataTableOutput(

            # name gene table output
            outputId="gene.table"

            # end gene table output panel definition
            )

          # generate gene table XLSX export button
          ,downloadButton(

            # name gene table XLSX export button
            outputId='gene.table.xlsx.export.button'

            # label gene table XLSX export button
            ,label=params$gene.table.xlsx.export.button.label

            # end gene table XLSX export button definition
            )

          # end main panel definition for heatmap tab panel
          )

        # end heatmap tab panel sidebar layout definition
        )

      # end heatmap tab panel definition
      )

    # end tabset panel definition
    )

  # end page definition
  ) %>%

# initialize user interface
shinyUI
