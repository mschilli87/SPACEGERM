# tomo-seq shiny app server script
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

# file:         server.R
# author(s):    Marcel Schilling <marcel.schilling@mdc-berlin.de>
# created:      2017-02-21
# last update:  2017-02-24
# license:      GNU Affero General Public License Version 3 (GNU AGPL v3)
# purpose:      define back end for tomo-seq shiny app


######################################
# change log (reverse chronological) #
######################################

# 2017-02-24: added license comment
#             added dynamic sample shift input panel assignment & corresponding user input support
# 2017-02-23: added user specified plot option selection
#             added user specified sample selection
#             replaced gene names output assignment by profile plot output assignment
# 2017-02-21: added gene names output assignment
#             initial version (empty template)

#############
# libraries #
#############

# get pipe operators
require(magrittr)


#############
# functions #
#############

# load functions
source("functions.R")


########
# data #
########

# load input data
source("data.R")


################
# shiny server #
################

# define shiny server function parameters
function(

  # app input list
  input

  # app output list
  ,output

  # end shiny server function parameter definition
  )

  # begin shiny server function definition
  {

    # assign sample shifts input panel output
    output$shifts.input<-

      # render sample shifts input panel
      renderUI(

        # take sample names to include in plot
        input$sample.names %>%

        # generate sample shifts input panel
        generate.sample.shifts.input

        # end sample shifts input panel rendering
        )

    # assign profile plot output
    output$profile.plot<-

      # render gene profiles plot
      renderPlot(

        # take tomo-seq data
        input.data$tomoseq.data %>%

        # generate profile plot
        generate.profile.plot(

          # plot profiles of genes specified by the user
          gene.names=input$gene.names

          # include sample specified by the user in plot
          ,sample.names=input$sample.names

          # set plot options specified by the user
          ,plot.options=input$plot.options

          # set sample shifts specified by the user
          ,sample.shifts=

            # take sample names of samples included in plot
            input$sample.names %>%

            # extract corresponding sample shifts specified by user
            get.sample.shifts(input)

          # end profile plot generation
          )

        # end profile plot rendering
        )

  # end shiny server function definition
  } %>%

# initialize shiny server
shinyServer