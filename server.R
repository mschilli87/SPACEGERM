#######################
# general information #
#######################

# file:         server.R
# author(s):    Marcel Schilling <marcel.schilling@mdc-berlin.de>
# created:      2017-02-21
# last update:  2017-02-23
# purpose:      define back end for tomo-seq shiny app


######################################
# change log (reverse chronological) #
######################################

# 2017-02-23: replaced gene names output assignment by profile plot output assignment
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

          # end profile plot generation
          )

        # end profile plot rendering
        )

  # end shiny server function definition
  } %>%

# initialize shiny server
shinyServer
