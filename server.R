#######################
# general information #
#######################

# file:         server.R
# author(s):    Marcel Schilling <marcel.schilling@mdc-berlin.de>
# created:      2017-02-21
# last update:  2017-02-21
# purpose:      define back end for tomo-seq shiny app


######################################
# change log (reverse chronological) #
######################################

# 2017-02-21: added gene names output assignment
#             initial version (empty template)


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

    # assign gene names output
    output$gene.names<-

      # render gene names input
      renderText(input$gene.names)

  # end shiny server function definition
  } %>%

# initialize shiny server
shinyServer
