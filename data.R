#######################
# general information #
#######################

# file:       data.R
# author(s):  Marcel Schilling <marcel.schilling@mdc-berlin.de>
# created:    2017-02-23
# purpose:    load input data for tomo-seq shiny app


######################################
# change log (reverse chronological) #
######################################

# 2017-02-23: initial version (double-sourcing check & tomo-seq data loading)


#############
# libraries #
#############

# tomo-seq data is provided as tibble
require(tibble)


##############
# parameters #
##############

# load parameter definitions
source("params.R")


########
# data #
########

# ensure input data are not loaded already
if(!exists("input.data"))

  # load input data
  input.data<-

    # store input data in a named list
    list(

      # load tomo-seq data from file
      tomoseq.data=readRDS(params$tomoseq.data.file)

      # load input data list definition
      )

