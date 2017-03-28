# tomo-seq shiny app data loading script
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

# file:         data.R
# author(s):    Marcel Schilling <marcel.schilling@mdc-berlin.de>
# created:      2017-02-23
# last update:  2017-03-28
# license:      GNU Affero General Public License Version 3 (GNU AGPL v3)
# purpose:      load input data for tomo-seq shiny app


######################################
# change log (reverse chronological) #
######################################

# 2017-03-28: added missing explicit magrittr loading
# 2017-02-24: added license comment
# 2017-02-23: added sample names extraction
#             initial version (double-sourcing check & tomo-seq data loading)


#############
# libraries #
#############

# tomo-seq data is provided as tibble
require(tibble)

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

# ensure input data are not loaded already
if(!exists("input.data"))

  # begin data loading
  {

    # load input data
    input.data<-

      # store input data in a named list
      list(

        # load tomo-seq data from file
        tomoseq.data=readRDS(params$tomoseq.data.file)

        # load input data list definition
        )

    # get sample names
    input.data$sample.names<-

      # take tomo-seq data
      input.data$tomoseq.data %$%

      # extract used sample names
      unique(sample.name)

  # end data loading
  }
