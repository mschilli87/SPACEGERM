# tomo-seq shiny app data loading script
# Copyright (C) 2017-2018  Marcel Schilling
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
# last update:  2018-04-04
# license:      GNU Affero General Public License Version 3 (GNU AGPL v3)
# purpose:      load input data for tomo-seq shiny app


######################################
# change log (reverse chronological) #
######################################

# 2018-04-04: added default shift/stretch file loading and formatting
# 2018-03-20: added slice width calculation
#             added gonad model loading
# 2017-04-18: fixed changelog comments (broken since 2017-03-29/baea8e9)
#             added gene type extraction
#             fixed copy/paste-error in comment
# 2017-04-13: replaced unique by distinct
# 2017-04-12: added missing changelog entry
#             made dplyr an explicit dependency
#             switched to tidy gene profile input data
# 2017-03-29: added gene profile loading (incl. sample descriptions & genotypes extraction)
# 2017-03-28: added missing explicit magrittr loading
# 2017-02-24: added license comment
# 2017-02-23: added sample names extraction
#             initial version (double-sourcing check & tomo-seq data loading)


#############
# libraries #
#############

# tomo-seq data is provided as tibble
require(tibble)

# get dlply
require(plyr)

# get %>% & distinct
# Note: dplyr must be loaded after plyr!
require(dplyr)


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
    input.data <- list(tomoseq.data = readRDS(params$tomoseq.data.file),
                       shift_stretch = readRDS(params$shift_stretch.file),
                       gene.profiles = readRDS(params$gene.profiles.file),
                       gonad.model = readRDS(params$gonad.model.file))

    # get slice width [% gonad arm]
    input.data$tomoseq.data %<>%
        group_by(sample.name, gene, transcript.name) %>%
        mutate(n.slices = n()) %>%
        ungroup %>%
        mutate(width.percent = 100 / n.slices) %>%
        select(-n.slices)

    # get sample names
    input.data$sample.names<-

      # take tomo-seq data
      input.data$tomoseq.data %>%

      # extract used sample names
      distinct(sample.name) %>%

      # convert single-column tibble to vector
      unlist %>%

      # drop names
      unname

    # get default sample shifts
    input.data$sample.shift.defaults <-
      data_frame(sample.name = input.data$sample.names) %>%
      left_join(input.data$shift_stretch) %>%
      mutate(
        shift.default = ifelse(is.na(shift.default),
                               params$sample.shifts.input.default,
                               shift.default),
        shift.default = ifelse(shift.default < params$sample.shifts.input.min,
                               params$sample.shifts.input.min, shift.default),
        shift.default = ifelse(shift.default > params$sample.shifts.input.max,
                               params$sample.shifts.input.max,
                               shift.default)) %$%
        setNames(shift.default, sample.name)

    # get default sample stretches
    input.data$sample.stretch.defaults <-
      data_frame(sample.name = input.data$sample.names) %>%
      left_join(input.data$shift_stretch) %>%
      mutate(stretch.default = ifelse(is.na(stretch.default),
                                params$sample.stretches.input.default,
                                stretch.default),
             stretch.default =
               ifelse(stretch.default < params$sample.stretches.input.min,
                      params$sample.stretches.input.min, stretch.default),
             stretch.default =
               ifelse(stretch.default > params$sample.stretches.input.max,
                      params$sample.stretches.input.max, stretch.default)) %$%
      setNames(stretch.default, sample.name)

    # get sample descriptions
    input.data$sample.descriptions<-

      # take gene profiles
      input.data$gene.profiles %>%

      # extract used sample descriptions
      distinct(sample.description) %>%

      # convert single-column tibble to vector
      unlist %>%

      # drop names
      unname

    # get genotypes per sample description
    input.data$genotypes<-

      # take gene profiles
      input.data$gene.profiles %>%

      # extract used genotype/sample description combinations
      distinct(sample.description,genotype) %>%

      # extract used genotypes per sample description
      dlply("sample.description",with,unique(genotype))

    # get gene types per sample description & genotype
    input.data$gene.types<-

      # take gene profiles
      input.data$gene.profiles %>%

      # extract used genotype/sample description/gene type combinations
      distinct(sample.description,genotype,gene.type) %>%

      # extract used genotypes/gene type combinations per sample description
      dlply("sample.description",distinct,genotype,gene.type) %>%

      # extract used gene types per sample description/genotype combination
      llply(dlply,"genotype",with,unique(gene.type))

  # end data loading
  }
