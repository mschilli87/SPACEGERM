# tomo-seq shiny app function definition script
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

# file:         functions.R
# author(s):    Marcel Schilling <marcel.schilling@mdc-berlin.de>
# created:      2017-02-23
# last update:  2017-03-29
# license:      GNU Affero General Public License Version 3 (GNU AGPL v3)
# purpose:      define functions for tomo-seq shiny app


######################################
# change log (reverse chronological) #
######################################

# 2017-03-29: added genotype input panel generation and heatmap functions
# 2017-03-19: added plot columns count support
# 2017-03-02: made single y-scale for all sub-plots optional
#             added fixed x-axis limits support
# 2017-02-24: added license comment
#             added sample shifts input panel generation & input extraction functions / added sample
#             based shift assignment / added shift support to profile plot function
# 2017-02-23: re-ordered plot options
#             made raw data points optional
#             made raw data lines optional
#             added per-sample smooth fits
#             made across-sample smooth fit optional
#             added logarithmic y-axis scaling
#             initial version (profile plot generation function)


#############
# libraries #
#############

# get pipe operators
require(magrittr)

# get ggplot for plotting
require(ggplot2)

# get log2_trans for logscale
require(scales)

# get llply & alply
require(plyr)

# get heatmaply
require(heatmaply)


##############
# parameters #
##############

# load parameter definitions
source("params.R")


#############
# functions #
#############

  ####################
  # helper functions #
  ####################

# get sample shift input names by sample name
get.sample.shift.inputId<-

  # define sample shift input name function
  function(

    # sample name to get sample shift input name for
    sample.name

    # end sample shift input name function parameter definition
    )

    # begin sample shift input name function definition
    {

      # take sample name to get sample shift input name for
      sample.name %>%

      # add sample shift input name prefix
      paste0("sample.shift.",.)

    # end sample shift input name function definition
    }

# generate sample shift input panel
generate.sample.shift.input<-

  # define sample shift input panel generation function
  function(

    # sample name of generate sample shift input panel for
    sample.name

    # end sample shift input panel generation function parameter definition
    )

    # begin sample shift input panel generation function definition
    {

      # take sample name of generate sample shift input panel for
      sample.name %>%

      # generate sample shift input panel
      sliderInput(

        # name sample shift input value
        inputId=get.sample.shift.inputId(.)

        # label sample shift input panel
        ,label=.

        # set sample shift input panel minimum value
        ,min=params$sample.shifts.input.min

        # set sample shift input panel maximum value
        ,max=params$sample.shifts.input.max

        # set sample shift input panel default value
        ,value=params$sample.shifts.input.default

        # set sample shift input panel value suffix
        ,post=params$sample.shifts.input.suffix

        # end sample shift input panel generation
        )

    # end sample shift input panel generation function definition
    }

# parse gene names input from single string to vector
parse.gene.names<-

  # define gene names parsing function
  function(

    # input gene names string
    gene.names.input

    # end gene names parsing function parameter definition
    )

    # begin gene names parsing function definition
    {

      # take single string gene names input
      gene.names.input %>%

      # split single string into vector
      strsplit(params$gene.names.input.separator) %>%

      # strsplit returns a list (single element in this case) --> convert to plain vector
      unlist

    # end gene names parsing function definition
    }


# filter data by sample names
filter.data.by.sample.names<-

  # define filter by sample names function
  function(

    # unfiltered data
    unfiltered.data

    # sample names to keep
    ,sample.names.to.keep

    # end filter by sample names function parameter definition
    )

    # begin filter by sample names function definition
    {

      # take unfiltered data
      unfiltered.data %>%

      # subset unfiltered data
      subset(

        # select data with sample names to keep
        sample.name %in% sample.names.to.keep

        # end data subsetting
        )

    # end filter by sample names function definition
    }


# filter data by gene names
filter.data.by.gene.names<-

  # define filter by gene names function
  function(

    # unfiltered data
    unfiltered.data

    # gene names to keep
    ,gene.names.to.keep

    # end filter by gene names function parameter definition
    )

    # begin filter by gene names function definition
    {

      # take unfiltered data
      unfiltered.data %>%

      # subset unfiltered data
      subset(

        # select data with sample names to keep
        gene.name %in% gene.names.to.keep

        # end data subsetting
        )

    # end filter by gene names function definition
    }

# convert named sample shift list (or NULL) to named (default) sample shift vector
parse.sample.shifts<-

  # define sample shift conversion function
  function(

    # list with sample shifts labeled with sample names (or NULL)
    shifts

    # sample names to include in output sample shift vector
    ,samples

    # end sample shift conversion function parameter definition
    )

    # begin sample shift conversion function definition
    {

      # if no sample shifts were specified, use default values
      if(is.null(shifts))

        # take default sample shift input value
        params$sample.shifts.input.default %>%

        # repeat for each sample specified
        rep(times=length(samples)) %>%

        # label default sample shifts with sample names
        setNames(samples)

      # if sample shifts were specified, convert list to vector
      else

        # take sample names to include in output sample shift vector
        samples %>%

        # extract corresponding sample shifts
        shifts[.] %>%

        # convert names list to vector
        unlist

    # end sample shift conversion function definition
    }

# add shift column to data
add.shift.column<-

  # define shift column addition function
  function(

    # data to add shift column to
    input.data

    # shifts to add to data labeled by sample name
    ,shifts.by.sample

    # end shift column addition function parameter definition
    )

    # begin shift column addition function definition
    {

      # take data to add shift column to
      input.data %>%

      # add column to data
      cbind(

        # label shift column
        shift=

          # assign shift by sample name
          shifts.by.sample[.$sample.name]

        # and column addition
        )

    # end shift column addition function definition
    }

# plot gene profiles
plot.profiles<-

  # define profile plot function
  function(

    # plot data
    plot.data

    # show raw data points by default
    ,raw.points=T

    # don't show raw data lines by default
    ,raw.lines=F

    # show per-sample smooth fits by default
    ,smooth.each=T

    # show across-sample smooth fit by default
    ,smooth.pooled=T

    # log-transform y-axis by default
    ,logscale=T

    # don't fix x-axis limits by default
    ,fix.xlim=F

    # use a single y-scale for all sub-plots by default
    ,single.y.scale=T

    # use default plot columns count by default
    ,ncols=params$ncols.plot.input.default

    # end profile plot function parameter definition
    )

    # begin profile plot function definition
    {

      # create profile plot
      profile.plot<-

        # take plot data
        plot.data %>%

        # generate plot object
        ggplot(

          # bind data variables to plot parameters
          aes(

            # position points on the x-axis according to the (shifted) center of the corresponding
            # slice (percent distal-to-proximal)
            x=percent.center+shift

            # on the y-axis, plot the gene abundance (counts per million reads mapped)
            ,y=cpm

            # group data by sample name & color points/lines accordingly
            ,color=sample.name

            # end data derived plot parameter definition
            )

          # end general plot parameter definition
          ) %>%

        # split plot into panels
        + facet_wrap(

            # generate one sub-plot per gene name
            ~gene.name

            # adjust the (max.) number of sub-plots to put next to each other
            ,ncol=ncols

            # set sub-plot scales
            ,scales=

              # take single y-scale parameter
              single.y.scale %>%

              # set sub-plot scales based on single y-scale parameter
              ifelse(

                # if single y-scale requested: fix both scales for sub-plots
                "fixed"

                # if single y-scale not requested: fix only x-scale for sub-plots
                ,"free_y"

                # end single y-scale parameter evaluation
                )

            # end sub-plot definition
            ) %>%

        # color non-data plot elements in black, white & shades of grey only
        + theme_bw(

            # adjust plot base font sizes
            base_size=params$profile.plot.fontsize.base

            # end general plot style parameter definition
            ) %>%

        # add caption to plot
        + ggtitle(

            # set plot title
            label=params$profile.plot.title

            # end plot caption parameter definition
            ) %>%

        # adjust x-axis labelling
        + xlab(

            # adjust x-axis label
            label=params$profile.plot.xlab

            # end x-axis label parameter definition
            ) %>%

        # adjust y-axis labelling
        + ylab(

            # adjust y-axis label
            label=params$profile.plot.ylab

            # end y-axis label parameter definition
            ) %>%

        # adjust color (i.e. sample name) palette/legend
        + scale_color_brewer(

            # adjust sample name -> color mapping
            palette=params$profile.plot.brewer.palette

            # adjust color legend label
            ,name=params$profile.plot.sample.legend.label

            # end color palette/legend parameter definition
            )

        # add raw data points if specified
        if(raw.points)

          # modify profile plot
          profile.plot %<>%

          # plot raw data as points
          + geom_point(

              # adjust data point size
              size=params$profile.plot.pointsize

              # end data point parameter definition
              )

        # add raw data lines if specified
        if(raw.lines)

          # modify profile plot
          profile.plot %<>%

          # plot per-slice data as lines
          + geom_line(

              # adjust data line style
              linetype=params$profile.plot.linetype.raw

              # adjust data line size
              ,size=params$profile.plot.linesize.each

              # end data line parameter definition
              )

        # add per-sample smooth fits if specified
        if(smooth.each)

          # modify profile plot
          profile.plot %<>%

          # fit smooth line across samples
          + geom_smooth(

              # adjust smoothing method
              ,method=params$profile.plot.smoothing.method

              # adjust smooth line style
              ,linetype=params$profile.plot.linetype.smooth

              # adjust smooth line size
              ,size=params$profile.plot.linesize.each

              # end smooth line parameter definition
              )

        # add across-sample smooth fit if specified
        if(smooth.pooled)

          # modify profile plot
          profile.plot %<>%

          # fit smooth line across samples
          + geom_smooth(

              # pool data across samples & adjust smooth line color
              color=params$profile.plot.color.smooth

              # adjust smoothing method
              ,method=params$profile.plot.smoothing.method

              # adjust smooth line style
              ,linetype=params$profile.plot.linetype.smooth

              # adjust smooth line size
              ,size=params$profile.plot.linesize.pooled

              # end smooth line parameter definition
              )

        # scale y-axis logarithmically if specified
        if(logscale)

          # modify profile plot
          profile.plot %<>%

          # log2-transform y-axis
          + scale_y_continuous(trans=log2_trans())

        # fix x-axis limits if specified
        if(fix.xlim)

          # modify profile plot
          profile.plot %<>%

          # fix x-axis limits
          + xlim(params$profile.plot.xlim)

        # return profile plot
        profile.plot

    # end profile plot function definition
    }


# only top variable genes
keep.top.genes<-

  # define top variable genes filter function
  function(

    # expression matrix (genes as rows)
    mat

    # maximum number of genes to keep
    ,nmax=

      #by default, keep default maximum number of genes
      params@heatmap.max.ngenes

    # end top variable genes filter function parameter definition
    )

    # begin top variable genes filter function definition
    {

      # calculate standard deviations
      sds<-

        # take expression matrix
        mat %>%

        #  calculate row standard deviations
        alply(1,sd,na.rm=TRUE) %>%

        #  convert list to vector
        unlist

      # set missing standard deviations to zero
      sds[is.na(sds)]<-0

      # keep only variable genes
      variable<-sds!=0
      mat<-mat[variable,]
      sds<-sds[variable]

      # take standard deviations
      sds %>%

      # sort from high to low
      order(decreasing=TRUE) %>%

      # keep up to specified maximum number of highest varying genes
      head(nmax) %>%

      # return corresponding expression matrix rows
      {mat[.,]}

    # end top variable genes filter function definition
    }


# calculate distance matrix of matrix columns
get.column.distances<-

  # define matrix columns distance calculation function
  function(

    # matrix to get column distance of
    column.matrix

    # end matrix columns distance calculation function parameter definition
    )

    # begin matrix columns distance calculation function definition
    {

      # take matrix to get column distance of
      column.matrix %>%

      # calculate (pairwise) column (Pearson) correlations
      cor(use="complete") %>%

      # convert column similarities to distances
      `-`(1,.) %>%

      # return column distance matrix
      as.dist

    # end matrix columns distance calculation function definition
    }


# cluster columns of matrix
cluster.columns<-

  # define matrix column clustering function
  function(

    # matrix to cluster columns of
    mat

    # end matrix column clustering function parameter definition
    )

    # begin matrix column clustering function definition
    {

      # take matrix to cluster columns of
      mat %>%

      # calculate columns distance matrix
      get.column.distances %>%

      # (hirarchically) cluster columns based on their pairwise distances
      hclust

    # end matrix column clustering function definition
    }


# generate heatmap treating columns as rows
heatmap.cols_as_rows<-

  # define columns-as-rows heatmap function
  function(

    # expression matrix (genes as columns)
    expr.matrix

    # end columns-as-rows heatmap function parameter definition
    )

    # begin columns-as-rows heatmap function definition
    {

      # take expression matrix (genes as columns)
      expr.matrix %>%

      # z-transform columns (genes)
      scale %>%

      # use genes as rows
      t %>%

      # generate heatmap
      heatmaply(

        # don't change column order
        ,Colv=FALSE

        # set row dendrogram
        ,Rowv=

          # take expression matrix (genes as columns)
          expr.matrix %>%

          # cluster columns (genes)
          cluster.columns
      )

    # end columns-as-rows heatmap function definition
    }


# generate heatmap clustering rows
heatmap.rows<-

  # define row-clustered heatmap function
  function(

    # expression matrix (genes as rows)
    expression.matrix

    # end row-clustered heatmap function parameter definition
    )

    # begin row-clustered heatmap function definition
    {

      # take expression matrix (genes as rows)
      expression.matrix %>%

      # use genes as columns
      t %>%

      # generate heatmap treating columns as rows
      heatmap.cols_as_rows

    # end row-clustered heatmap function definition
    }


  ####################
  # server functions #
  ####################

    #############################
    # data extraction functions #
    #############################

# extract sample shift inputs from inputs list
get.sample.shifts<-

  # define sample shift input extraction function
  function(

    # sample names of samples to extract sample shift inputs for
    sample.names

    # input list to extract sample shift inputs from
    ,inputs

    # end sample shift input extraction function parameter definition
    )

    # begin sample shift input extraction function definition
    {

      # take sample names of samples to extract sample shift inputs for
      sample.names %>%

      # get corresponding sample shift input names
      get.sample.shift.inputId %>%

      # extract sample shift inputs from input list by input names
      llply(. %>% inputs[[.]]) %>%

      # label sample shift inputs with sample names
      setNames(sample.names)

    # end sample shift input extraction function definition
    }


    #########################
    # input panel functions #
    #########################

# sample shifts input panel generation
generate.sample.shifts.input<-

  # define sample shifts input panel generation function
  function(

    # sample names to include in sample shifts input panel
    sample.names

    # end sample shifts input panel generation function parameter definition
    )

    # begin sample shifts input panel generation function
    {

      # take sample names to include in sample shifts input panel
      sample.names %>%

      # generate sample shift input panel per sample
      llply(

        # generate sample shift input panel
        generate.sample.shift.input

        # end sample shift input panel generate per sample
        )

    # end sample shifts input panel generation function
    }


# generate genotype input panel
generate.genotype.input<-

  # define genotype input panel generation function
  function(

    # sample description to generate genotype input panel for
    sample.description

    # end genotype input panel generation function parameter definition
    )

    # begin genotype input panel generation function definition
    {

      # take sample description to generate genotype input panel for
      sample.description %>%

      # get corresponding available genotypes
      input.data$genotypes[[.]] %>%

      # generate genotype input panel
      selectInput(

        # name genotype input
        inputId="genotype"

        # label genotype input panel
        ,label=params$genotype.input.label %>%

          # make label 3rd level header
          h3

        # set choices for genotype input panel
        ,choices=.

        # set default selection for genotype input panel
        ,selected=params$genotype.input.default

        # end genotype input panel generation
        )

    # end genotype input panel generation function definition
    }


    ##################
    # plot functions #
    ##################

# profile plot generation function
generate.profile.plot<-

  # define filter by gene names function
  function(

    # tomoseq data to plot
    tomoseq.data

    # gene names of genes to include in plot
    ,gene.names

    # sample names of samples to include in plot
    ,sample.names=

      # derive default sample name to include in plot on given tomo seq data
      tomoseq.data %$%

      # identify all sample names provided data for
      unique(sample.name)

    # plot options to use
    ,plot.options=

      # by default, don't use any plot option
      character(0)

    # plot columns count to use
    ,ncols.plot=

      # by default, use default plot columns count
      params$ncols.plot.input.default

    # sample shift input list
    ,sample.shifts=

      # use default sample shift input values for all samples by default
      NULL

    # end profile plot generation parameter definition
    )

    # begin profile plot generation function definition
    {

      # take tomo-seq data
      tomoseq.data %>%

      # filter tomo-seq data by sample names
      filter.data.by.sample.names(

        # specify sample names of samples to keep
        sample.names.to.keep=

          # take sample names to include in plot
          sample.names %>%

          # drop last sample names if given sample number exceeds maximum
          head(params$profile.plot.max.nsamples)

        # end data filtering by sample names
        ) %>%

      # filter tomo-seq data by gene names
      filter.data.by.gene.names(

        # specify gene names of genes to keep
        gene.names.to.keep=

          # take gene names to include in plot
          gene.names %>%

          # convert single string input to vector
          parse.gene.names

        # end data filtering by gene names
        ) %>%

      # add shift column to data
      add.shift.column(

        # specify shifts to add by sample
        shifts.by.sample=

          # take sample shift
          sample.shifts %>%

          # convert named list (or NULL) to named vector
          parse.sample.shifts(

            # specify samples to convert sample shifts for
            samples=sample.names

            # end sample shift conversion
            )

        # end addition of shift column
        ) %>%

      # plot profiles based on filtered tomo-seq data
      plot.profiles(

        # show raw data points if specified in plot options
        raw.points="raw.points" %in% plot.options

        # show raw data lines if specified in plot options
        ,raw.lines="raw.lines" %in% plot.options

        # show per-sample smooth fits if specified in plot options
        ,smooth.each="smooth.each" %in% plot.options

        # show across-sample smooth fit if specified in plot options
        ,smooth.pooled="smooth.pooled" %in% plot.options

        # log-transform y-axis if specified in plot options
        ,logscale="logscale" %in% plot.options

        # fix x-axis limits if specified in plot options
        ,fix.xlim="fix.xlim" %in% plot.options

        # use a single y-scale for all sub-plots if specified in plot options
        ,single.y.scale="single.y.scale" %in% plot.options

        # use plot columns count specified
        ,ncols=ncols.plot

        # end profile plotting
        )

    # end profile plot generation function definition
    }


# heatmap generation function
generate.heatmap<-

  # define heatmap generation function
  function(

    # gene profiles to plot
    gene.profiles

    # sample description to generate heatmap for
    ,sample.description

    # genotype to generate heatmap for
    ,genotype

    # maximum number of genes to include in heatmap
    ,max.genes=

      # by default, use default maximum number of genes
      params$heatmap.max.ngenes

    # end heatmap generation function parameter definition
    )

    # begin heatmap generation function definition
    {

      # take gene profile
      gene.profiles %>%

      # extract gene profile for specified sample description
      `[[`(sample.description) %>%

      # extract gene profile for specified genotype
      `[[`(genotype) %>%

      # keep specified number of top varying genes
      keep.top.genes(max.genes) %>%

     heatmap.rows

    # end heatmap generation function definition
    }
