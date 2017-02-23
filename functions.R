#######################
# general information #
#######################

# file:         functions.R
# author(s):    Marcel Schilling <marcel.schilling@mdc-berlin.de>
# created:      2017-02-23
# last update:  2017-02-23
# purpose:      define functions for tomo-seq shiny app


######################################
# change log (reverse chronological) #
######################################

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

            # position points on the x-axis according to the center of the corresponding slice
            # (percent distal-to-proximal)
            x=percent.center

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

            # adjust the (max.) number of sub-plots to put underneath each other
            ,nrow=params$profile.plot.nrow

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

        # return profile plot
        profile.plot

    # end profile plot function definition
    }


  ####################
  # server functions #
  ####################

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

        # end profile plotting
        )

    # end profile plot generation function definition
    }
