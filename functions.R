# tomo-seq shiny app function definition script
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

# file:         functions.R
# author(s):    Marcel Schilling <marcel.schilling@mdc-berlin.de>
# created:      2017-02-23
# last update:  2018-04-03
# license:      GNU Affero General Public License Version 3 (GNU AGPL v3)
# purpose:      define functions for tomo-seq shiny app


######################################
# change log (reverse chronological) #
######################################

# 2018-04-03: parameterized span to use for smoothing
#             parameterized number of points to impute for smoothing
# 2018-03-21: added support for CPM / cell abundance unit
# 2018-03-20: added support for slice width bars
#             added support for gonad arm model plot
# 2018-02-27: added usage of smooth fit span and number of points parameters
# 2018-01-05: switched from identifying transcripts by name to isoform number (for color assignment)
#             (identied as part after the last dot (".") in the transcript name; this enables
#             re-using the same color for different genes)
# 2017-10-23: switched from encoding samples by color and isoforms by shape/linetype to the opposite
#             added support for isoform-level profiles (using shape & linetype)
# 2017-10-20: added support for filtering by expression level type (gene/isoform profiles?) to
#             account for inclusion of isoform-specific CPM data (not used by now)
# 2017-05-29: added sample stretches input panel generation & input extraction functions / added
#             sample based stretch assignment / added stretch support to profile plot function
# 2017-05-24: added support for non-positive y-axis minumum for linearly scaled gene profiles
# 2017-05-23: added gene filtering by peak CPM minimum function
# 2017-05-22: added support for overwriting y-axis limits with user specified values
# 2017-05-17: made string-matching for normalization scheme determination fixed (i.e. non-regex)
#             replaced log.transform heatmap option by abundance measure (CPM, log2(1 + CPM) or
#             log10(1 + CPM))
#             replaced row.scaling heatmap option by row normalization scheme (scaling, centering or
#             none)
# 2017-04-19: added Euclidean distance as possible alternative to "1 - Pearson's r" metric for
#             gene clustering
# 2017-04-18: added gene type input panel generation and gene type based filtering functions
#             added removal of non-varying genes for heatmap
# 2017-04-13: added cpm.mean/min/max/lfc & percent.min/max to gene table annotation
#             added gene table annotation (gene.name & cpm.sd) support
# 2017-04-12: switched (back) from gene rank to count based filtering
#             moved gene rank based filtering out of heatmap function
#             moved genotype based filtering out of heatmap function
#             moved sample description based filtering out of heatmap function
#             moved gene list based filtering out of heatmap function
#             switched to tidy gene profile input data & rank based gene filtering
# 2017-04-11: added gene list file import functions
# 2017-04-10: added gene table XLSX export functions
# 2017-04-07: simplified gene cluster assignment extraction as suggested by Alicia Schep
# 2017-04-06: added gene table generation functions (incl. cluster assignment)
# 2017-04-05: added gene expression data log-transformation support
#             made row-scaling for heatmap optional
#             added gene cluster summary profiles
#             added gene cluster count support
#             fixed typo in default parameter definition
#             switched heatmap generation from heatmaply to iheatmapr
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

# get llply
require(plyr)

# get tibbles
require(tibble)

# get filter
require(dplyr)

# get geom_circle
require(ggforce)

# get plot_grid
require(cowplot)

# get acast
require(reshape2)

# get interactive complex heatmap framework
require(iheatmapr)

# get write.xlsx & read.xlsx
require(xlsx)


##############
# parameters #
##############

# load input data
source("data.R")

# load parameter definitions
source("params.R")


#############
# functions #
#############

  ####################
  # helper functions #
  ####################

# generate y-axis minimum input panel
generate.manual.ymin.input<-

  # define y-axis minimum input panel generation function
  function(

    # plot options selected by the user
    plot.options

    # end y-axis minimum input panel generation function parameter definition
    )

    # begin y-axis minimum input panel generation function definition
    {

      # assign input panel if manual y-axis limits plot option selected by the user
      if("set.ylim" %in% plot.options){

        # initiate minimal y-axis minimum
        min.value=params$manual.ymin.input.min

        # adjust minimal y-axis minimum based on logscale plot option
        if("logscale" %in% plot.options)

            # don't allow non-positive y-axis minimum for logarithmically scaled plots
            min.value %<>% max(1)

        # assign input panel
        manual.ymin.input<-

          # define y-axis minimum input panel
          numericInput(

            # name y-axis minimum input
            inputId="manual.ymin"

            # label y-axis minimum input panel
            ,label=params$manual.ymin.input.label %>%

              # make label 3rd level header
              h3

            # set minimal value for y-axis minimum input panel
            ,min=min.value

            # set maximal value for y-axis minimum input panel
            ,max=params$manual.ymin.input.max

            # set default value for y-axis minimum input panel to minimum
            ,value=min.value

            # end y-axis minimum input panel definition
            )

      # hide input panel if manual y-axis limits plot option not selected by the user
      } else {

        # assign placeholder
        manual.ymin.input<-

          # use empty HTML object as valid placeholder for input panel
          HTML("")

        # end input panel/placeholder assignment
        }

      # return input panel/placeholder for rendering
      return(manual.ymin.input)

    # end y-axis minimum input panel generation function definition
    }


# generate y-axis maximum input panel
generate.manual.ymax.input<-

  # define y-axis maximum input panel generation function
  function(

    # plot options selected by the user
    plot.options

    # manual y-axis limit set by the user
    ,ymin

    # end y-axis maximum input panel generation function parameter definition
    )

    # begin y-axis maximum input panel generation function definition
    {

      # assign input panel if manual y-axis limits plot option selected by the user
      if("set.ylim" %in% plot.options){

        # assign input panel
        manual.ymax.input<-

          # define y-axis maximum input panel
          numericInput(

            # name y-axis maximum input
            inputId="manual.ymax"

            # label y-axis maximum input panel
            ,label=params$manual.ymax.input.label %>%

              # make label 3rd level header
              h3

            # set minimal value for y-axis maximum input panel
            ,min=ymin+1

            # set maximal value for y-axis maximum input panel
            ,max=params$manual.ymax.input.max

            # set default value for y-axis maximum input panel to maximum
            ,value=params$manual.ymax.input.max

            # end y-axis maximum input panel definition
            )

      # hide input panel if manual y-axis limits plot option not selected by the user
      } else {

        # assign placeholder
        manual.ymax.input<-

          # use empty HTML object as valid placeholder for input panel
          HTML("")

        # end input panel/placeholder assignment
        }

      # return input panel/placeholder for rendering
      return(manual.ymax.input)

    # end y-axis maximum input panel generation function definition
    }

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

# get sample stretch input names by sample name
get.sample.stretch.inputId<-

  # define sample stretch input name function
  function(

    # sample name to get sample stretch input name for
    sample.name

    # end sample stretch input name function parameter definition
    )

    # begin sample stretch input name function definition
    {

      # take sample name to get sample stretch input name for
      sample.name %>%

      # add sample stretch input name prefix
      paste0("sample.stretch.",.)

    # end sample stretch input name function definition
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

# generate sample stretch input panel
generate.sample.stretch.input<-

  # define sample stretch input panel generation function
  function(

    # sample name of generate sample stretch input panel for
    sample.name

    # end sample stretch input panel generation function parameter definition
    )

    # begin sample stretch input panel generation function definition
    {

      # take sample name of generate sample stretch input panel for
      sample.name %>%

      # generate sample stretch input panel
      sliderInput(

        # name sample stretch input value
        inputId=get.sample.stretch.inputId(.)

        # label sample stretch input panel
        ,label=.

        # set sample stretch input panel minimum value
        ,min=params$sample.stretches.input.min

        # set sample stretch input panel maximum value
        ,max=params$sample.stretches.input.max

        # set sample stretch input panel default value
        ,value=params$sample.stretches.input.default

        # set sample stretch input panel value suffix
        ,post=params$sample.stretches.input.suffix

        # end sample stretch input panel generation
        )

    # end sample stretch input panel generation function definition
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


# filter data by expression level type (gene/isoform profiles?)
filter.data.by.expression.level<-

  # define filter by expression level type (gene/isoform profiles?) function
  function(

    # unfiltered data
    unfiltered.data

    # use isoform-level expression estimates?
    ,isoform.level

    # end filter by expression level type (gene/isoform profiles?) function parameter definition
    )

    # begin filter by expression level type (gene/isoform profiles?) function definition
    {

      # take unfiltered data
      unfiltered.data %>%

      # subset unfiltered data
      subset(

        # select data isoform-level data for isoform-level profiles
        is.na(transcript.name) != isoform.level

        # end data subsetting
        )

    # end filter by expression level type (gene/isoform profiles?) function definition
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

# convert named sample stretch list (or NULL) to named (default) sample stretch vector
parse.sample.stretches<-

  # define sample stretch conversion function
  function(

    # list with sample stretches labeled with sample names (or NULL)
    stretches

    # sample names to include in output sample stretch vector
    ,samples

    # end sample stretch conversion function parameter definition
    )

    # begin sample stretch conversion function definition
    {

      # if no sample stretches were specified, use default values
      if(is.null(stretches))

        # take default sample stretch input value
        params$sample.stretches.input.default %>%

        # repeat for each sample specified
        rep(times=length(samples)) %>%

        # label default sample stretches with sample names
        setNames(samples)

      # if sample stretches were specified, convert list to vector
      else

        # take sample names to include in output sample stretch vector
        samples %>%

        # extract corresponding sample stretches
        stretches[.] %>%

        # convert names list to vector
        unlist

    # end sample stretch conversion function definition
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

# add stretch column to data
add.stretch.column<-

  # define stretch column addition function
  function(

    # data to add stretch column to
    input.data

    # stretches to add to data labeled by sample name
    ,stretches.by.sample

    # end stretch column addition function parameter definition
    )

    # begin stretch column addition function definition
    {

      # take data to add stretch column to
      input.data %>%

      # add column to data
      cbind(

        # label stretch column
        stretch=

          # assign stretch by sample name
          stretches.by.sample[.$sample.name]

        # and column addition
        )

    # end stretch column addition function definition
    }

# generate gonad arm model plot
plot.model <- function(model.data)
  model.data$fit.radii %>%
  ggplot(aes(dist.to.dtc.um * 100 / model.data$len.gonad.arm.um,
             radius.gonad.um * 100 / model.data$len.gonad.arm.um)) %>%
  + geom_circle(data = model.data$germ.layers %>%
                       mutate(center.dist.um =
                                (start.dist.um + end.dist.um) / 2,
                              n.circles =
                                as.integer(diam.gonad.um /
                                             model.data$diam.germ.cell.um)) %>%
                       group_by(germ.layer) %>%
                       do(data_frame(center.dist.um = .$center.dist.um,
                                     n.circles = .$n.circles,
                                     circle = 1:.$n.circles)) %>%
                       ungroup %>%
                       mutate(y.center = (circle -.5 - n.circles / 2) *
                                           model.data$diam.germ.cell.um,
                              radius.um = model.data$diam.germ.cell.um / 2),
                aes(x0 = center.dist.um * 100 / model.data$len.gonad.arm.um,
                    y0 = y.center * 100 / model.data$len.gonad.arm.um,
                    r = radius.um * 100 / model.data$len.gonad.arm.um),
                inherit.aes = FALSE, color = "grey") %>%
  + geom_circle(data = model.data$oocytes %>%
                       mutate(center.dist.um =
                                (start.dist.um + end.dist.um) / 2,
                              radius.um = diam.um / 2),
                aes(x0 = center.dist.um * 100 / model.data$len.gonad.arm.um,
                    y0 = 0, r = radius.um * 100 / model.data$len.gonad.arm.um),
                inherit.aes = FALSE, color = "grey") %>%
  + geom_segment(data = model.data$zones %>%
                        filter(end.dist.um != model.data$len.gonad.arm.um),
                 aes(x = end.dist.um * 100 / model.data$len.gonad.arm.um,
                     xend = end.dist.um * 100 / model.data$len.gonad.arm.um,
                     y = -2500 / model.data$len.gonad.arm.um,
                     yend = 2500 / model.data$len.gonad.arm.um),
                 linetype = "dotted") %>%
  + geom_line() %>%
  + geom_line(data=model.data$fit.radii %>%
                     mutate(radius.gonad.um = -radius.gonad.um)) %>%
  + geom_segment(
      data = model.data$fit.radii[c(1, nrow(model.data$fit.radii)),],
      aes(xend = dist.to.dtc.um * 100 / model.data$len.gonad.arm.um,
          yend = -radius.gonad.um * 100 / model.data$len.gonad.arm.um)) %>%
  + annotate("segment",
             x = (model.data$len.gonad.arm.um - 60) * 100 /
                   model.data$len.gonad.arm.um,
             xend = (model.data$len.gonad.arm.um - 10) * 100 /
                      model.data$len.gonad.arm.um,
             y = -3000 / model.data$len.gonad.arm.um,
             yend = -3000 / model.data$len.gonad.arm.um, size = 1) %>%
  + annotate("text",
             x = (model.data$len.gonad.arm.um - 35) * 100 /
                   model.data$len.gonad.arm.um,
             y = -3500 / model.data$len.gonad.arm.um, vjust = 1,
             label = "50 μm") %>%
  + expand_limits(y = -5000 / model.data$len.gonad.arm.um) %>%
  + coord_fixed() %>%
  + theme_void()

# globally define gonad arm model plot
model.plot <- plot.model(input.data$gonad.model)

dist2n.cells <- function(end.dist.um, start.dist.um = 0,
                         cell.data = input.data$gonad.model$cells){

  # keep only relevant cell data
  cell.data %<>%
    filter(center.dist.um + radius.um >= start.dist.um &
           center.dist.um - radius.um <= end.dist.um)

  # catch zero-cells case
  if(!nrow(cell.data)) return(0)

  # weight cells by contained fraction of their total volume
  cell.data %>%
    mutate(width.um = pmin(center.dist.um + radius.um, end.dist.um) -
                      pmax(center.dist.um - radius.um, start.dist.um),
           vol.um3 = (pi * width.um^2) / 3 * (3 * radius.um - width.um),
           fraction = vol.um3 / (4 / 3 * pi * radius.um^3)) %$%
    sum(n.cells * fraction)
}

unit2col <- function(unit.name)
  ifelse(unit.name == "CPM / cell", "cpm.per.cell", "cpm")

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

    # set y-axis limits automatically by default
    ,y.limits=NULL

    # use default plot columns count by default
    ,ncols=params$ncols.plot.input.default

    # plot gene level estimates by default
    ,isoform.level=F,

    show.slice.width = TRUE,

    abundance.unit = params$abundance.unit.default,
    n.points.smooth = params$smoothing.n.input.default,
    span.smooth = params$smoothing.span.input.default,
    show.model = TRUE)

    # begin profile plot function definition
    {

      # create profile plot
      profile.plot<-

        # take plot data
        plot.data %>%
        mutate(percent.center = percent.center * stretch + shift,
               percent.start = percent.center - (width.percent / 2 * stretch),
               percent.end = percent.center + (width.percent / 2 * stretch),
               dist.start.um = percent.start / 100 *
                                 input.data$gonad.model$len.gonad.arm.um,
               dist.end.um = percent.end / 100 *
                               input.data$gonad.model$len.gonad.arm.um,
               tx_color = isoform.level %>%
                          ifelse(list(transcript.name %>%
                                      sub(".*[.]", "isoform ", .)
                                     ),
                                 list("black")) %>%
                          unlist) %>%
        group_by(dist.start.um, dist.end.um) %>%
        mutate(n.cells = dist2n.cells(dist.end.um[1], dist.start.um[1])) %>%
        ungroup %>%
        mutate(cpm.per.cell = cpm / n.cells) %>%
        ggplot(aes_string(x = "percent.center", y = unit2col(abundance.unit),
                          color = "tx_color", linetype = "sample.name")) %>%

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
        + ylab(label = params$profile.plot.ylab %>%
                       paste(paste0("[", abundance.unit, "]"),
                             sep = ifelse((nchar(abundance.unit) > 4) &
                                            show.model, "\n", " "))) %>%

        # adjust linetype (i.e. sample name) legend
        + scale_linetype_discrete(

            # adjust linetype legend label
            name=params$profile.plot.sample.legend.label

            # end linetype legend parameter definition
            )

        # use color for isoform-level profiles
        if(isoform.level)

          # modify profile plot
          profile.plot %<>%

          # adjust color (i.e. isoform) palette/legend
          + scale_color_brewer(

              # adjust sample name -> color mapping
              palette=params$profile.plot.brewer.palette

              # adjust color legend label
              ,name=params$profile.plot.isoform.legend.label

              # end color palette/legend parameter definition
              )

        # keep gene-level profiles black & white
        else

          # modify profile plot
          profile.plot %<>%

          # adjust color (i.e. isoform) palette/legend
          + scale_color_grey(

              # hide color legend
              guide="none"

              # end color palette/legend parameter definition
              )

        if(show.slice.width)
          profile.plot %<>%
          + geom_errorbarh(aes(xmin = percent.start, xmax = percent.end))

        # add raw data points if specified
        if(raw.points)

          # modify profile plot
          profile.plot %<>%

          # plot raw data as points
          + geom_point(

              # group data by sample name & shape points accordingly
              aes(shape=sample.name)

              # adjust data point size
              ,size=params$profile.plot.pointsize

              # end data point parameter definition
              ) %>%

          # adjust shape (i.e. sample name) legend
          + scale_shape_discrete(

              # adjust shape legend label
              name=params$profile.plot.sample.legend.label

              # end shape legend parameter definition
              )

        # add raw data lines if specified
        if(raw.lines)

          # modify profile plot
          profile.plot %<>%

          # plot per-slice data as lines
          + geom_line(

              # adjust data line size
              size=params$profile.plot.linesize.each

              # end data line parameter definition
              )

        # add per-sample smooth fits if specified
        if(smooth.each)

          # modify profile plot
          profile.plot %<>%

          # fit smooth line across samples
          + geom_smooth(

              # adjust smoothing method
              ,method = params$profile.plot.smoothing.method,

              span = span.smooth,
              n = n.points.smooth,
              size = params$profile.plot.linesize.each

              # end smooth line parameter definition
              )

        # add across-sample smooth fit if specified
        if(smooth.pooled)

          # modify profile plot
          profile.plot %<>%

          # fit smooth line across samples
          + geom_smooth(

              # pool data across samples & adjust smooth line color
              linetype=params$profile.plot.linetype.smooth

              # adjust smoothing method
              ,method=params$profile.plot.smoothing.method,

              span = span.smooth,
              n = n.points.smooth,
              size = params$profile.plot.linesize.pooled

              # end smooth line parameter definition
              )

        # scale y-axis logarithmically if specified
        if(logscale){

          # only log2-transform y-axis if no y-axis limits specified
          if(is.null(y.limits)) {

            # modify profile plot
            profile.plot %<>%

            # log2-transform y-axis
            + scale_y_continuous(trans=log2_trans())

          # log2-transform y-axis & overwrite y-axis limits if specified
          } else {

            # modify profile plot
            profile.plot %<>%

            # adjust y-axis
            + scale_y_continuous(

                # log2-transform y-axis
                trans=log2_trans()

                # overwrite y-axis limits
                ,limits=y.limits

                # end y-axis adjustment
                )

          # end logarithmic y-axis adjustment
          }

        # adjust linear y-axis if specified
        } else {

          # overwrite y-axis limits if specified
          if(!is.null(y.limits)) {

            # modify profile plot
            profile.plot %<>%

            # overwrite y-axis limits
            + ylim(y.limits)

          # end overwriting of y-axis limits
          }

        # end linear y-axis adjustment
        }

        # fix x-axis limits if specified
        if(fix.xlim)

          # modify profile plot
          profile.plot %<>%

          # fix x-axis limits
          + xlim(params$profile.plot.xlim)

        if(show.model)
          profile.plot %<>%
            plot_grid(model.plot, ncol = 1, align = "v", axis = "lr")

        # return profile plot
        profile.plot

    # end profile plot function definition
    }


# extract column from table file by column index
get.column.from.file<-

  # define column from file extraction function
  function(

    # single row data.frame with name of table file to extract column from in datapath column
    input.file

    # index of table column to extract
    ,column.index

    # index of table sheet to extract column from
    ,sheet.index=1

    # end column from file extraction function parameter definition
    )

    # begin column from file extraction function definition
    {

      # take file to extract column from in datapath column
      input.file %$%

      # read column from XLSX file
      read.xlsx(

        # get input file path from input data.frame
        file=datapath

        # read only specified table sheet
        ,sheetIndex=sheet.index

        # read only specified table column
        ,colIndex=column.index

        # end XLSX file reading
        ) %>%

      # convert single column data.frame to vector
      unlist %>%

      # strip off rownames before returning column data
      unname

    # end column from file extraction function definition
    }


# extract gene list from table file
get.genes.from.file<-

  # define gene list from file extraction function
  function(

    # single row data.frame with name of table file to extract gene list from in datapath column
    table.file

    # index of gene list column within table file
    ,gene.column=1

    # end gene list from file extraction function parameter definition
    )

    # begin gene list from file extraction function definition
    {

      # take table file to extract gene list from
      table.file %>%

      # extract gene list column from table by specified column index
      get.column.from.file(gene.column) %>%

      # convert gene list to characters before returning them
      as.character

    # end gene list from file extraction function definition
    }


# extract CPM matrix (genes as rows)
get.cpm.matrix<-

  # define CPM matrix extraction function
  function(

    # expression data (gene, percent & cpm.fit columns)
    expression.data

    # name of CPM column
    ,cpm.colname="cpm.fit"

    # end CPM matrix extraction function parameter definition
    )

    # begin CPM matrix extraction function definition
    {

      # take expression data
      expression.data %>%

      # convert tidy table to matrix
      acast(

        # use genes as rows, percent as columns
        gene~percent

        # fill matrix with corresponding CPM values
        ,value.var=cpm.colname

        # end tidy table to matrix conversion
        )

    # end CPM matrix extraction function definition
    }


# calculate distance matrix of matrix columns using "1 - Pearson's r" metric
get.column.distances.pearson<-

  # define matrix columns "1 - Pearson's r" distance calculation function
  function(

    # matrix to get column distances of
    column.matrix

    # end matrix columns "1 - Pearson's r" distance calculation function parameter definition
    )

    # begin matrix columns "1 - Pearson's r" distance calculation function definition
    {

      # take matrix to get column distances of
      column.matrix %>%

      # calculate (pairwise) column (Pearson) correlations
      cor(use="complete") %>%

      # convert column similarities to distances
      `-`(1,.) %>%

      # return column distance matrix
      as.dist

    # end matrix columns "1 - Pearson's r" distance calculation function definition
    }


# normalize matrix columns
normalize.columns<-

  # define matrix column normalization function
  function(

    # matrix to normalize columns of
    column.matrix

    # should columns be z-transformed (FALSE: centering only)
    ,z.transform = TRUE

    # end matrix column normalization function parameter definition
    )

    # begin matrix column normalization function definition
    {

      # take matrix to normalize columns of
      column.matrix %>%

      # scale columns (z-transform if specified, otherwise center only)
      scale(scale = z.transform)

    # end matrix column normalization function definition
    }


# normalize matrix rows
normalize.rows<-

  # define matrix row normalization function
  function(

    # matrix to normalize rows of
    row.matrix

    # should rows be z-transformed (FALSE: centering only)
    ,scaling = TRUE

    # end matrix row normalization function parameter definition
    )

    # begin matrix row normalization function definition
    {

      # take matrix to normalize rows of
      row.matrix %>%

      # treat rows as columns
      t %>%

      # normalize columns (scale if specified)
      normalize.columns(z.transform = scaling) %>%

      # treat columns as rows
      t

    # end matrix row normalization function definition
    }


# calculate distance matrix of matrix rows using "1 - Pearson's r" metric
get.row.distances.pearson<-

  # define matrix rows "1 - Pearson's r" distance calculation function
  function(

    # matrix to get row distances of
    row.matrix

    # end matrix rows "1 - Pearson's r" distance calculation function parameter definition
    )

    # begin matrix rows "1 - Pearson's r" distance calculation function definition
    {

      # take matrix to get row distances of
      row.matrix %>%

      # treat rows as columns
      t %>%

      # calculate column "1 - Pearson's r" distances
      get.column.distances.pearson

    # end matrix rows "1 - Pearson's r" distance calculation function definition
    }


# calculate distance matrix of matrix rows using "Euclidean distance" metric
get.row.distances.euclidean<-

  # define matrix rows Euclidean distance calculation function
  function(

    # matrix to get row distances of
    row.matrix

    # end matrix rows Euclidean distance calculation function parameter definition
    )

    # begin matrix rows Euclidean distance calculation function definition
    {

      # take matrix to get row distances of
      row.matrix %>%

      # calculate Euclidean distances of rows
      dist(method="euclidean")

    # end matrix rows Euclidean distance calculation function definition
    }


# calculate distance matrix of matrix rows
get.row.distances<-

  # define named list of matrix rows distance calculation functions
  list(

    # add matrix rows "1 - Pearson's r" distance calculation function
    `1 - Pearson's r`=get.row.distances.pearson

    # add matrix rows Euclidean distance calculation function
    ,`Euclidean distance`=get.row.distances.euclidean

    # end definition of  named list of matrix rows distance calculation functions
    )


# log2-transform abundance estimates
transformation.log2p1<-

  # define log2-transformation function
  function(

    # abundance estimates to log2-transform
    abundance

    # pseudocount to add to abundance values before log2-transformation
    ,pseudocount=1

    # end log2-transformation function parameter definition
    )

    # end log2-transformation function definition
    {

      # take abundance estimates to log2-transform
      abundance %>%

      # add pseudocount
      `+`(pseudocount) %>%

      # log2-transform abundance estimates
      log2

    # end log2-transformation function definition
    }


# log10-transform abundance estimates
transformation.log10p1<-

  # define log10-transformation function
  function(

    # abundance estimates to log10-transform
    abundance

    # pseudocount to add to abundance values before log10-transformation
    ,pseudocount=1

    # end log10-transformation function parameter definition
    )

    # end log10-transformation function definition
    {

      # take abundance estimates to log10-transform
      abundance %>%

      # add pseudocount
      `+`(pseudocount) %>%

      # log10-transform abundance estimates
      log10

    # end log10-transformation function definition
    }


# get transformation to apply to CPM abundances to obtain specified measure
get.transformation<-

  # define transformation assignment function
  function(

    # abundance measure to transform CPM to
    measure=

      # by default, use default abundance measure scheme
      params$abundance.measure.input.default

    # end transformation assignment function parameter definition
    )

    # begin transformation assignment function definition
    {

      # check if log2 transformed CPM were requested
      if(grepl("log2(",measure,fixed=T))

        # return log2 transform function
        return(transformation.log2p1)

      # check if log10 transformed CPM were requested
      if(grepl("log10(",measure,fixed=T))

        # return log10 transform function
        return(transformation.log10p1)

      # if no transformation was requested, 'transform' CPM by identity function
      return(identity)

    # end transformation assignment function definition
    }


# generate heatmap clustering rows
heatmap.rows<-

  # define row-clustered heatmap function
  function(

    # expression matrix (genes as columns)
    expression.matrix

    # transformation to apply to abundance
    ,transformation=

      # by default, use default abundance measure scheme
      params$abundance.measure.input.default %>%

      # get corresponding transformation function
      get.transformation

    # number of clusters to cluster rows into
    ,nclust=

      # by default, use default number of gene clusters
      params$nclust.genes.input.default

    # normalization scheme to use for rows
    ,row.norm=

      # by default, use default row normalization scheme
      params$row.normalization.input.default

    # distance metric to cluster rows by
    ,dist.metric=

      # by default, use default distance metric
      params$distance.metric.input.default

    # color values to use for heatmap tiles
    ,color.values=

      # by default, use default color values
      params$heatmap.color.values

    # side of heatmap to place dendrogram at
    ,dendrogram.side=

      # by default, place dendrogram at default side of heatmap
      params$heatmap.dendrogram.side

    # plotly scatter mode to use for row cluster summary
    ,summary.mode=

      # by default, use default mode
      params$heatmap.summary.mode

    # line type to use for row cluster summary
    ,summary.line=

      # by default, use list of default plotly line attributes
      # see https://plot.ly/javascript/reference/#scatter-line for options
      list(

        # by default, modify line with to non-plotly default
        width=

          # by default, use default summary line width
          params$heatmap.summary.linewidth

        # end default plotly line attributes list definition for summary line
        )

    # font to use for x-axis labels
    ,xlab.font=

      # by default, use list of default plotly font attributes
      # see https://plot.ly/javascript/reference/#layout-font for options
      list(

        # by default, modify font size to non-plotly default
        size=

          # by default, use default x-axis label font size
          params$heatmap.xlab.fontsize

        # end default plotly font attributes list definition for x-axis labels
        )

    # font to use for y-axis labels
    ,ylab.font=

      # by default, use list of default plotly font attributes
      # see https://plot.ly/javascript/reference/#layout-font for options
      list(

        # by default, modify font size to non-plotly default
        size=

          # by default, use default y-axis label font size
          params$heatmap.ylab.fontsize

        # end default plotly font attributes list definition for x-axis labels
        )

    # end row-clustered heatmap function parameter definition
    )

    # begin row-clustered heatmap function definition
    {

      # take expression matrix
      expression.matrix %<>%

      # apply specified transformation
      transformation

      # normalize rows if specified
      if(row.norm != "none")

        # take expression matrix
        expression.matrix %<>%

        # normalize rows (i.e. genes)
        normalize.rows(

          # determine if scaling should be performed (or centering only)
          scaling =

            # take row normalization scheme
            row.norm %>%

            # search for "scaling" substring
            grepl("scaling",.,fixed=T)

          # end row normalization
          )

      # take expression matrix (genes as columns)
      expression.matrix %>%

      # generate heatmap
      main_heatmap(

        # set color values to use for heatmap tiles
        colors=

          # use specified color values
          color.values

        # end main heatmap definition
        ) %>%

      # cluster heatmap rows (i.e. genes)
      add_row_clustering(

        # cluster rows (i.e. genes) into as many clusters as specified
        k=nclust

        # set distance function for row clustering
        ,clust_dist=

          # use matrix rows distance calculation function
          get.row.distances[[dist.metric]]

        # set side of heatmap to place row dendrogram
        ,side=

          # place row dendrogram on the specified side
          dendrogram.side

        # end heatmap row clustering parameter definition
        ) %>%

      # add summary profiles
      add_col_summary(

        # show summary profile per row group
        groups=

          # group rows by cluster assignment
          TRUE

        # set summary mode
        ,mode=

          # use specified summary mode
          summary.mode

        # set line type
        ,line=

          # use specified summary line type
          summary.line

        # end summary profiles addition
        ) %>%

      # label heatmap columns
      add_col_labels(

        # set column label font
        font=

          # use specified x-axis label font
          xlab.font

        # end heatmap column label parameter definition
        ) %>%

      # label heatmap rows
      add_row_labels(

        # set row label font
        font=

          # use specified y-axis label font
          ylab.font

        # end heatmap row label parameter definition
        )

    # end row-clustered heatmap function definition
    }


# extract row cluster assignment from heatmap
# see https://github.com/AliciaSchep/iheatmapr/issues/2#issuecomment-292393474
get.row.clusters<-

  # define row cluster assignment extraction function
  function(

    # heatmap to extract row cluster assignment from
    heatmap

    # end row cluster assignment extraction function parameter definition
    )

    # begin row cluster assignment extraction function definition
    {

      # take heatmap to extract row cluster assignment from
      heatmap %>%

      # extract heatmap plot list
      plots %>%

      # extract heatmap row cluster plot
      `[[`("Row<br>Clusters") %>%

      # extract heatmap row cluster data
      (iheatmapr:::get_data)

    # end row cluster assignment extraction function definition
    }


# save table to XLSX file of specified name
save.xlsx<-

  # define XLSX export function
  function(

    # table to save to XLSX
    input.table

    # output file name
    ,output.file

    # other parameters to be passed on to write.xlsx
    ,...

    # end XLSX export function parameter definition
    )

    # begin XLSX export function definition
    {

      # define temporary output file name
      # see http://stackoverflow.com/a/21388005/2451238
      tmp.xlsx<-

        # take specified output file name
        output.file %>%

        # append XLSX extension
        paste0(".xlsx")

      # take table to save to XLSX
      input.table %>%

      # convert to data.frame for write.xlsx
      as.data.frame %>%

      # write data.frame to XLSX file
      write.xlsx(

        # set output file name
        file=

          # use temporary output file name
          tmp.xlsx

        # pass on other parameters
        ,...

        # end writign of data.frame to XLSX file
        )

      # take temporary output XLSX file name
      tmp.xlsx %>%

      # rename temporary XLSX output file to specified output file name
      file.rename(output.file)

    # end XLSX export function definition
    }


  ####################
  # server functions #
  ####################


    ##################################
    # parameter extraction functions #
    ##################################

# extract file name to use for gene table XLSX export
get.gene.table.xlsx.name<-

  # define gene table XLSX file name export function
  function(

    # This function doesn't have any parameters as it is returning a constant.

    # end gene table XLSX file name export function parameter definition
    )

    # begin gene table XLSX file name export function definition
    {

      # return gene table XLSX file name parameter
      params$gene.table.xlsx.name

    # end gene table XLSX file name export function definition
    }


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

# extract sample stretch inputs from inputs list
get.sample.stretches<-

  # define sample stretch input extraction function
  function(

    # sample names of samples to extract sample stretch inputs for
    sample.names

    # input list to extract sample stretch inputs from
    ,inputs

    # end sample stretch input extraction function parameter definition
    )

    # begin sample stretch input extraction function definition
    {

      # take sample names of samples to extract sample stretch inputs for
      sample.names %>%

      # get corresponding sample stretch input names
      get.sample.stretch.inputId %>%

      # extract sample stretch inputs from input list by input names
      llply(. %>% inputs[[.]]) %>%

      # label sample stretch inputs with sample names
      setNames(sample.names)

    # end sample stretch input extraction function definition
    }


# extract gene annotation from gene profiles table
get.gene.annotation<-

  # define gene annotation extract function
  function(

    # gene profiles table to extract annotation from
    gene.profiles

    # end gene annotation extract function parameter definition
    )

    # begin gene annotation extract function definition
    {

      # take gene profiles table to extract annotation from
      gene.profiles %>%

      # extract gene annotation columns
      distinct(

        # get gene ID
        gene

        # add gene name
        ,gene.name

        # add gene CPM mean
        ,cpm.mean

        # add gene CPM standard deviation
        ,cpm.sd

        # add gene CPM minimum
        ,cpm.min

        # add gene CPM maximum
        ,cpm.max

        # add gene CPM minimum to maximum log(2) fold-change
        ,cpm.lfc

        # add (median) position (percent distal-to-proximal) where minumum gene CPM was observed
        ,percent.min

        # add (median) position (percent distal-to-proximal) where maximum gene CPM was observed
        ,percent.max

        # extract gene annotation columns
        )

    # end gene annotation extract function definition
    }


    ############################
    # data filtering functions #
    ############################

# only keep matrix rows specified in input file (if any)
filter.data.by.genes.file<-

  # define matrix filtering by gene list file function
  function(

    # data to subset by gene list file
    unfiltered.data

    # single row data.frame with gene list file name in datapath column (or NULL)
    ,genes.file

    # end matrix filtering by gene list file function parameter definition
    )

    # begin matrix filtering by gene list file function definition
    {

      # if no gene list file was specified, return unfiltered data
      if(is.null(genes.file)) unfiltered.data

      # filter data based on gene list if applicable
      else

        # take list file file
        genes.file %>%

        # extract gene list from file
        get.genes.from.file %>%

        # subset data based on gene list
        {filter(unfiltered.data,gene %in% .)}

    # end data filtering by gene list file function definition
    }


# filter data by sample description
filter.data.by.sample.description<-

  # define filter by sample description function
  function(

    # unfiltered data
    unfiltered.data

    # sample description to keep
    ,sample.description.to.keep

    # end filter by sample description function parameter definition
    )

    # begin filter by sample description function definition
    {

      # take unfiltered data
      unfiltered.data %>%

      # subset unfiltered data
      filter(

        # select data with sample description to keep
        sample.description==sample.description.to.keep

        # end data subsetting
        )

    # end filter by sample description function definition
    }


# filter data by genotype
filter.data.by.genotype<-

  # define filter by genotype function
  function(

    # unfiltered data
    unfiltered.data

    # genotype to keep
    ,genotype.to.keep

    # end filter by genotype function parameter definition
    )

    # begin filter by genotype function definition
    {

      # take unfiltered data
      unfiltered.data %>%

      # subset unfiltered data
      filter(

        # select data with genotype to keep
        genotype==genotype.to.keep

        # end data subsetting
        )

    # end filter by genotype function definition
    }


# filter data by gene type
filter.data.by.gene.type<-

  # define filter by gene type function
  function(

    # unfiltered data
    unfiltered.data

    # gene type to keep
    ,gene.type.to.keep

    # end filter by gene type function parameter definition
    )

    # begin filter by gene type function definition
    {

      # take unfiltered data
      unfiltered.data %>%

      # subset unfiltered data
      filter(

        # select data with gene type to keep
        gene.type==gene.type.to.keep

        # end data subsetting
        )

    # end filter by gene type function definition
    }


# filter data by minimum peak CPM
filter.data.by.min.cpm.max<-

  # define filter by minimum peak CPM function
  function(

    # unfiltered data
    unfiltered.data

    # peak CPM lower threshold
    ,min.cpm.max

    # end filter by minimum peak CPM function parameter definition
    )

    # begin filter by minimum peak CPM function definition
    {

      # take unfiltered data
      unfiltered.data %>%

      # subset unfiltered data
      filter(

        # select data exceeding peak CPM lower threshold
        cpm.max>=min.cpm.max

        # end data subsetting
        )

    # end filter by minimum peak CPM function definition
    }


# only top variable genes
keep.top.genes<-

  # define top variable genes filter function
  function(

    # expression matrix (genes as rows)
    unfiltered.data

    # maximum number of genes to keep
    ,nmax=

      #by default, keep default maximum number of genes
      params$heatmap.nmax.genes

    # end top variable genes filter function parameter definition
    )

    # begin top variable genes filter function definition
    {

      # take unfiltered data
      unfiltered.data %>%

      # isolate gene/standard-deviation pairs
      distinct(gene,cpm.sd) %>%

      # drop non-varying genes
      filter(cpm.sd>0) %>%

      # sort genes by standard deviation
      arrange(desc(cpm.sd)) %>%

      # keep up to specified maximum number of highest varying genes
      head(nmax) %>%

      # extract gene profiles of top varying genes
      left_join(unfiltered.data)

    # end top variable genes filter function definition
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


# sample stretches input panel generation
generate.sample.stretches.input<-

  # define sample stretches input panel generation function
  function(

    # sample names to include in sample stretches input panel
    sample.names

    # end sample stretches input panel generation function parameter definition
    )

    # begin sample stretches input panel generation function
    {

      # take sample names to include in sample stretches input panel
      sample.names %>%

      # generate sample stretch input panel per sample
      llply(

        # generate sample stretch input panel
        generate.sample.stretch.input

        # end sample stretch input panel generate per sample
        )

    # end sample stretches input panel generation function
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


# generate gene type input panel
generate.gene.type.input<-

  # define gene type input panel generation function
  function(

    # sample description to generate gene type input panel for
    sample.description

    # genotype to generate gene type input panel for
    ,genotype

    # end gene type input panel generation function parameter definition
    )

    # begin gene type input panel generation function definition
    {

      # take sample description to generate gene type input panel for
      sample.description %>%

      # get corresponding available gene types
      input.data$gene.types[[.]] %>%

      # get gene types for genotype to generate gene type input panel for
      `[[`(genotype) %>%

      # generate gene type input panel
      selectInput(

        # name gene type input
        inputId="gene.type"

        # label gene type input panel
        ,label=params$gene.type.input.label %>%

          # make label 3rd level header
          h3

        # set choices for gene type input panel
        ,choices=.

        # set default selection for gene type input panel
        ,selected=params$gene.type.input.default

        # end gene type input panel generation
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

    # manual y-axis limits specified by the user
    ,manual.ylim=

      # by default, set y-axis limits automatically
      NULL

    # plot columns count to use
    ,ncols.plot=

      # by default, use default plot columns count
      params$ncols.plot.input.default

    # sample shift input list
    ,sample.shifts=

      # use default sample shift input values for all samples by default
      NULL

    # sample stretch input list
    ,sample.stretches=

      # use default sample stretch input values for all samples by default
      NULL

    # plot isoform-specific profiles?
    ,per.isoform=

      # by default, plot gene level profiles
      FALSE,

    unit = params$abundance.unit.default,
    smoothing.n = params$smoothing.n.input.default,
    smoothing.span = params$smoothing.span.input.default)

    # begin profile plot generation function definition
    {

      # if manual y-axis limits plot option not selected by the user
      if(!("set.ylim" %in% plot.options))

        # discard any manual y-axis limits specified
        manual.ylim<-NULL

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

      # filter tomo-seq data by expression level type (gene/isoform profiles?)
      filter.data.by.expression.level(

        # specify wether to use isoform-level expression estimates
        isoform.level=

          # take expression level type to use for plot
          per.isoform

        # end data filtering by expression level type (gene/isoform profiles?)
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

      # add stretch column to data
      add.stretch.column(

        # specify stretches to add by sample
        stretches.by.sample=

          # take sample stretch
          sample.stretches %>%

          # convert named list (or NULL) to named vector
          parse.sample.stretches(

            # specify samples to convert sample stretches for
            samples=sample.names

            # end sample stretch conversion
            )

        # end addition of stretch column
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

        # use y-axis limits specified
        ,y.limits=manual.ylim

        # use plot columns count specified
        ,ncols=ncols.plot

        # plot isoform-specific profiles if specified
        ,isoform.level=per.isoform,

        abundance.unit = unit,
        n.points.smooth = smoothing.n,
        span.smooth = smoothing.span,
        show.slice.width = "show.slice.width" %in% plot.options,

        show.model =
          ("show.model" %in% plot.options) & ("fix.xlim" %in% plot.options) &
            ((ncols.plot == 1) | (length(parse.gene.names(gene.names)) == 1))

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

    # number of clusters to cluster genes into
    ,nclust.genes=

      # by default, use default number of gene clusters
      params$nclust.genes.input.default

    # abundance measure to use for gene profiles
    ,abundance.measure=

      # by default, use default normalization scheme
      params$abundance.measure.input.default

    # normalization scheme to use for gene profiles
    ,row.normalization=

      # by default, use default normalization scheme
      params$row.normalization.input.default

    # distance metric to use for clustering
    ,distance.metric=

      # by default, use default distance metric for clustering
      params$distance.metric.input.default

    # end heatmap generation function parameter definition
    )

    # begin heatmap generation function definition
    {

      # take gene profile
      gene.profiles %>%

      # extract CPM matrix (genes as rows)
      get.cpm.matrix %>%

      # generate heatmap clustering rows (i.e. genes)
      heatmap.rows(

        # transform abundance to obtain the measure specified
        transformation=get.transformation(abundance.measure)

        # cluster genes into as many clusters as specified
        ,nclust=nclust.genes

        # normalize gene profiles using the scheme specified
        ,row.norm=row.normalization

        # cluster genes according to distance metric specified
        ,dist.metric=distance.metric

        # end row-clustered heatmap generation
        )

    # end heatmap generation function definition
    }


    ##########################
    # table output functions #
    ##########################

# generate gene table
generate.gene.table<-

  # define gene table generation function
  function(

    # gene clustered heatmap to generate gene table for
    gene.heatmap

    # gene annotation used to annotate gene table
    ,annotation=

      # by default, don't add any annotation
      data_frame(gene=character(0))

    # end gene table generation function parameter definition
    )

    # begin gene table generation function definition
    {

      # take gene clustered heatmap to generate gene table for
      gene.heatmap %>%

      # get row (i.e gene) cluster assignment
      get.row.clusters %>%

      # convert named vector to tibble
      data_frame(gene=names(.),cluster=.) %>%

      # merge gene table with given annotation
      left_join(annotation)

    # end gene table generation function definition
    }


# generate gene table option list
generate.gene.table.options<-

  # define gene table option list generation function
  function(

    # (default) number of genes to show per page
    ngenes=

      # by default, use default (default) number of genes to show per page
      params$gene.table.ngenes

    # end gene table option list generation function parameter definition
    )
    # begin gene table option list generation function definition
    {

      # define table output option list
      # see https://datatables.net/reference/option/ for options
      list(

        # set (initial) number of rows per page
        pageLength=

          # use specified (default) number of genes to show per page
          ngenes

        # end  table output option list definition
        )

    # end gene table option list generation function definition
    }


    #########################
    # data export functions #
    #########################

# export gene table to XLSX
save.gene.table.xlsx<-

  # define gene table XLSX export function
  function(

    # gene table to export to XLSX
    gene.table

    # file name to save XLSX ouput to
    ,output.xlsx

    # end gene table XLSX export function parameter definition
    )

    # begin gene table XLSX export function definition
    {

      # take gene table to export to XLSX
      gene.table %>%

      # export gene table to XLSX
      save.xlsx(

        # set output file name for XLSX export
        output.file=

          # use specified XLSX output file name
          output.xlsx

        # don't include row names as a column in output XLSX file
        # (gene table should be in tidy format)
        ,row.names=FALSE

        # end gene table XLSX export
        )

    # end gene table XLSX export function definition
    }
