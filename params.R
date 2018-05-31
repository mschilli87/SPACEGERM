# SPACEGERM shiny app parameter definition script
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

# file:         params.R
# author(s):    Marcel Schilling <marcel.schilling@mdc-berlin.de>
# created:      2017-02-21
# last update:  2018-05-31
# license:      GNU Affero General Public License Version 3 (GNU AGPL v3)
# purpose:      define parameters for SPACEGERM shiny app


######################################
# change log (reverse chronological) #
######################################

# 2018-05-31: increased smooth fit span default to match publication
#             removed slice data RDS (replaced by SQLite database)
# 2018-05-30: replaced shift/stretch RDS input by SQLite database
# 2018-05-17: adjusted tab titles for publication
#             added app subtitle
#             replaced require by library
# 2018-05-16: renamed app for publication
# 2018-04-23: added 3D model fixed expression range plot option
# 2018-04-16: renamed y-axis limits inputs to expression range inputs
# 2018-04-13: added 3D model gene name input default parameter
#             added 3D model gene name input panel label parameter
#             removed un-needed parameters left over from sample shift/stretch input panels
#             added default ymin/max values
#             adjusted sample description and genotype defaults
#             added 3D model colorbar parameters
#             added 3D model tab & plot parameters
# 2018-04-12: disabled single-y-scale by default
# 2018-04-10: added plot option to show/hide dropout slices
#             added plot option to show/hide smoothing standard error
# 2018-04-09: disabled gonad arm model by default
# 2018-04-05: adjusted defaults
#             adjusted (panel/plot/axes/legend) labels
# 2018-04-04: added default shift/stretch input file
# 2018-04-03: added smoothing span input panel parameters
#             added smoothing point count input panel parameters
# 2018-03-20: added slice width bar plot option
#             added gonad arm model input file & plot options
# 2018-02-28: increased smooth fit span as agreed with Filippos
# 2018-02-27: adjusted smooth fit parameters according to Filippos
#             added smooth fit span and number of points parameters
# 2017-10-23: replaced smooth fit color parameter by linetype parameter (color now used for
#             isoforms)
#             added isoform level input panel parameters / added isoform legend parameter for
#             profile plot / removed linetype parameters for profile plot (now used for isoforms)
# 2017-05-29: added sample stretch input panel parameters
# 2017-05-24: added minimum y-axis minimum parameter
# 2017-05-23: added minimum peak CPM input panel parameters
# 2017-05-22: added set.ylim option & manual ymin/ymax input panel parameters
# 2017-05-17: replaced heatmap options input panel parameters by abundance measure input panel
#             parameters
#             added row normalization input panel parameters / removed row.scaling heatmap option
# 2017-04-19: added distance metric input panel parameters
# 2017-04-18: added gene type input panel parameters
# 2017-04-12: (re-)replaced maximum gene rank by maximum number of genes parameter
#             replaced maximum number of genes by maximum gene rank parameter
# 2017-04-11: added gene list file import panel parameters
# 2017-04-10: added gene table XLSX export button parameters
# 2017-04-06: fixed capitalization of heatmap option input panel label
#             added gene table output panel parameters (gene.table.ngenes only)
# 2017-04-05: added log.transform plot option
#             added heatmap option input panel parameters (row.scaling only)
#             added gene cluster summary profile parameters
#             added gene cluster count input panel parameters
#             added heatmap dendrogram.side x/ylab.fontsize & color.values parameters
# 2017-03-29: added tab title, sample description & genotype input & heatmap output panel parameters
# 2017-03-19: replaced profile.plot.nrow parameter with plot columns count input panel parameters
# 2017-03-02: added single.y.scale plot option
#             added fix.xlim plot option & fixed x-axis limits parameter
#             added missing blank lines (cosmetics)
# 2017-02-24: added license comment
#             added sample shift input panel parameters
# 2017-02-23: re-ordered plot options
#             added raw.points plot option
#             added raw.lines plot option
#             added smooth.each plot option
#             added smooth.pooled plot option
#             added plot option input panel parameters (logscale only)
#             added sample names input panel parameters
#             added double-sourcing check / structured parameter definition / added profile plot
#             parameters
# 2017-02-21: added gene names input panel parameters
#             initial version (app title only)


#############
# libraries #
#############

# get viridis color palette for heatmap
library(viridis)


##############
# parameters #
##############

# ensure parameters are not defined already
if(!exists("params"))

  # define parameters
  params<-

    # store parameters in a named list
    list(


  ##################
  # user interface #
  ##################

    ##########################
    # titel panel parameters #
    ##########################

      # title of the app
      app.title = "SPACEGERM",
      app.subtitle.md =
        paste("**Spa**tial ***C****aenorhabditis* ***e****legans*",
              "**g**ermline **e**xpression of m**R**NA & **m**iRNA"),

    ###########################
    # gene profiles tab panel #
    ###########################

      # title of gene profiles
      gene.profiles.tab.title = "Expression profiles"


      #################
      # sidebar panel #
      #################

        ##########################
        # gene names input panel #
        ##########################

      # label of gene names input panel
      ,gene.names.input.label="Gene/miRNA names"

      # gene names input separator
      ,gene.names.input.separator=" "

      # default value of gene names input panel
      ,gene.names.input.default=paste("rpl-17"
                                     ,"iff-1"
                                     ,"perm-2"
                                     ,"perm-4"
                                     )

      # placeholder of gene names input panel
      ,gene.names.input.placeholder="enter gene names to plot (space separated)"


        #############################
        # isoform level input panel #
        #############################

      # label of isoform level input panel
      ,isoform.level.input.label="Expression level"

      # possible choices for isoform level input panel
      ,isoform.level.input.choices=c(`gene level estimates`=FALSE
                                    ,`isoform level estimates`=TRUE
                                    )

      # default choice for isoform level input panel
      ,isoform.level.input.default=FALSE,


        ##############################
        # abundance unit input panel #
        ##############################

       # label of abundance unit input panel
       abundance.unit.input.label = "Expression measure",

       # choices for abundance unit input panel
       abundance.unit.input.choices = c("CPM", "CPM / cell"),

       # default selection of abundance unit input panel
       abundance.unit.input.default = "CPM"


        ############################
        # sample names input panel #
        ############################

      # label of sample names input panel
      ,sample.names.input.label = "Samples",

      # default selection of sample names input panel
      sample.names.input.default =
        paste0("N2_mRNA_", rep(c("A", "P"), each = 3), 1:3),


        ############################
        # plot options input panel #
        ############################

      # label of plot option input panel
      plot.options.input.label = "Plot options"

      # plot options
      ,plot.options =

        c(`show raw data points`="raw.points",
          `show raw data lines`="raw.lines",
          `show per-sample smooth fits (LOESS)`="smooth.each",
          `show across-sample smooth fit (LOESS)`="smooth.pooled",
          `scale y-axis logarithmically (log2)`="logscale",
          `fix x-axis limits`="fix.xlim",
          `use single y-scale for all sub-plots`="single.y.scale",
          `manually set expression range limits` = "set.exprlim",
          `show slice width bars`="show.slice.width",
          `show gonad arm model (if single column & fixed x-axis limits)`=
            "show.model",
          `show smoothing standard error` = "show.smoothing.se",
          `show dropout slices (<10K pseudoaligned reads)` = "show.dropouts"),

      # default selection of plot option input panel
       plot.options.input.default =
         c("raw.points", "smooth.pooled", "fix.xlim", "show.smoothing.se", "show.dropouts"),


        #################################
        # expression range input panels #
        #################################

      manual.exprmin.input.label="Expression range minimum",
      manual.exprmin.input.min = -10^4,
      manual.exprmin.input.max = 10^6 - 1,
      manual.exprmin.input.default = 0,
      manual.exprmax.input.label = "Expression range maximum",
      manual.exprmax.input.max = 10^6,
      manual.exprmax.input.default = 10^4,


        ##################################
        # plot columns count input panel #
        ##################################

      # label of plot columns count input panel
      ncols.plot.input.label = "# plot columns"

      # minimum value of plot columns count input panel
      ,ncols.plot.input.min=1

      # default value of plot columns count input panel
      ,ncols.plot.input.default = 2,


        #####################################
        # smoothing point count input panel #
        #####################################

      smoothing.n.input.label = "# points to impute for smoothing",

      smoothing.n.input.min = 1,

      smoothing.n.input.default = 20,


        ##############################
        # smoothing span input panel #
        ##############################

      smoothing.span.input.label = "Span to use for smoothing",
      smoothing.span.input.min = 0,
      smoothing.span.input.max = 1,
      smoothing.span.input.default = .4,


        #############################
        # sample shifts input panel #
        #############################

      sample.shifts.input.min = -50,
      sample.shifts.input.max = 50,
      sample.shifts.input.default = 0,


        ################################
        # sample stretches input panel #
        ################################

      sample.stretches.input.min = .5,
      sample.stretches.input.max = 2,
      sample.stretches.input.default = 1,


      ##############
      # main panel #
      ##############

        #############################
        # profile plot output panel #
        #############################

      # There are currently no parameters for the plot output panel.
      # See the plot parameters below to adjust the plot itself.


    #####################
    # heatmap tab panel #
    #####################

      # title of heatmap tab panel
      heatmap.tab.title = "Spatial clustering"


      #################
      # sidebar panel #
      #################

        ##################################
        # sample description input panel #
        ##################################

      # label of sample names input panel
      ,sample.description.input.label = "Sample description",

      # default selection of sample names input panel
      sample.description.input.default = "gonad",


        ########################
        # genotype input panel #
        ########################

      # label of genotype input panel
      genotype.input.label = "Genotype",

      # default selection of genotype input panel
      genotype.input.default = "N2",


        #########################
        # gene type input panel #
        #########################

      # label of gene type input panel
      gene.type.input.label = "Gene type"

      # default selection of gene type input panel
      ,gene.type.input.default="mRNA"


        ##################################
        # gene cluster count input panel #
        ##################################

      # label of gene cluster count input panel
      ,nclust.genes.input.label="# gene clusters"

      # minimum value of gene cluster count input panel
      ,nclust.genes.input.min=1

      # default value of gene cluster count input panel
      ,nclust.genes.input.default=5


        #################################
        # abundance measure input panel #
        #################################

      # label of abundance measure input panel
      ,abundance.measure.input.label="Expression measure"

      # choices for abundance measure input panel
      ,abundance.measure.input.choices=c("CPM"
                                        ,"log2(1 + CPM)"
                                        ,"log10(1 + CPM)"
                                        )

      # default selection of abundance measure input panel
      ,abundance.measure.input.default="CPM"


        #################################
        # row normalization input panel #
        #################################

      # label of row normalization input panel
      ,row.normalization.input.label="Row normalization"

      # choices for row normalization input panel
      ,row.normalization.input.choices=c("scaling (z-scores)"
                                        ,"centering only"
                                        ,"none"
                                        )

      # default selection of row normalization input panel
      ,row.normalization.input.default="scaling (z-scores)"


        ###############################
        # distance metric input panel #
        ###############################

      # label of distance metric input panel
      ,distance.metric.input.label="Distance metric"

      # choices for distance metric input panel
      ,distance.metric.input.choices=c("1 - Pearson's r"
                                      ,"Euclidean distance"
                                      )

      # default selection of distance metric input panel
      ,distance.metric.input.default="1 - Pearson's r"


        ###############################
        # gene list file import panel #
        ###############################

      # label of gene list file import panel
      ,gene.list.import.label="Gene list"

      # MIME type to accept for gene list file import
      ,gene.list.file.import.mime.accept="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"

      # label to use for gene list file import button
      ,gene.list.file.import.button.label="import gene list"

      # placeholder to use for gene list file import
      ,gene.list.file.import.placeholder="use most variable genes"


        #################################
        # minimum peak CPM input panel  #
        #################################

      # label of minimum peak CPM input panel
      ,min.cpm.max.input.label="Peak CPM minimum"

      # minimum value of minimum peak CPM input panel
      ,min.cpm.max.input.min=0

      # maximum value of minimum peak CPM input panel
      ,min.cpm.max.input.max=10^6

      # default value of minimum peak CPM input panel
      ,min.cpm.max.input.default=0


      ##############
      # main panel #
      ##############

        ########################
        # heatmap output panel #
        ########################

      # There are currently no parameters for the heatmap output panel.
      # See the plot parameters below to adjust the plot itself.


        ###########################
        # gene table output panel #
        ###########################

      # (default) number of genes to include in gene table (per page)
      ,gene.table.ngenes=10


        #################################
        # gene table XLSX export button #
        #################################

      # label to use for gene table XLSX export button
      ,gene.table.xlsx.export.button.label="Export gene table as XLSX"

      # file name to use for gene table XLSX export
      ,gene.table.xlsx.name = "gene_table.xlsx",
      model3d.tab.title = "Virtual in-situ hybridization (vISH)",
      gene3d.input.label = "Gene/miRNA name",
      gene3d.input.default = "rpl-17",
      plot.options3d =
        c(`manually set expression range limis` = "set.exprlim"),
      plot.options3d.input.default = c(),


  ##############
  # data paths #
  ##############

      data.sqlite = "data.sqlite",

      # (relative) file path of Rds file with gene profiles
      gene.profiles.file = "gene.profiles.Rds",

      # (relative) file path of Rds file with gene profiles
      gonad.model.file = "gonad.model.Rds"


  ################
  # plot options #
  ################

    ################
    # profile plot #
    ################

      #########
      # input #
      #########

      # maximum number of samples to include in profile plot
      # Data will be colored by sample.
      # Thus, the maximum number of samples limits the valid choices of color palettes (see
      # profile.plot.brewer.palette parameter).
      ,profile.plot.max.nsamples=8


      ##############
      # processing #
      ##############

        #############
        # smoothing #
        #############

      # data smoothing method of profile plot
      ,profile.plot.smoothing.method="loess"


      ##############
      # appearance #
      ##############

        ########
        # text #
        ########

      # base font size of profile plot
      ,profile.plot.fontsize.base=18


        ##########
        # points #
        ##########

      # point size of profile plot
      ,profile.plot.pointsize = 2,
      profile.plot.dropouts.alpha = .5,


        #########
        # lines #
        #########

      # per-sample data line size of profile plot
      profile.plot.linesize.each = 1

      # smooth data line size of profile plot
      ,profile.plot.linesize.pooled=2

      # smooth data linetype of profile plot
      ,profile.plot.linetype.smooth="solid"


        ##########
        # colors #
        ##########

      # brewer palette of profile plot
      # See http://colorbrewer2.org for available palettes.
      # This palette will be used to color data by sample.
      # Thus, the palette must support at least as many data classes as the maximum number of
      # samples to be included (see profile.plot.max.nsamples parameter).
      # The palette should be qualitative, as there is no order in samples.
      # Also, it should be unpaired, as there are no meaningful sample pairs.
      # Dark2 and Set2 are the only non-paired qualitative palettes that are labelled 'colorblind
      # safe' for up to three data classes.
      # The Dark2 palette supports up to eigth data classes.
      ,profile.plot.brewer.palette="Dark2"


        ########
        # axes #
        ########

      # x-axis limits for fixed limit x-axis profile plot [%]
      ,profile.plot.xlim=c(0,100)


      ##########
      # labels #
      ##########

      # title of profile plot
      ,profile.plot.title="Spatial gene expression"

      # x-axis label of profile plot
      ,profile.plot.xlab="Position [% distal-to-proximal]"

      # y-axis label of profile plot
      ,profile.plot.ylab="Expression"

      # sample legend label of profile plot
      ,profile.plot.sample.legend.label="Sample"

      # isoform legend label of profile plot
      ,profile.plot.isoform.legend.label="Isoform"


    ###########
    # heatmap #
    ###########

      #########
      # input #
      #########

      # maximum number of genes to include in heatmap
      ,heatmap.nmax.genes=500


      ##############
      # appearance #
      ##############

        ##########
        # layout #
        ##########


      # side of heatmap to place dendrogram at
      ,heatmap.dendrogram.side="right"


        ########
        # text #
        ########

      # font size to use for heatmap x-axis labels
      ,heatmap.xlab.fontsize=8

      # font size to use for heatmap y-axis labels
      ,heatmap.ylab.fontsize=8


        ########
        # lines #
        ########

      # plotly scatter mode to use for gene cluster summary profiles
      # see https://plot.ly/javascript/reference/#scatter-mode for options
      ,heatmap.summary.mode="lines"

      # line width to use for gene cluster summary profiles
      ,heatmap.summary.linewidth=5


        ##########
        # colors #
        ##########

      # color values to use for heatmap tiles
      ,heatmap.color.values=

        # use viridis colors with heatmaply default settings
        viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "viridis"),


  plot.title.model3d = "Gonad arm model",
  colorscale.model3d = "Reds",
  colorlab.model3d = "Expression [CPM]",
  dplab.model3d = "distal-proximal [μm from DTC]",
  lrlab.model3d = "left-right [μm from center]",
  dvlab.model3d = "dorsal-ventral [μm from center]" ,
  eye.model3d.dp = 0,
  eye.model3d.lr = -3,
  eye.model3d.dv = 1.5)
