# tomo-seq shiny app parameter definition script
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

# file:         params.R
# author(s):    Marcel Schilling <marcel.schilling@mdc-berlin.de>
# created:      2017-02-21
# last update:  2017-04-18
# license:      GNU Affero General Public License Version 3 (GNU AGPL v3)
# purpose:      define parameters for tomo-seq shiny app


######################################
# change log (reverse chronological) #
######################################


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
require(viridis)


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
      app.title="tomo-seq"


    ###########################
    # gene profiles tab panel #
    ###########################

      # title of gene profiles
      ,gene.profiles.tab.title="Gene profiles"


      #################
      # sidebar panel #
      #################

        ##########################
        # gene names input panel #
        ##########################

      # label of gene names input panel
      ,gene.names.input.label="Gene names"

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


        ############################
        # sample names input panel #
        ############################

      # label of sample names input panel
      ,sample.names.input.label="Samples"

      # default selection of sample names input panel
      ,sample.names.input.default=c("ad_PG50um"
                                   ,"ad_PG20um"
                                   ,"ad_PG_wt_50"
                                   ,"ad_PG_wt_50_2"
                                   )


        ############################
        # plot options input panel #
        ############################

      # label of plot option input panel
      ,plot.options.input.label="Plot options"

      # plot options
      ,plot.options=c(`show raw data points`="raw.points"
                     ,`show raw data lines`="raw.lines"
                     ,`show per-sample smooth fits (LOESS)`="smooth.each"
                     ,`show across-sample smooth fit (LOESS)`="smooth.pooled"
                     ,`scale y-axis logarithmically (log2)`="logscale"
                     ,`fix x-axis limits`="fix.xlim"
                     ,`use single y-scale for all sub-plots`="single.y.scale"
                     )

      # default selection of plot option input panel
      ,plot.options.input.default=c("raw.points"
                                   ,"smooth.each"
                                   ,"smooth.pooled"
                                   ,"logscale"
                                   ,"single.y.scale"
                                   )


        ##################################
        # plot columns count input panel #
        ##################################

      # label of plot columns count input panel
      ,ncols.plot.input.label="# plot columns"

      # minimum value of plot columns count input panel
      ,ncols.plot.input.min=1

      # default value of plot columns count input panel
      ,ncols.plot.input.default=2


        #############################
        # sample shifts input panel #
        #############################

      # label of sample shift input panel
      ,sample.shifts.input.label="Sample shifts"

      # minimum value of sample shift input panel
      ,sample.shifts.input.min=-50

      # maximum value of sample shift input panel
      ,sample.shifts.input.max=50

      # default value of sample shift input panel
      ,sample.shifts.input.default=0

      # value suffix of sampleshift input panel
      ,sample.shifts.input.suffix="%"


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
      ,heatmap.tab.title="Heatmap"


      #################
      # sidebar panel #
      #################

        ##################################
        # sample description input panel #
        ##################################

      # label of sample names input panel
      ,sample.description.input.label="Sample description"

      # default selection of sample names input panel
      ,sample.description.input.default="anterior gonad"


        ########################
        # genotype input panel #
        ########################

      # label of genotype input panel
      ,genotype.input.label="Genotype"

      # default selection of genotype input panel
      ,genotype.input.default="wild type"


        #########################
        # gene type input panel #
        #########################

      # label of gene type input panel
      ,gene.type.input.label="Gene type"

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


        ###############################
        # heatmap options input panel #
        ###############################

      # label of heatmap option input panel
      ,heatmap.options.input.label="Heatmap options"

      # heatmap options
      ,heatmap.options=c(`log-transform gene expression (log2 1+CPM)`="log.transform"
                        ,`row-normalize expression matrix (z-scores)`="row.scaling"
                        )

      # default selection of heatmap option input panel
      ,heatmap.options.input.default=c("row.scaling"
                                      )


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
      ,gene.table.xlsx.name="gene_table.xlsx"


  ##############
  # data paths #
  ##############

      # (relative) file path of Rds file with tomo-seq data
      ,tomoseq.data.file="tomoseq.data.Rds"

      # (relative) file path of Rds file with gene profiles
      ,gene.profiles.file="gene.profiles.Rds"


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
      ,profile.plot.pointsize=2


        #########
        # lines #
        #########

      # per-sample data line size of profile plot
      ,profile.plot.linesize.each=1

      # raw data line type of profile plot
      ,profile.plot.linetype.raw="dashed"

      # smooth data line size of profile plot
      ,profile.plot.linesize.pooled=2

      # smooth data linetype of profile plot
      ,profile.plot.linetype.smooth="solid"


        ##########
        # colors #
        ##########

      # smooth data color of profile plot
      ,profile.plot.color.smooth="black"

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
      ,profile.plot.title="Spatial gene abundance"

      # x-axis label of profile plot
      ,profile.plot.xlab="distal-to-proximal position [%]"

      # y-axis label of profile plot
      ,profile.plot.ylab="gene abundance [CPM]"

      # sample legend label of profile plot
      ,profile.plot.sample.legend.label="sample"


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
        viridis(n=256,alpha=1,begin=0,end=1,option="viridis")

      # end parameter list definition
      )
