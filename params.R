#######################
# general information #
#######################

# file:         params.R
# author(s):    Marcel Schilling <marcel.schilling@mdc-berlin.de>
# created:      2017-02-21
# last update:  2017-02-23
# purpose:      define parameters for tomo-seq shiny app


######################################
# change log (reverse chronological) #
######################################

# 2017-02-23: added sample names input panel parameters
#             added double-sourcing check / structured parameter definition / added profile plot
#             parameters
# 2017-02-21: added gene names input panel parameters
#             initial version (app title only)


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

    ##############
    # main panel #
    ##############

      #############################
      # profile plot output panel #
      #############################

      # There are currently no parameters for the plot output panel.
      # See the plot parameters below to adjust the plot itself.


  ##############
  # data paths #
  ##############

      # (relative) file path of Rds file with tomo-seq data
      ,tomoseq.data.file="tomoseq.data.Rds"


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

      # raw data line size of profile plot
      ,profile.plot.linesize.raw=1

      # raw data line type of profile plot
      ,profile.plot.linetype.raw="dashed"

      # smooth data line size of profile plot
      ,profile.plot.linesize.smooth=2

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


      ##########
      # layout #
      ##########

      # row count of profile plot
      ,profile.plot.nrow=2


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

      # end parameter list definition
      )
