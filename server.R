# tomo-seq shiny app server script
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

# file:         server.R
# author(s):    Marcel Schilling <marcel.schilling@mdc-berlin.de>
# created:      2017-02-21
# last update:  2017-10-23
# license:      GNU Affero General Public License Version 3 (GNU AGPL v3)
# purpose:      define back end for tomo-seq shiny app


######################################
# change log (reverse chronological) #
######################################

# 2017-10-23: added user specified expression level (gene/isoform profiles?)
# 2017-10-17: replaced renderPlotly by (new) renderIheatmap
# 2017-05-29: added dynamic sample stretch input panel assignment & corresponding user input support
# 2017-05-23: added filtering of genes by peak CPM minimum specified by the user
# 2017-05-22: added support for user specified y-axis limits
# 2017-05-17: replaced user specified heatmap options by user specified abundance measure
#             added user specified row normalization choice
# 2017-04-19: added user specified distance metric choice
# 2017-04-18: added user specified gene type filtering
# 2017-04-13: added gene table annotation
# 2017-04-12: switched (back) from gene rank to count based filtering
#             moved gene rank based filtering out of heatmap function
#             moved genotype based filtering out of heatmap function
#             moved sample description based filtering out of heatmap function
#             moved gene list based filtering out of heatmap function
# 2017-04-11: added user specified gene list file input
# 2017-04-10: added gene table XLSX export button assignment
# 2017-04-06: added gene table output assignment
# 2017-04-05: added user specified heatmap option selection
#             added user specified gene cluster count
# 2017-03-29: added dynamic genotype input panel & heatmap output panel assignment
# 2017-03-19: added user specified plot columns count
# 2017-02-24: added license comment
#             added dynamic sample shift input panel assignment & corresponding user input support
# 2017-02-23: added user specified plot option selection
#             added user specified sample selection
#             replaced gene names output assignment by profile plot output assignment
# 2017-02-21: added gene names output assignment
#             initial version (empty template)

#############
# libraries #
#############

# get pipe operators
require(magrittr)

# get renderIheatmap
require(iheatmapr)


#############
# functions #
#############

# load functions
source("functions.R")


########
# data #
########

# load input data
source("data.R")


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

    # assign y-axis minimum input panel output
    output$manual.ymin.input<-

      # render y-axis minimum input panel
      renderUI(

        # take plot options selected by the user
        input$plot.options %>%

        # generate y-axis minimum input panel
        generate.manual.ymin.input

        # end y-axis minimum input panel rendering
        )

    # assign y-axis maximum input panel output
    output$manual.ymax.input<-

      # render y-axis maximum input panel
      renderUI(

        # take plot options selected by the user
        input$plot.options %>%

        # generate y-axis maximum input panel
        generate.manual.ymax.input(

          # pass on y-axis minimum specified by the user
          ymin=input$manual.ymin

          # end y-axis maximum input panel generation
          )

        # end y-axis maximum input panel rendering
        )

    # assign sample shifts input panel output
    output$shifts.input<-

      # render sample shifts input panel
      renderUI(

        # take sample names to include in plot
        input$sample.names %>%

        # generate sample shifts input panel
        generate.sample.shifts.input

        # end sample shifts input panel rendering
        )

    # assign sample stretches input panel output
    output$stretches.input<-

      # render sample stretches input panel
      renderUI(

        # take sample names to include in plot
        input$sample.names %>%

        # generate sample stretches input panel
        generate.sample.stretches.input

        # end sample stretches input panel rendering
        )

    # assign profile plot output
    output$profile.plot<-

      # render gene profiles plot
      renderPlot(

        # take tomo-seq data
        input.data$tomoseq.data %>%

        # generate profile plot
        generate.profile.plot(

          # plot profiles of genes specified by the user
          gene.names=input$gene.names

          # include sample specified by the user in plot
          ,sample.names=input$sample.names

          # set plot options specified by the user
          ,plot.options=input$plot.options

          # set manual y-axis limits specified by the user
          ,manual.ylim=c(input$manual.ymin,input$manual.ymax)

          # set plot columns count specified by the user
          ,ncols.plot=input$ncols.plot

          # set sample shifts specified by the user
          ,sample.shifts=

            # take sample names of samples included in plot
            input$sample.names %>%

            # extract corresponding sample shifts specified by user
            get.sample.shifts(input)

          # set sample stretches specified by the user
          ,sample.stretches=

            # take sample names of samples included in plot
            input$sample.names %>%

            # extract corresponding sample stretches specified by user
            get.sample.stretches(input)

          # set expression level specified by the user
          ,per.isoform=input$isoform.level

          # end profile plot generation
          )

        # end profile plot rendering
        )

    # assign genotype input panel output
    output$genotype.input<-

      # render genotype input panel
      renderUI(

        # take sample description to use for heatmap
        input$sample.description %>%

        # generate genotype input panel
        generate.genotype.input

        # end genotype input panel rendering
        )

    # assign gene type input panel output
    output$gene.type.input<-

      # render gene type input panel
      renderUI(

        # take sample description to use for heatmap
        input$sample.description %>%

        # generate gene type input panel for genotype to use for heatmap
        generate.gene.type.input(input$genotype)

        # end gene type input panel rendering
        )

    # assign gene list filtered gene profiles
    gene.profiles.filtered.gene.list<-

      # re-calculate gene list filtered gene profiles when necessary
      reactive(

        # take gene profiles
        input.data$gene.profiles %>%

        # extract gene profiles for genes in gene list file specified by the user
        filter.data.by.genes.file(input$gene.list.file)

        # end gene list filtered gene profiles re-calculation
        )

    # assign sample description filtered gene profiles
    gene.profiles.filtered.sample.description<-

      # re-calculate sample description filtered gene profiles when necessary
      reactive(

        # take gene list filtered gene profiles
        gene.profiles.filtered.gene.list() %>%

        # extract gene profiles for sample description specified by the user
        filter.data.by.sample.description(input$sample.description)

        # end sample description filtered gene profiles re-calculation
        )

    # assign genotype filtered gene profiles
    gene.profiles.filtered.genotype<-

      # re-calculate genotype filtered gene profiles when necessary
      reactive(

        # take sample description filtered gene profiles
        gene.profiles.filtered.sample.description() %>%

        # extract gene profiles for genotype specified by the user
        filter.data.by.genotype(input$genotype)

        # end genotype filtered gene profiles re-calculation
        )

    # assign gene type filtered gene profiles
    gene.profiles.filtered.gene.type<-

      # re-calculate gene type filtered gene profiles when necessary
      reactive(

        # take genotype filtered gene profiles
        gene.profiles.filtered.genotype() %>%

        # extract gene profiles for gene type specified by the user
        filter.data.by.gene.type(input$gene.type)

        # end gene type filtered gene profiles re-calculation
        )

    # assign minimum peak CPM filtered gene profiles
    gene.profiles.filtered.min.cpm.max<-

      # re-calculate minimum peak CPM filtered gene profiles when necessary
      reactive(

        # take gene type filtered gene profiles
        gene.profiles.filtered.gene.type() %>%

        # extract gene profiles for genes with peak CPM above minimum specified by the user
        filter.data.by.min.cpm.max(input$min.cpm.max)

        # end minimum peak CPM filtered gene profiles re-calculation
        )

    # assign top gene profiles
    gene.profiles.top<-

      # re-calculate top gene profiles when necessary
      reactive(

        # take minimum peak CPM filtered gene profiles
        gene.profiles.filtered.min.cpm.max() %>%

        # keep top varying genes
        keep.top.genes

        # end top gene profiles re-calculation
        )

    # assign heatmap object
    heatmap.object<-

      # re-calculate heatmap when necessary
      reactive(

        # take top gene profiles
        gene.profiles.top() %>%

        # generate heatmap
        generate.heatmap(

          # cluster genes into as many clusters as specified by the user
          nclust.genes=input$nclust.genes

          # set abundance measure specified by the user
          ,abundance.measure=input$abundance.measure

          # set row normalization specified by the user
          ,row.normalization=input$row.normalization

          # set distance metric specified by the user
          ,distance.metric=input$distance.metric

          # end heatmap generation
          )

        # end heatmap re-calculation
        )

    # assign heatmap output
    output$heatmap<-

      # render heatmap
      renderIheatmap(

        # take heatmap object
        heatmap.object()

        # end heatmap rendering
        )

    # assign gene annotation
    gene.annotation<-

      # re-calculate gene annotation when necessary
      reactive(

        # take gene profiles
        gene.profiles.top() %>%

        # extract gene annotation
        get.gene.annotation

        # end gene annotation re-calculation
        )

    # assign gene table object
    gene.table.object<-

      # re-calculate gene table when necessary
      reactive(

        # take heatmap object
        heatmap.object() %>%

        # generate gene table
        generate.gene.table(

          # annotate gene table with annotation extracted from gene profile table
          annotation=gene.annotation()

          # end gene table generation
          )

        # end gene table re-calculation
        )

    # assign gene table output
    output$gene.table<-

      # render gene table
      renderDataTable(

        # take gene table object
        gene.table.object()

        # set table rendering option
        ,options=

          # generate gene table options using helper function
          generate.gene.table.options

        # end gene table rendering
        )

    # assign gene table XLSX export button
    output$gene.table.xlsx.export.button<-

      # generate file download for gene table XLSX export
      downloadHandler(

        # set file name for file download
        filename=

          # use file name defined for gene table XLSX export
          get.gene.table.xlsx.name

        # set file content for file download
        ,content=

          # define gene table XLSX export file content generation function
          function(

            # out file name specified by downloadHandler
            file

            # end gene table XLSX export file content generation function parameter definition
            )

            # begin gene table XLSX export file content generation function definition
            {

              # use gene table object
              gene.table.object() %>%

              # save gene table to XLSX
              save.gene.table.xlsx(

                # set output file name for gene table XLSX
                output.xlsx=

                  # use file name specified by downloadHandler
                  file

                # end gene table XLSX export
                )

            # end gene table XLSX export file content generation function definition
            }

        # end file download generation for gene table XLSX export
        )

  # end shiny server function definition
  } %>%

# initialize shiny server
shinyServer
