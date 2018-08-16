# SPACEGERM shiny app server script
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

# file:         server.R
# author(s):    Marcel Schilling <marcel.schilling@mdc-berlin.de>
# created:      2017-02-21
# last update:  2018-08-16
# license:      GNU Affero General Public License Version 3 (GNU AGPL v3)
# purpose:      define back end for SPACEGERM shiny app


######################################
# change log (reverse chronological) #
######################################

# 2018-08-16: added support user specified location measure
# 2018-05-31: added support for gene profiles passed in as database query
#             added support for slice data passed in as database query
# 2018-05-17: replaced require by library
# 2018-05-16: renamed app for publication
# 2018-04-23: added 3D expression range inputs
# 2018-04-16: renamed y-axis limits inputs to expression range inputs
# 2018-04-13: added 3D model gene name default
#             replaced user specified sample shifts by defaults
#             replaced user specified sample description for heatmap by default
#             added user specified smoothing span for 3D model
#             added 3D model gene/genotype selection and CPM fitting support / fixed indentation
# 2018-04-09: removed sample stretches input
# 2018-04-03: added user specified smoothing span
#             added user specified smoothing point count
# 2018-03-21: added user specified abundance unit
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
library(magrittr)

# get renderIheatmap
library(iheatmapr)
library(plotly)


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
function(input, output, session){
  updateSelectizeInput(session, 'gene3d', choices = input.data$genes.name,
                       selected = params$gene3d.input.default, server = TRUE)

  output$manual.exprmin.input <-
    renderUI(generate.manual.exprmin.input(input$plot.options))
  output$manual.exprmax.input <-
    renderUI(generate.manual.exprmax.input(input$plot.options,
                                           exprmin = input$manual.exprmin))

  # cache model plot
  model.plot <- reactive(input.data$gonad.model %>%
                         plot.model(input$location.measure))

  # assign profile plot output
  output$profile.plot<-
    renderPlot({
      if(input$gene.names %in% "chicken!"){
        profile.plot <-
          ggdraw() +
            draw_image(paste0("http://www.factroom.ru/facts/wp-content/uploa",
                              "ds/2013/12/319-620x411.jpg")) +
            draw_label("There is no gene called chicken!", colour = "white",
                       fontface = "bold", y = .3, size = 32)
      } else {
        profile.plot <-
         input.data$slice.data %>%
         generate.profile.plot(
           gene.names = input$gene.names,
           sample.names = input$sample.names,
           plot.options = input$plot.options,
           manual.exprlim = c(input$manual.exprmin, input$manual.exprmax),
           ncols.plot = input$ncols.plot,
           per.isoform = input$isoform.level,
           unit = input$abundance.unit,
           location = input$location.measure,
           smoothing.n = input$smoothing.n,
           smoothing.span = input$smoothing.span,
           model2d = model.plot())}
      profile.plot})

  # assign genotype input panel output
  output$genotype.input <-
    renderUI(generate.genotype.input(params$sample.description.input.default))

  output$genotype3d.input <-
    renderUI(generate.genotype.input(params$sample.description.input.default,
                                     id = "genotype3d"))

  # assign gene type input panel output
  output$gene.type.input <-
    renderUI(
      generate.gene.type.input(params$sample.description.input.default,
                               input$genotype))
  output$manual.exprmin3d.input <-
    renderUI(generate.manual.exprmin.input(input$plot.options3d,
                                           id = "manual.exprmin3d"))
  output$manual.exprmax3d.input <-
    renderUI(generate.manual.exprmax.input(input$plot.options3d,
                                           exprmin = input$manual.exprmin3d,
                                           id = "manual.exprmax3d"))

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
      filter.data.by.sample.description(params$sample.description.input.default))

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
      collect %>%

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

  output$model3d <-
    renderPlotly(
      plot.model3d(
        outline = input.data$gonad.model$outline,
        cpm.fit = input.data$slice.data %>%
                  filter(gene.name == input$gene3d,
                         genotype == input$genotype3d) %>%
                  collect %>%
                  fit.cpm(model.length = max(input.data$gonad.model$outline$dp),
                          smoothing.span = input$span3d),
        plot.options = input$plot.options3d,
        manual.exprlim = c(input$manual.exprmin3d, input$manual.exprmax3d)))

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
