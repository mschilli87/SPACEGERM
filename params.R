#######################
# general information #
#######################

# file:         params.R
# author(s):    Marcel Schilling <marcel.schilling@mdc-berlin.de>
# created:      2017-02-21
# last update:  2017-02-21
# purpose:      define parameters for tomo-seq shiny app


######################################
# change log (reverse chronological) #
######################################

# 2017-02-21: added gene names input panel parameters
#             initial version (app title only)


##############
# parameters #
##############

# define parameters
params<-

  # store parameters in a named list
  list(

    # title of the app
    app.title="tomo-seq"

    # label of gene names input panel
    ,gene.names.input.label="Gene names"

    # default value of gene names input panel
    ,gene.names.input.default=paste("rpl-17"
                                   ,"iff-1"
                                   ,"perm-2"
                                   ,"perm-4"
                                   )

    # placeholder of gene names input panel
    ,gene.names.input.placeholder="enter gene names to plot (space separated)"

    # end parameter list definition
    )
