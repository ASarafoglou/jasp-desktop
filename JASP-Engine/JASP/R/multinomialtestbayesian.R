#
# Copyright (C) 2017 University of Amsterdam
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

MultinomialTestBayesian <- function(dataset = NULL, options, perform = "run",
               callback = function(...) 0,  ...) {

  # First, we load the variables into the R environment
  factor <- NULL
  asnum <- NULL
  if (options$factor != "") {
    factor <- options$factor
    if (options$counts != "") {
      asnum <- options$counts
      if (options$exProbVar != "") {
        asnum <- c(asnum, options$exProbVar)
      }
    }
  }



  if (is.null(dataset)) {
    if (perform == "run") {
      dataset <- .readDataSetToEnd(columns.as.numeric=asnum,
                                   columns.as.factor=factor,
                                   exclude.na.listwise=NULL)
    } else {
      dataset <- .readDataSetHeader(columns.as.numeric=asnum,
                                    columns.as.factor=factor)
    }
  } else {
    dataset <- .vdf(dataset, columns.as.numeric=asnum,
                    columns.as.factor=factor)
  }

  results <- list() # Initialise results object


  # Then, we retrieve the state and initialise the output objects
  state <- .retrieveState()

  bayesMultResults <- NULL # result of the Bayesian multinomial test
  bayesMultTable <- NULL # chi-square table
  bayesDescriptivesTable <- NULL # expected versus observed
  bayesDescriptivesPlot <- NULL # barplot of factor levels
  bayesResultsPlot <- NULL # posterior estimates and CI


  # Then, we can fill the output objects with old info if its option did not
  # change.
  if (!is.null(state)) {
    diff <- .diff(options, state$options) # a list of TRUE/FALSE

    if (is.list(diff)){

      if (!any(diff[["factor"]], diff[["counts"]],
               diff[["credibleIntervalInterval"]],
               diff[["hypothesis"]], diff[["exProbVar"]],
               diff[["expectedProbs"]])){

        bayesMultResults <- state[["bayesMultResults"]]

        # the following depend on bayesMultResults so in same if statement
        if (!any(diff[["credibleInterval"]],
                 diff[["countProp"]])) {
          bayesDescriptivesTable <- state[["bayesDescriptivesTable"]]
        }

        if (!any(diff[["bayesDescriptivesPlotcredibleInterval"]],
                 diff[["countProp"]], diff[["plotWidth"]],
                 diff[["plotHeight"]])) {
          bayesDescriptivesPlot <- state[["bayesDescriptivesPlot"]]
        }

      }

      #... etcetera
      # TODO
    }

  }

  # Meta information
  results[["title"]] <- "Bayesian Multinomial Test"
  results[[".meta"]] <- list(list(name = "bayesMult", type = "table"),
                             list(name = "bayesDescriptivesTable", type = "table"),
                             list(name = "bayesDescriptivesPlot", type = "image"))

  # chi-square Table
  # Generate results
  if (is.null(bayesMultResults)) {
    bayesMultResults <- .bayesMultinomialTest(dataset, options, factor, perform)
  }

  results[["bayesMult"]] <- .bayesMultTable(bayesMultResults, options, perform)


  # Descriptives Table
  if (options[["descriptives"]]) {
    # Generate descriptives table
    if (is.null(bayesDescriptivesTable)) {
      bayesDescriptivesTable <- .multinomialDescriptives(bayesMultResults, factor, options, perform)
    }

    results[["bayesDescriptivesTable"]] <- bayesDescriptivesTable

  } else {

    results[["bayesDescriptivesTable"]] <- NULL

  }

  # Multinomial Descriptives Plot
  if (options[["bayesDescriptivesPlot"]]) {
    # Generate descriptives plots
    if (is.null(bayesDescriptivesPlot)) {
      bayesDescriptivesPlot <- .multinomialbayesDescriptivesPlot(bayesMultResults, options, perform)
    }

    plotPath <- list(bayesDescriptivesPlot$data) # for keep later

    results[["bayesDescriptivesPlot"]] <- bayesDescriptivesPlot

  } else {

    results[["bayesDescriptivesPlot"]] <- NULL
    plotPath <- list()

  }


  if (perform == "run") {

    state <- list()
    state[["options"]] <- options
    state[["bayesMultResults"]] <- bayesMultResults
    state[["bayesDescriptivesTable"]] <- bayesDescriptivesTable
    state[["bayesDescriptivesPlot"]] <- bayesDescriptivesPlot

    return(list(results=results, status="complete", state=state,
                keep = plotPath))

  } else {

    return(list(results=results, status="inited", state=state,
                keep = plotPath))

  }



}

# Run Bayesian multinomial test and return object
.bayesMultiomialTest <- function(dataset, options, factor, perform){

  bayesMultResults <- NULL

  if (perform == "run" && !is.null(factor)) {
    # first determine the hypotheses
    f <- dataset[[.v(factor)]]
    f <- f[!is.na(f)]
    nlev <- nlevels(f)

    if (options$counts != ""){
      # convert to "regular" factor
      c <- dataset[[.v(options$counts)]]
      f <- rep(f, c)
    }

    # Create table in order of appearance of the dataset
    t <- table(f)
    val <- t[unique(f)]

    # prior concentration (based on user input)
    if(options$priorConcentration == ''){
      a <- rep(1, nlev)
    } else {
      a <- options$priorConcentration
    }
    
    # create a named list with Bayesian multinomial result object

    # Multinomial Test
    if(options$hypothesis == 'bayesMultinomialTest' | options$hypothesis == 'bayesMultExpectedProbTest'){
      
        # catch warning message and append to object if necessary
        bmr  <- NULL
        warn <- NULL
        
        # Extract hypothesis
        if(options$hypothesis == 'bayesMultinomialTest'){
          hyp     <- rep(1/nlev, nlev)
        } else{
          hyp <- .generateExpectedProbs(dataset, options, nlevels)
        }
        
        # check if probabilities sum to 1
        sumsToOne <- sum(hyps) == 1
        if(sumsToOne == FALSE){
          warn <- 'Probabilites were rescaled to sum to 1.'
        }
        
        # calculate Savage-Dickey BF
        beta.xa <- prod(gamma(val + a))/gamma(sum(val + a)) # Beta(x + a)
        beta.a  <- prod(gamma(a))/gamma(sum(a))             # Beta(a)
        prior     <- (1/beta.a)  * (prod(hyp^(a - 1)))
        posterior <- (1/beta.xa) * (prod(hyp^((val + a) - 1)))
        BF        <- prior/posterior
        bmr       <- list(nlevels = length(h),
                          BF      = BF)
        bmr[["warn"]] <- warn
        bayesMultResults <- bmr
        
    } else if(options$hypothesis == 'bayesMultOrderConstrTest'){
      
      factorlevels <- levels(f)
      # step 1: produce adjacency matrix and perform check (is check needed?)
      A <- .toAmat(hyp, factorlevels)
      A.check <- .isAsyclic(A)
      # step 2: truncated sampling
      post.samples <- .truncatedSampling(options$a, A, niter = options$niter)
      # step 3: prepare samples for bridge sampling
      # step 4: bridge sampling and calculation of Bayes factor
      
    }
  }

  # return the output object
  return(bayesMultResults)
}

# Transform Bayesian multinomial test object into table for JASP
# bayesMultResults = list(H1 = obj, H2 = obj, ....)
.bayesMultTable <- function(bayesMultResults, options, perform){
  table <- list()
  footnotes <- .newFootnotes()
  table[["title"]] <- "Multinomial Test"

  # include fields
  fields <- list(
    list(name="case", title="", type="string", combine=TRUE),
    list(name="chisquare", title="\u03C7\u00B2", type = "integer", format = "sf:4;dp:3"),
    list(name="df", title="df", type="integer"),
    list(name="p", title="p", type="number", format="dp:3;p:.001")
    )

  # include footnotes

  table[["schema"]] <- list(fields = fields)

  message <- list()

  for(r in 1:length(bayesMultResults)){

    if(!is.null(bayesMultResults[[r]][["warn"]])) {

      message[[r]] <- bayesMultResults[[r]][["warn"]]
      .addFootnote(footnotes, symbol="<em>Note.</em>", text=message)
    }
  }

  table[["footnotes"]] <- as.list(footnotes)

  # fill in results one row at a time
  if (!is.null(bayesMultResults)){

  for(r in 1:length(bayesMultResults)){
      table[["data"]][[r]] <- list(case = names(bayesMultResults)[r],
                                   chisquare = bayesMultResults[[r]][["statistic"]][["X-squared"]],
                                   df = bayesMultResults[[r]][["parameter"]][["df"]],
                                   p = bayesMultResults[[r]][["p.value"]])

      if (options$VovkSellkeMPR){
        for (row in 1:length(table[["data"]])){
          table[["data"]][[row]][["VovkSellkeMPR"]] <- .VovkSellkeMPR(table[["data"]][[row]][["p"]])
        }
      }
      table[["status"]] <- "complete"
    }

  } else {
    # init state?
    data <- list()

    if(is.null(bayesMultResults[[r]])){
      htables <- ""
    }
    # for (h in htables){
    #   if (options$VovkSellkeMPR){
    #     data[[length(data) + 1]] <- list(case=h, chisquare=".", df=".", p=".",
    #                                      VovkSellkeMPR=".", lowerCI=".",
    #                                      upperCI=".")
    #   } else {
    #     data[[length(data) + 1]] <- list(case=h, chisquare=".", df=".", p=".",
    #                                      lowerCI=".", upperCI=".")
    #   }
    #  }

    table[["data"]] <- data
  }

  return(table)
}

# Create multinomial descriptives table
.multinomialDescriptivesTable <- function(bayesMultResults, factor, options, perform) {
  if (options[["countProp"]]=="descCounts"){
    numberType = list(type="integer")
  } else {
    numberType = list(type="number", format="sf:4;dp:3")
  }

  # Expected vs. Observed table
  table <- list("title" = "Descriptives table")


  if (is.null(factor)){
    # If we have no variable init table with generic name

    fields <- list(
      list(name="factor", title="Factor", type = "string"),
      c(list(name="observed", title="Observed"), numberType),
      c(list(name="expected", title="Expected"), numberType)
    )
    if (options$credibleInterval){
      interval <- 100 * options$credibleIntervalInterval
      title <- paste0(interval, "% Confidence Interval")
      fields[[length(fields)+1]] <- list(name="lowerCI",
                                         title="Lower",
                                         type = "number",
                                         format = "sf:4;dp:3",
                                         overTitle = title)
      fields[[length(fields)+1]] <- list(name="upperCI",
                                        title="Upper",
                                        type = "number",
                                        format = "sf:4;dp:3",
                                        overTitle = title)
    }
    rows <- list(list(factor = ".", observed = ".", expected = "."))


  } else if (perform != "run") {
    # If we have a variable then init table with factor name

    fields <- list(
      list(name="factor", title=factor, type = "string"),
      c(list(name="observed", title="Observed"), numberType),
      c(list(name="expected", title="Expected"), numberType)
    )
    if (options$credibleInterval){
      interval <- 100 * options$credibleIntervalInterval
      title <- paste0(interval, "% Confidence Interval")
      fields[[length(fields)+1]] <- list(name="lowerCI",
                                         title="Lower",
                                         type = "number",
                                         format = "sf:4;dp:3",
                                         overTitle = title)
      fields[[length(fields)+1]] <- list(name="upperCI",
                                        title="Upper",
                                        type = "number",
                                        format = "sf:4;dp:3",
                                        overTitle = title)
    }
    rows <- list(list(factor = ".", observed = ".", expected = "."))

  } else {

    # now  we want to create the full table

    # First we create the correct columns
    fields <- list(
      list(name="factor", title=factor, type = "string"),
      c(list(name="observed", title="Observed"), numberType)
    )

    nms <- names(bayesMultResults)

    if (length(nms) == 1) {
      fields[[length(fields)+1]] <- c(list(name="expected",
                                           title = paste0("Expected: ", nms)),
                                      numberType)
    } else {
      for (i in 1:length(nms)) {
        fields[[length(fields)+1]] <- c(list(name=nms[i],
                                             title = nms[i]),
                                        numberType,
                                        overTitle = "Expected")
      }
    }

    if (options$credibleInterval){
      interval <- 100 * options$credibleIntervalInterval
      title <- paste0(interval, "% Confidence Interval")
      fields[[length(fields)+1]] <- list(name="lowerCI",
                                         title="Lower",
                                         type = "number",
                                         format = "sf:4;dp:3",
                                         overTitle = title)
      fields[[length(fields)+1]] <- list(name="upperCI",
                                        title="Upper",
                                        type = "number",
                                        format = "sf:4;dp:3",
                                        overTitle = title)
    }

    # Then we fill the columns with the information
    if (options[["countProp"]]=="descCounts"){
      div <- 1
    } else {
      div <- sum(bayesMultResults[[1]][["observed"]])
    }

    tableFrame <- data.frame(factor = names(bayesMultResults[[1]][["observed"]]),
                             observed = as.integer(bayesMultResults[[1]][["observed"]])/div,
                             stringsAsFactors = FALSE)


    for (r in bayesMultResults){
      tableFrame <- cbind(tableFrame, as.integer(r[["expected"]])/div)
    }

    if (length(nms) == 1) {
      colnames(tableFrame)[-(1:2)] <- "expected"
    } else {
      colnames(tableFrame)[-(1:2)] <- nms
    }

    # Add credibleInterval
    if (options$credibleInterval){
      n <- sum(bayesMultResults[[1]][["observed"]])
      # make a list of cis
      ci <- lapply(bayesMultResults[[1]][["observed"]], function(l) {
        bt <- binom.test(l,n,conf.level = options$credibleIntervalInterval)
        return(bt$conf.int * n) # on the count scale
      })

      # add these to the tableFrame
      ciDf <- t(data.frame(ci))
      colnames(ciDf) <- c("lowerCI", "upperCI")
      tableFrame <- cbind(tableFrame, ciDf/div)
    }

    rows <- list()

    for (i in 1:nrow(tableFrame)){
      rows[[i]] <- as.list(tableFrame[i,])
    }
    table[["status"]] <- "complete"

  }

  table[["schema"]] <- list(fields = fields)
  table[["data"]] <- rows

  return(table)
}

# Create multinomial descriptives plot, histogram
.multinomialbayesDescriptivesPlot <- function(bayesMultResults, options, perform) {

  # init output object
  bayesDescriptivesPlot <- list("title" = "Descriptives plot")
  bayesDescriptivesPlot[["width"]] <- options$plotWidth
  bayesDescriptivesPlot[["height"]] <- options$plotHeight
  bayesDescriptivesPlot[["custom"]] <- list(width = "plotWidth",
                                      height = "plotHeight")

  if (perform == "run"){
    # Generate the plot

    # Counts or props
    if (options[["countProp"]]=="descCounts"){
      div <- 1
      yname <- "Observed counts"
    } else {
      div <- sum(bayesMultResults[[1]][["observed"]])
      yname <- "Observed proportions"
    }

    # Observed values
    obs <- bayesMultResults[[1]][["observed"]]/div

    # Calculate confidence interval
    cl <- options$bayesDescriptivesPlotcredibleInterval
    n <- sum(bayesMultResults[[1]][["observed"]])
    ci <- lapply(bayesMultResults[[1]][["observed"]], function(l) {
      bt <- binom.test(l,n,conf.level = cl)
      return(bt$conf.int * n) # on the count scale
    })
    ciDf <- data.frame(t(data.frame(ci)))/div
    colnames(ciDf) <- c("lowerCI", "upperCI")


    # Define custom y axis function
    base_breaks_y <- function(x){
      b <- pretty(c(0,x))
      d <- data.frame(x=-Inf, xend=-Inf, y=min(b), yend=max(b))
      list(ggplot2::geom_segment(data=d, ggplot2::aes(x=x, y=y,
                                                      xend=xend, yend=yend),
                                 size = 0.75,
                                 inherit.aes=FALSE),
           ggplot2::scale_y_continuous(breaks=b))
    }


    # Create plot
    p <- ggplot2::ggplot(mapping = ggplot2::aes(x = factor(names(obs),
                                                           levels = names(obs)),
                                                y = as.numeric(obs))) +
      ggplot2::geom_bar(stat = "identity", size = 0.75, colour="black", fill = "grey") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=ciDf$lowerCI, ymax = ciDf$upperCI),
                             size = 0.75, width = 0.3) +
      base_breaks_y(ciDf$upperCI) +
      ggplot2::xlab(options$factor) +
      ggplot2::ylab(yname) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size = 18),
        panel.grid.major = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_text(size = 18, vjust=0.1),
        axis.title.y = ggplot2::element_text(size = 18, vjust=0.9),
        axis.text.x = ggplot2::element_text(size = 15),
        axis.text.y = ggplot2::element_text(size = 15),
        panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
        plot.background = ggplot2::element_rect(fill = "transparent", colour = NA),
        legend.background = ggplot2::element_rect(fill = "transparent", colour = NA),
        panel.border = ggplot2::element_blank(),
        axis.line =  ggplot2::element_blank(),
        legend.key = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_line(size = 0.5),
        axis.ticks.length = grid::unit(3, "mm"),
        axis.ticks.margin = grid::unit(1,"mm"),
        plot.margin = grid::unit(c(0.1, 0.1, 0.6, 0.6), "cm"),
        legend.position = "none")


    # create plot object
    content <- .writeImage(width = options$plotWidth,
                           height = options$plotHeight,
                           plot = p, obj = TRUE)

    bayesDescriptivesPlot[["convertible"]] <- TRUE
    bayesDescriptivesPlot[["obj"]] <- content[["obj"]]
    bayesDescriptivesPlot[["data"]] <- content[["png"]]
    bayesDescriptivesPlot[["status"]] <- "complete"

  } else {
    bayesDescriptivesPlot[["data"]] <- ""
  }

  return(bayesDescriptivesPlot)

}

# Parse expected probabilities/counts
.generateExpectedProbs <- function(dataset, options, nlevels){
  # This function returns a vector of expected probabilities.

  if (options$exProbVar != ""){
    # use only exProbVar
    eDf <- dataset[[.v(options$exProbVar)]]

    if (nlevels != nrow(eDf)){
      stop("Expected counts do not match number of levels of factor!")
    }

    return(eDf)

  } else {

    stop("No expected counts entered!")

  }

}

# Helper functions for the multinomial order constrained hypothesis test
# Produce adjacency matrix
.toAmat <- function(OR, factorlevels){
  # input: character vector with order constrained hypothesis
  # inputOrderRestrictions <- c('A', '<', 'B', '=', 'C', '=', 'H', '=', 
  #                             'I', '<', 'E',',', 'F', '<', 'D', '<', 'X', 
  #                             '=', 'Z')
  # output: adjacency matrix
  
  otherConstraintsIndex <- which(stringr::str_detect(OR, '[,=]+'))
  ORsplit <- OR[-otherConstraintsIndex]
  orderConstraintsIndex <- which(stringr::str_detect(ORsplit, '[<]')) 
  
  # create object which encodes order constraints
  ORvector <- NULL
  lower <- c(1, orderConstraintsIndex + 1)
  upper <- c(orderConstraintsIndex - 1, length(ORsplit))
  for(i in 1:length(lower)){
    ORvector[i] <- stringr::str_c(ORsplit[lower[i]:upper[i]], collapse = ',')
  }
  
  # create object of vertices names
  verticesNames <- stringr::str_c(OR, collapse = '')
  verticesNames <- unlist(strsplit(verticesNames, '[,|<]+'))
  verticesNames <- stringr::str_replace_all(verticesNames, "=", ",")    
  # include all factor levels in verticesNames
  for(i in seq_along(factorlevels)){
    isElement <- any(grepl(factorlevels[i], verticesNames) == TRUE)
    if(isElement == FALSE) verticesNames <- c(verticesNames, factorlevels[i])
  }
  
  # construct and fill in adjacency matrix
  nVertices <- length(verticesNames)
  adjM <- matrix(0, ncol = nVertices, nrow = nVertices)
  colnames(adjM) <- rownames(adjM) <- verticesNames
  
  # for every factor level
  for(i in seq_along(factorlevels)){
    # if factor level is order restricted
    source <- factorlevels[i]
    isOrderRestricted <- any(grepl(factorlevels[i], ORvector) == TRUE)
    if(isOrderRestricted){
      # find its sink
      sourceIndex <- grep(factorlevels[i], verticesNames)
      sink <- unlist(strsplit(ORvector[grep(factorlevels[i], ORvector) + 1], ','))
      if(!any(is.na(sink))){
        sinkIndex <- sapply(sink, function(x) grep(x, verticesNames))
        # and fill out adjacency matrix
        adjM[cbind(sourceIndex, sinkIndex)] <- 1
      } 
    }
  }
  return(adjM)
}

# Check adjacency matrix for transitivity
.isAsyclic <- function(amat) {
  
  # input: adjacency matrix
  # output: length 0 vector if it's acyclic (transitivity holds)
  #         Vector of vertex names if intransitive
  # algorithm based on
  # https://math.stackexchange.com/questions/513288/test-for-acyclic-graph-property-based-on-adjacency-matrix
  
  v <- ncol(amat)
  error <- FALSE
  
  for (i in 1:v) {
    
    d <- diag(amat %^% i)
    if (any(d != 0)) {
      error <- colnames(amat)[which(d!=0)]
    } 
    error <- stringr::str_c(error, collapse = ',')
  }
  
  
  if(error != FALSE){
    output <- paste('Error. Order restriction not possible in:', error, sep = ' ')
  } else {
    output <- TRUE
  }
  
  return(output)
}

# Do truncated sampling
.truncatedSampling <- function(a, A = matrix(NA, ncol=length(a)), niter = 1e6) {
  # input: a = prior or posterior Dirichlet paramters (can be either)
  #        A = Adjacency matrix containing order restrictions
  # output: posterior samples from truncated dirichlet distribution
  
  nburnin <- 10 
  
  allParameterNames <- strsplit(colnames(A), ',')
  nEquals           <- sapply(allParameterNames, length)
  equalities        <- which(stringr::str_count(colnames(A)) != 1)
  # postSamples       <- matrix(ncol=length(unlist(allParameterNames)), nrow = niter + nburnin)
  # colnames(postSamples) <- unlist(allParameterNames)
  postSamples <- matrix(ncol=ncol(A), nrow = (niter + nburnin))
  colnames(postSamples) <- colnames(A)
  # starting values of Gibbs Sampler
  k  <- ncol(A)
  z  <- rgamma(k, a, 1)
  iteration <- 0
  
  for(iter in 1:(niter+nburnin)){
    
    for(i in 1:k){
      
      # 1. if parameter is unrestricted
      #    (no entry in ith row or column of A),
      #    sample from unrestricted gammadistribution
      if(sum(A[,i])==0 & sum(A[i,])==0){
        z[i] <- rgamma(1, a[i], 1)
        
        # 2. if parameter is has order constraints
        #    sample from truncated gamma distribution
      } else {
        
        v <- runif(1, 0, exp(-z[i]))
        
        # Lo: lower bound
        Lo <- 0
        # check for lower bound
        if(1 %in% A[,i]){
          # check for accompaning lower bound (in column)
          smallerValue <- which(A[,i] == 1)
          Lo           <- max(z[smallerValue])
        }
        # Hi: upper bound
        Hi <- -log(v)
        # check for a upper bound (in row)
        if(1 %in% A[i,]) {
          # check for accompaning upper bound (denoted as 1 in ith column of A)
          biggerValue <- which(A[i,] == 1)
          Hi <- min(z[biggerValue], Hi)
        }
        
        z[i] <- (runif(1)*(Hi^a[i] - Lo^a[i]) + Lo^a[i])^{1/a[i]}
      }
    }
    
    #  zAll <- rep(z, nEquals)
    
    # 4. transform Gammato Dirichlet samples
    # postSamples[iter,] <- zAll/sum(zAll)
    postSamples[iter,] <- z/sum(z)
    
    # show progress
    if (iter %% 100000 == 0) {
      iteration <- iteration + 1
      print(iteration)
    }
  }
  postSamples <- t(apply(postSamples, 1, function(x) x/nEquals))
  postSamples <- postSamples[-(1:nburnin), ]
  return(postSamples)
}

# To Do: Move parameters to real line 
# To Do: Perform bridge sampling 

.multinomialPosteriorPlot <- function(factor, options, run, medianPlusHDI, nlevels){
  # plots the posterior parameter estimates against either the prior distribution
  # or against the posterior distribution of the encompassing model

  # Initialisation plot
  .initSplitPlot <- function(){
    plot(1, type='n', xlim=0:1, ylim=0:1, bty='n', axes=FALSE, xlab="",
         ylab="")
    axis(2, at=0:1, labels=FALSE, cex.axis= 1.4, ylab="")
    mtext(text = 'Parameter estimates', side = 1, cex=1.5, line = 3)
  }

  # Plot
  posteriorPlot <- list()
  
  if (run == FALSE){
    
    image <- .writeImage(width = options$plotWidth, height = options$plotHeight,
                         plot = .initPosteriorPlot, obj = FALSE)
    posteriorPlot[["data"]] <- image[["png"]]
    
  } else {
  
  model1 <- medianPlusHDI[['Model1']]
  posterior <- medianPlusHDI[['Posterior']]
  level.names <- levels(factor)

  x.breaks <- pretty(c(min(model1$lower, posterior$lower), max(model1$upper, posterior$upper)), 5)
  d.y <- data.frame(y=-Inf, yend=-Inf, x=1, xend=nlevels)
  d.x <- data.frame(x=-Inf, xend=-Inf, y=x.breaks[1], yend=x.breaks[length(x.breaks)])
  
  p <- ggplot2::ggplot(posterior, aes(x = 1:nlevels, 
                             y = posterior$median, 
                          ymin = posterior$lower, 
                          ymax = posterior$upper
                          )) + 
    ggplot2::scale_y_continuous('Parameter estimates' , breaks = x.breaks) +
    ggplot2::scale_x_continuous(name = '', breaks = 1:nlevels, labels = level.names, limits = c(1,nlevels)) +
    ggplot2::coord_flip() + 
    ggplot2::geom_segment(data=d.x, ggplot2::aes(x=x, y=y, xend=xend, yend=yend),
                          size = 0.75,
                          inherit.aes=FALSE) +
    ggplot2::geom_segment(data=d.y, ggplot2::aes(x=x, y=y, xend=xend, yend=yend),
                          size = 0.75,
                          inherit.aes=FALSE) +
    ggplot2::geom_ribbon(aes(ymin=model1$lower, ymax=model1$upper), alpha = .9, fill = 'lightgrey') +
    ggplot2::geom_ribbon(aes(ymin=posterior$lower, ymax=posterior$upper), alpha = .6, fill = 'grey36') +
    ggplot2::geom_pointrange() + 
    ggplot2::geom_point(shape = 21, color = "black", fill = "grey", size = 4.5) +
    ggplot2::theme_classic(base_size = 20) + 
    ggplot2::theme(axis.line = element_blank())
  
  image <- .writeImage(width = options$plotWidth,
                       height = options$plotHeight,
                       plot = p)
  
  posteriorPlot[["data"]] <- image[["png"]]
  posteriorPlot[["obj"]] <- image[["obj"]]
  posteriorPlot[["convertible"]] <- TRUE
  posteriorPlot[["status"]] <- "complete"
  
  }
 
    return(posteriorPlot)
}
# calculation of median and highest density interval
.multinomialCredibleIntervalPlusMedian <- function(ciInterval = .95, a, counts, nlevels, 
                                                   model1.samples = NULL, post.samples = NULL) {
  # input : user input CI interval (default = .95)
  # a     : user input prior concentration as specified by user
  # counts: user input counts per category
  # post.samples: if user specified order constrained hypothesis
  # model1.samples: if user specified order constrained hypothesis
  
  output <- list(Model1 = data.frame(lower = NA, median = NA, upper = NA), 
                 Posterior = data.frame(lower = NA, median = NA, upper = NA))
  lower <- (1 - ciInterval) / 2
  upper <- 1 - lower
  N     <- sum(counts)
  A     <- sum(a)
  
  # Model 1 (e.g., prior samples or posterior samples from encompassing model)
  # equality constrained hypothesis: median and CI are calculated analytically
  if(is.null(model1.samples)){
    for(i in 1:nlevels){
      d <- qbeta(c(lower, .5, upper), a[i] + counts[i] , N - counts[i] + A - a[i])
      output[['Model1']][i,] <- d
    }
  # order-constrained hypothesis: median and CI are based on posterior samples  
  } else {
    for(i in 1:nlevels){
      d <- quantile(model1.samples[,i], c(lower, 0.5, upper))
      output[['Model1']][i,] <- d
    }
  }
  
  # Posterior
  # equality constrained hypothesis: median and CI are calculated analytically
  if(is.null(post.samples)){
    for(i in 1:nlevels){
      d <- qbeta(c(lower, .5, upper), a[i] + counts[i] , N - counts[i] + A - a[i])
      output[['Posterior']][i,] <- d
    }
    # order-constrained hypothesis: median and CI are based on posterior samples  
  } else {
    for(i in 1:nlevels){
      d <- quantile(post.samples[,i], c(lower, 0.5, upper))
      output[['Posterior']][i,] <- d
    }
  }
  return(output)
}






