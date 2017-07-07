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

MultinomialTestBayesian <- function(dataset = NULL, options, factor, perform = "run",
                                    callback = function(...) 0,  ...){

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

  if(!is.factor(options$factor)){
    options$factor <- factor(options$factor)
  }

  results <- list() # Initialise results object


  # Then, we retrieve the state and initialise the output objects
  state <- .retrieveState()

  bayesianMultinomResults <- NULL # result of the chi-square test
  bayesMultinomTable <- NULL # chi-square table
  bayesMultinomEstimatesTable <- NULL # expected, observed, posterior median and credible intervals
  bayesMultinomEstimatesPlot <- NULL # barplot of factor levels
  bayesMultinomMarginalsPlot <- NULL # list of prior and posterior plots for marginal betas

  # TODO
  # Then, we can fill the output objects with old info if its option did not
  # change.

  # Meta information
  results[["title"]] <- "Bayesian Multinomial Test"
  results[[".meta"]] <- list(list(name = "bayesMultinomTable", type = "table"),
                             list(name = "bayesMultinomEstimatesTable", type = "table"),
                             list(name = "bayesMultinomEstimatesPlot", type = "image"),
                             list(name = "bayesMultinomMarginalPlots", type = "image"))

  # Multinomial Table
  # Generate results
  if (is.null(bayesMultinomResults)){
    bayesMultinomResults <- .bayesMultinomTest(dataset, options, perform)
  }

  results[["bayesMultinomTable"]] <- .bayesMultinomTable(bayesMultinomResults, options, perform)


  # Estimates Table
  if (options[["estimatesTable"]]) {
    # Generate estimates table
    if (is.null(estimatesTable)) {
      estimatesTable <- .bayesMultinomEstimatesTable(bayesMultinomResults, factor, options, perform)
    }

    results[["bayesMultinomEstimatesTable"]] <- estimatesTable

  } else {

    results[["bayesMultinomEstimatesTable"]] <- NULL

  }

  # Multinomial Estimates Plot
  if (options[["estimatesPlot"]]){
    # Generate estimates plot
    if (is.null(esimatesPlot)) {
      estimatesPlot <- .bayesMultinomEstimatesPlot(bayesMultinomResults, options, perform)
    }

    plotPathEstimates <- list(estimatesPlot$data) # for keep later

    results[["estimatesPlot"]] <- estimatesPlot

  } else {

    results[["estimatesPlot"]] <- NULL
    plotPathEstimates <- list()

  }

  # Marginal Plots
  if (options[["marginalPlots"]]){
    # Generate estimates plot
    if (is.null(marginalPlots)) {
      marginalPlots <- .bayesMultinomMarginalPlots(bayesMultinomResults, options, perform)
    }

    plotPathMarginals <- list(marginalPlots$data) # for keep later

    results[["marginalPlots"]] <- marginalPlots

  } else {

    results[["marginalPlots"]] <- NULL
    plotPathMarginals <- list()

  }


  if (perform == "run") {

    state <- list()
    state[["options"]] <- options
    state[["bayesMultinomResults"]] <- bayesMultinomResults
    state[["estimatesTable"]] <- estimatesTable
    state[["estimatesPlot"]] <- estimatesPlot
    state[["marginalPlots"]] <- marginalPlots

    return(list(results=results, status="complete", state=state,
                keep = list(plotPathEstimates, plotPathMarginals)))

  } else {

    return(list(results=results, status="inited", state=state,
                keep = plotPath))

  }

}

# Run bayesian multinomial Test and return Object
.bayesMultinomTest <- function(dataset, options, factor, perform){

  bayesMultinomResults <- NULL

  if(perform == "run" && !is.null(factor)){
    # extract number of categories and counts
    f <- dataset[[.v(factor)]]
    f <- f[!is.na(f)]
    nlev <- nlevels(f)

    # Create table in order of appearance of the dataset
    t <- table(f)
    val <- t[unique(f)]
    if (options$counts != ""){
      val <- dataset[[.v(options$counts)]]
    }

    # calculate Bayes factor
    BF <- .bayesMultinomBF(dataset, options, perform, nlev, val)

    # add marginal credible intervals and median
    ciPlusMedian <- .bayesMultinomCredibleIntervalsPlusMedian(options, nlev, val)

    # create a named list with Bayes factor and individual credible intervals

    bayesMultinomResults <- list(factorName = factor, BF = BF, ciPlusMedian = ciPlusMedian)
  }

  # return the out object
  return(bayesMultinomResults)
}

# Create multinomial Bayes factor table
.bayesMultinomTable(bayesMultinomResults, options, perform){

  table <- list()
  table[["title"]] <- "Bayesian Multinomial Test"

  # set Bayes factor title
  if (options$bayesFactorType == "BF10"){
    bf.title <- "BF\u2081\u2080"
  } else if (options$bayesFactorType == "BF01"){
    bf.title <- "BF\u2080\u2081"
  } else if (options$bayesFactorType == "LogBF10"){
    bf.title <- "Log(\u2009\u0042\u0046\u2081\u2080\u2009)"
  }

  # include fields
  fields <- list(
    list(name="case", title="", type="string", combine=TRUE),
    list(name="BF", title="BF", type="number", format="sf:4;dp:3", title = bf.title)
    )

  # fill in result
  if (!is.null(bayesMultinomResults)){
    table[["data"]] <- list(case = bayesMultinomResults[["factorName"]],
                            BF   = bayesMultinomResults[["BF"]])
  } else {
    # init state?
    data <- list()
    table[["data"]] <- data
  }

  return(table)
}

# TO DO: Create multinomial estimates table
.bayesMultinomEstimatesTable(dataset, options, perform){

      ciInterval <- options$ciInterval

}

# TO DO: Plots for summary

# TO DO: individual plots

# Calculate Savage-Dickey Bayes factor
.bayesMultinomBF <- function(dataset, options, perform, nlevels, counts){

  # extract prior vector alpha and theta
  alpha <- options$prior

  if (options$exProbVar != ""){
    theta <- dataset[[.v(options$exProbVar)]]
  } else {
    theta <- rep(1/nlevels, nlevels)
  }

  # calculate Savage-Dickey Bayes factor
  betaXPlusA <- prod(gamma(counts + a))/gamma(sum(counts + alpha))
  betaA      <- prod(gamma(alpha))/gamma(sum(a))
  priorHeight      <- (1/betaA)  * (prod(theta^(alpha - 1)))
  posteriorHeight  <- (1/betaXPlusA) * (prod(theta^((counts + alpha) - 1)))

  if (options$bayesFactorType == "BF10"){
    BF <- priorHeight/posteriorHeight
  } else if (options$bayesFactorType == "BF01") {
    BF <- posteriorHeight/priorHeight
  } else if (options$bayesFactorType == "LogBF10") {
    BF <- (log(priorHeight/posteriorHeight))
  }

  return(BF)
}

# Calculate individual credible intervals and median
.bayesMultinomCredibleIntervalsPlusMedian <- function(options, nlevels, counts){

  # extract prior vector alpha and theta
  alpha <- options$prior
  ciInterval <- options$credibleIntervalInterval

  # Calculate individual credible intervals and median
  ciPlusMedian <- list()

  for(lev in 1:nlevels){
    lower <- (1 - ciInterval) / 2
    upper <- 1 - lower

    # beta distribution parameters
    a <- alpha[lev] + counts[lev]
    b <- sum(alpha + counts) - a

    # make a list of cis and median
    ciPlusMedian[[lev]] <- qbeta(c(lower, .5, upper), a, b)
    names(ciPlusMedian[[lev]]) <- c("lowerCI, median, upperCI")
  }

  return(ciPlusMedian)

  }


}

# Check function
dataset <- data.frame(id     = 1:5,
                      nation = c("GER", "HUN", "NL", "UK", "USA"),
                      food   = c(34, 82, 90, 45, 67))
options <- list(factor = "nation", counts = "food", a = rep(1, 5), ciInterval = 0.95)
