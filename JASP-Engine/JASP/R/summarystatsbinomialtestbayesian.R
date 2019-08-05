#
# Copyright (C) 2013-2018 University of Amsterdam
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

SummaryStatsBinomialTestBayesian <- function(jaspResults, dataset = NULL, options, ...) {
  
  # Reading in a datafile is not necessary
  # Error checking is not necessary
  
  # Compute the results
  summaryStatsBinomialResults <- .summaryStatsBinomialComputeResults(jaspResults, options)
  
  # Output tables and plots
  .summaryStatsBinomialTableMain(jaspResults, options, summaryStatsBinomialResults)
  .summaryStatsBinomialPlot(     jaspResults, options, summaryStatsBinomialResults)
  
  return()
}

# Results functions ----
.summaryStatsBinomialComputeResults <- function(jaspResults, options) {
  
  # Take results from state if possible
  if (!is.null(jaspResults[["stateSummaryStatsBinomialResults"]])) 
    return(jaspResults[["stateSummaryStatsBinomialResults"]]$object)
  
  # This will be the object that we fill with results
  results        <- list(hypothesisList = list(),
                         binomTable     = list(),
                         binomPlot      = list())
    
  # Run the binomial test
  a         <- options$betaPriorParamA
  b         <- options$betaPriorParamB
  successes <- options$successes
  failures  <- options$failures
  n         <- successes + failures
  theta0    <- options$testValue
  
  # Extract hypothesis
  hypothesisList <- .hypothesisType.summarystats.binomial(hypothesis = options$hypothesis, theta0, bayesFactorType = options$bayesFactorType)
  hypothesis      <- hypothesisList$hypothesis
  
  # Conduct frequentist and Bayesian binomial test
  pValue <- stats::binom.test(x = successes, n = n, p = theta0, alternative = hypothesis)$p.value
  BF10   <- .bayesBinomialTest(counts = successes, n = n, theta0 = theta0, hypothesis = hypothesis, a = a, b = b)
  
  BFlist <- list(BF10    = BF10,
                 BF01    = 1/BF10,
                 LogBF10 = log(BF10))
  
  # Add results to results object
  results[["hypothesisList"]] <- hypothesisList
  results[["binomTable"]] <- list(
    successes = successes,
    failures  = failures,
    theta0    = theta0,
    BF        = BFlist[[options$bayesFactorType]],
    pValue    = pValue
  )
  results[["binomPlot"]] <- list(
    a         = a,
    b         = b,
    successes = successes,
    n         = n,
    theta0    = theta0,
    BF        = BFlist
  )
  
  # Save results to state
  defaultOptions <- c("successes", "failures", "betaPriorParamA", "betaPriorParamB", "testValue", "hypothesis", "bayesFactorType")
  jaspResults[["stateSummaryStatsBinomialResults"]] <- createJaspState(results)
  jaspResults[["stateSummaryStatsBinomialResults"]]$dependOn(defaultOptions)
  
  # Return results object
  return(results)
}

# Main table ----
.summaryStatsBinomialTableMain <- function(jaspResults, options, summaryStatsBinomialResults){
  if (!is.null(jaspResults[["bayesianBinomialTable"]])) return()
  
  tableResults <- summaryStatsBinomialResults[["binomTable"]]
  
  # extract important parameters
  theta0         <- tableResults$theta0
  hypothesisList <- summaryStatsBinomialResults[["hypothesisList"]]
  hypothesis     <- hypothesisList$hypothesis
  
  # create table and state dependencies
  bayesianBinomialTable <- createJaspTable("Bayesian Binomial Test")
  bayesianBinomialTable$dependOn(optionsFromObject = jaspResults[["stateSummaryStatsBinomialResults"]])
  bayesianBinomialTable$position <- 1
  
  # set title for different Bayes factor types
  bfTitle <- hypothesisList$bfTitle

  # set table citations and footnote message for different hypothesis types
  bayesianBinomialTable$addCitation(.summaryStatsCitations[c("Jeffreys1961", "OHagan2004", "Haldane1932")])
  
  message <- hypothesisList$message
  if (!is.null(message)) bayesianBinomialTable$addFootnote(message)
  
  bayesianBinomialTable$addColumnInfo(name = "successes", title = "Successes" , type = "integer")
  bayesianBinomialTable$addColumnInfo(name = "failures" , title = "Failures"  , type = "integer")
  bayesianBinomialTable$addColumnInfo(name = "theta0"   , title = "Test value", type = "number", format = "sf:4;dp:3")
  bayesianBinomialTable$addColumnInfo(name = "BF"       , title = bfTitle     , type = "number", format = "sf:4;dp:3")
  bayesianBinomialTable$addColumnInfo(name = "pValue"   , title = "p"         , type = "number", format = "sf:4;dp:3")
  
  errorMessageTable <- NULL
  
  if (theta0 == 1 && hypothesis == "greater") {
    
    errorMessageTable <- "Cannot test the hypothesis that the test value is greater than 1."
    
  } else if (theta0 == 0 && hypothesis == "less") {
    
    errorMessageTable <- "Cannot test the hypothesis that the test value is less than 0."
  }
  
  jaspResults[["bayesianBinomialTable"]] <- bayesianBinomialTable
  
  if (!is.null(errorMessageTable))
    .quitAnalysis(errorMessageTable)
  
  # extract rows from tableResults
  bayesianBinomialTable$addRows(tableResults)
}

# Prior and Posterior plot ----
.summaryStatsBinomialPlot <- function(jaspResults, options, summaryStatsBinomialResults) {
  
  plotResults    <- summaryStatsBinomialResults[["binomPlot"]]
  hypothesisList <- summaryStatsBinomialResults[["hypothesisList"]]
  hypothesis     <- hypothesisList$hypothesis
  bfSubscripts   <- hypothesisList$bfSubscripts

  # extract parameters needed for prior and posterior plot
  a         <- plotResults$a
  b         <- plotResults$b
  successes <- plotResults$successes
  n         <- plotResults$n
  theta0    <- plotResults$theta0
  BF10      <- plotResults$BF[["BF10"]]
  
  # Prior and posterior plot
  if(options$plotPriorAndPosterior) {
    quantiles       <- .credibleIntervalPlusMedian(credibleIntervalInterval = .95, a, b, successes, n, hyp = hypothesis, theta0 = theta0)
    medianPosterior <- quantiles$ci.median
    CIlower         <- quantiles$ci.lower
    CIupper         <- quantiles$ci.upper
    ppCri           <- c(CIlower, CIupper)
    dfLinesPP       <- .dfLinesPP(dataset = NULL, a = a, b = b, hyp = hypothesis, theta0 = theta0, counts = successes, n = n)
    dfPointsPP      <- .dfPointsPP(dataset = NULL, a = a, b = b, hyp = hypothesis, theta0 = theta0, counts = successes, n = n)
    xName           <- expression(paste("Population proportion ", theta))
    
    if(options$plotPriorAndPosteriorAdditionalInfo){
      p <- JASPgraphs::PlotPriorAndPosterior(dfLines = dfLinesPP, dfPoints = dfPointsPP, xName = xName, BF01 = 1/BF10,
                                             CRI = ppCri, median = medianPosterior, drawCRItxt = TRUE, bfSubscripts = bfSubscripts)
    } 
    else {
      p <- JASPgraphs::PlotPriorAndPosterior(dfLines = dfLinesPP, dfPoints = dfPointsPP, xName = xName, bfSubscripts = bfSubscripts)
    }
    
    # create JASP object
    plot <- createJaspPlot(
      title       = "Prior and Posterior",
      width       = 530,
      height      = 400,
      plot        = p,
      aspectRatio = 0.7
    )
    plot$position <- 2
    plot$dependOn(optionsFromObject = jaspResults[["stateSummaryStatsBinomialResults"]], 
                  options           = c("plotPriorAndPosterior, plotPriorAndPosteriorAdditionalInfo"))
    jaspResults[["priorPosteriorPlot"]] <- plot
  }
}

# helper functions
.hypothesisType.summarystats.binomial <- function(hypothesis_option, theta0, bayesFactorType) {
  if (hypothesis_option == "notEqualToTestValue") {
    
    hypothesis_for_common_functions   <- "twoSided"
    hypothesis                        <- "two.sided"
    message <- paste0("Proportions tested against value: ", theta0, ".")
    
  } else if (hypothesis_option == "greaterThanTestValue") {
    
    hypothesis_for_common_functions   <- "plusSided"
    hypothesis                        <- "greater"
    message <- paste0("For all tests, the alternative hypothesis specifies that the proportion is greater than ", theta0, ".")
    
  } else if (hypothesis_option == "lessThanTestValue") {
    
    hypothesis_for_common_functions   <- "minSided"
    hypothesis                        <- "less"
    message <- paste0("For all tests, the alternative hypothesis specifies that the proportion is less than ", theta0, ".")
    
  }
  
  bfSubscripts <- .setBFsubscripts.summarystats(hypothesis_for_common_functions)
  bfTitle      <- .getBayesfactorTitle.summarystats(bayesFactorType, hypothesis_for_common_functions)
  
  hypothesisList <- list(hypothesis    = hypothesis,
                          message      = message,
                          bfSubscripts = bfSubscripts,
                          bfTitle      = bfTitle)
  
  return(hypothesisList)
}