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
  results        <- list(hypothesis = NULL,
                         binomTable = list(),
                         binomPlot  = list())
    
  # Run the binomial test
  a         <- options$betaPriorParamA
  b         <- options$betaPriorParamB
  successes <- options$successes
  failures  <- options$failures
  n         <- successes + failures
  theta0    <- options$testValue
  
  if (options$hypothesis == "notEqualToTestValue") 
    hypothesis <- "two.sided"
  else if (options$hypothesis == "greaterThanTestValue")
    hypothesis <- "greater"
  else
    hypothesis <- "less"
  
  # conduct frequentist and Bayesian binomial test
  pValue <- stats::binom.test(x = successes, n = n, p = theta0, alternative = hypothesis)$p.value
  BF10   <- .bayesBinomialTest(counts = successes, n = n, theta0 = theta0, hypothesis = hypothesis, a = a, b = b)
  
  BFlist <- list(BF10    = BF10,
                 BF01    = 1/BF10,
                 LogBF10 = log(BF10))
  
  # Add results to results object
  results[["hypothesis"]] <- hypothesis
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
  theta0       <- tableResults$theta0
  hypothesis   <- summaryStatsBinomialResults[["hypothesis"]]
  
  # create table and state dependencies
  bayesianBinomialTable <- createJaspTable("Bayesian Binomial Test")
  bayesianBinomialTable$dependOn(optionsFromObject = jaspResults[["stateSummaryStatsBinomialResults"]])
  bayesianBinomialTable$position <- 1
  
  # set title for different Bayes factor types
  if (options$bayesFactorType == "BF01") {
    
    if (tableResults$hypothesis == "two.sided") 
      bf.title <- "BF\u2080\u2081"
    else if (tableResults$hypothesis == "greater")
      bf.title <- "BF\u2080\u208A"
    else if (tableResults$hypothesis == "less")
      bf.title <- "BF\u2080\u208B"
    
  } else if (options$bayesFactorType == "BF10") {
    
    if (hypothesis == "two.sided")
      bf.title <- "BF\u2081\u2080"
    else if (hypothesis == "greater")
      bf.title <- "BF\u208A\u2080"
    else if (hypothesis == "less")
      bf.title <- "BF\u208B\u2080"
    
  } else if (options$bayesFactorType == "LogBF10") {
    
    if (hypothesis == "two.sided")
      bf.title <- "Log(\u0042\u0046\u2081\u2080)"
    else if (hypothesis == "greater")
      bf.title <-"Log(\u0042\u0046\u208A\u2080)"
    else if (hypothesis == "less")
      bf.title <- "Log(\u0042\u0046\u208B\u2080)"
    
  }
  
  # set table citations and footnote message for different hypothesis types
  bayesianBinomialTable$addCitation(.summaryStatsCitations[c("Jeffreys1961", "OHagan2004", "Haldane1932")])
  
  if (hypothesis == "two.sided") {
    
    message <- paste0("Proportions tested against value: ", theta0, ".")
    bayesianBinomialTable$addFootnote(message)
    
  } else if (hypothesis == "greater") {
    
    note <- "For all tests, the alternative hypothesis specifies that the proportion is greater than "
    message <- paste0(note, theta0, ".")
    bayesianBinomialTable$addFootnote(message)
    
  } else if (hypothesis == "less") {
    
    note <- "For all tests, the alternative hypothesis specifies that the proportion is less than "
    message <- paste0(note, theta0, ".")
    bayesianBinomialTable$addFootnote(message)
    
  }
  
  bayesianBinomialTable$addColumnInfo(name = "successes", title = "Successes" , type = "integer")
  bayesianBinomialTable$addColumnInfo(name = "failures" , title = "Failures"  , type = "integer")
  bayesianBinomialTable$addColumnInfo(name = "theta0"   , title = "Test value", type = "number", format = "sf:4;dp:3")
  bayesianBinomialTable$addColumnInfo(name = "BF"       , title = bf.title    , type = "number", format = "sf:4;dp:3")
  bayesianBinomialTable$addColumnInfo(name = "pValue"   , title = "p"         , type = "number", format = "sf:4;dp:3")
  
  # set common error messages in results table
  if (theta0 == 1 && hypothesis == "greater") {
    
    errorMessageTable <- "Cannot test the hypothesis that the test value is greater than 1."
    
  } else if (theta0 == 0 && hypothesis == "less") {
    
    errorMessageTable <- "Cannot test the hypothesis that the test value is less than 0."
  }
  
  jaspResults[["bayesianBinomialTable"]] <- bayesianBinomialTable
  
  # extract rows from summaryStatsBinomialResults
  bayesianBinomialTable$addRows(tableResults)
}

# Prior and Posterior plot ----
.summaryStatsBinomialPlot <- function(jaspResults, options, summaryStatsBinomialResults) {
  
  plotResults <- summaryStatsBinomialResults[["binomPlot"]]
  hypothesis  <- summaryStatsBinomialResults[["hypothesis"]]
  
  if (hypothesis == "two.sided") {
    bfSubscripts <- "BF[1][0]"
  }
  else if (hypothesis == "greater"){
    bfSubscripts <- "BF['+'][0]"
  }
  else if (hypothesis == "less"){
    bfSubscripts <- "BF['-'][0]"
  }
  
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
    dfLinesPP       <- .dfLinesPP( a = a, b = b, hyp = hypothesis, theta0 = theta0, counts = successes, n = n)
    dfPointsPP      <- .dfPointsPP(a = a, b = b, hyp = hypothesis, theta0 = theta0, counts = successes, n = n)
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

.summaryStatsCitations <- c(
  "Jeffreys1961" = "Jeffreys, H. (1961). Theory of Probability. Oxford, Oxford University Press.",
  "OHagan2004"   =   "O'Hagan, A., & Forster, J. (2004). Kendall's advanced theory of statistics vol. 2B: Bayesian inference (2nd ed.). London: Arnold.",
  "Haldane1932"  =   "Haldane, J. B. S. (1932). A note on inverse probability. Mathematical Proceedings of the Cambridge Philosophical Society, 28, 55-61."
)