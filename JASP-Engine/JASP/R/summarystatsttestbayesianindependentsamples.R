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

SummaryStatsTTestBayesianIndependentSamples <- function(jaspResults, dataset = NULL, options, ...) {
  
  # Reading in a datafile is not necessary
  # Error checking is not necessary
  
  # Compute the results
  summaryStatsIndSamplesResults <- .summaryStatsIndSamplesComputeResults(jaspResults, options)
  
  # # Output tables and plots
  .summaryStatsIndSamplesTableMain(jaspResults, options, summaryStatsIndSamplesResults)
  .summaryStatsIndSamplesPriorPosteriorPlot(    jaspResults, options, summaryStatsIndSamplesResults)
  # .summaryStatsIndSamplesRobustnessPlot(    jaspResults, options, summaryStatsIndSamplesResults)
  
  return()
}

# Results functions ----
.summaryStatsIndSamplesComputeResults <- function(jaspResults, options) {
  
  # Take results from state if possible
  if (!is.null(jaspResults[["stateSummaryStatsIndSamplesResults"]])) 
    return(jaspResults[["stateSummaryStatsIndSamplesResults"]]$object)
  
  # This will be the object that we fill with results
  results        <- list(hypothesis_list          = list(),
                         ttestTable               = list(),
                         ttestPriorPosteriorPlot  = list(),
                         ttestRobustnessPlot      = list())
  
  # extract hypothesis
  hypothesis_list <- .hypothesisType.summarystats.ttest.independent(options$hypothesis)
  hypothesis      <- hypothesis_list$hypothesis
  
  # Conduct Bayesian independent samples t-test
  ttestResults <- .generalSummaryTtestBF(options = options, paired = FALSE)
  BF10         <- ttestResults$bf
  BFlist       <- list(BF10    = BF10,
                       BF01    = 1/BF10,
                       LogBF10 = log(BF10))
  
  # Add results to results object
  results[["hypothesis_list"]] <- hypothesis_list
  results[["ttestTable"]] <- list(
    t        = options$tStatistic,
    n1       = options$n1Size,
    n2       = options$n2Size,
    BF       = BFlist[[options$bayesFactorType]],
    error    = ttestResults$properror,
    pValue   = ttestResults$pValue[[hypothesis]]
  )
  # results[["ttestPriorPosteriorPlot"]] <- list(
  # t        = options$tStatistic,
  # n1       = options$n1Size,
  # n2       = options$n2Size,
  # oneSided = hypothesis_list[["oneSided"]],
  # BF       = BFlist[[options$bayesFactorType]],
  # BFH1H0   = BFlist[["BF10"]]
  # rscale   = ,
  # delta    =
  # )
  # results[["ttestRobustnessPlot"]] <- list(
  #   a         = a,
  #   b         = b,
  #   successes = successes,
  #   n         = n,
  #   theta0    = theta0,
  #   BF        = BFlist
  # )
  
  # Save results to state
  defaultOptions <- c("tStatistic", "n1Size", "n2Size", "hypothesis", "bayesFactorType", # standard entries
                      "priorWidth", "effectSizeStandardized",                            # default prior
                      "informativeCauchyLocation", "informativeCauchyScale",             # informed cauchy priors
                      "informativeNormalMean", "informativeNormalStd",                   # informed normal priors
                      "informativeTLocation", "informativeTScale", "informativeTDf"      # informed t-distribution
                      )
  jaspResults[["stateSummaryStatsIndSamplesResults"]] <- createJaspState(results)
  jaspResults[["stateSummaryStatsIndSamplesResults"]]$dependOn(defaultOptions)
  
  # Return results object
  return(results)
}

# Main table ----
.summaryStatsIndSamplesTableMain <- function(jaspResults, options, summaryStatsIndSamplesResults){
  if (!is.null(jaspResults[["indSamplesTTestTable"]])) return()
  
  tableResults    <- summaryStatsIndSamplesResults[["ttestTable"]]
  
  # extract important parameters that are dependent on the hypothesis selected
  hypothesis_list <- summaryStatsIndSamplesResults[["hypothesis_list"]]
  hypothesis      <- hypothesis_list$hypothesis
  message         <- hypothesis_list$message
  bf.title        <- .setBFtitle.summarystats(options$bayesFactorType, hypothesis)
  
  # create table and state dependencies
  indSamplesTTestTable <- createJaspTable("Bayesian Independent Samples T-Test")
  indSamplesTTestTable$dependOn(optionsFromObject = jaspResults[["stateSummaryStatsIndSamplesResults"]])
  indSamplesTTestTable$position <- 1
  
  # set footnote messages table citations  
  if (!is.null(message)) indSamplesTTestTable$addFootnote(message)
  
  if (options$effectSizeStandardized == "default") {
    
    indSamplesTTestTable$addCitation(.summaryStatsCitations[c("MoreyRounder2015", "RounderEtAl2009")])
    
  } else if (options$effectSizeStandardized == "informative") {
    
    indSamplesTTestTable$addCitation(.summaryStatsCitations[c("GronauEtAl2017")])
    
  }
  
  indSamplesTTestTable$addColumnInfo(name = "t"        , title = "t"         , type = "number", format = "sf:4;dp:3")
  indSamplesTTestTable$addColumnInfo(name = "n1"       , title = "n\u2081"   , type = "number")
  indSamplesTTestTable$addColumnInfo(name = "n2"       , title = "n\u2082"   , type = "number")
  indSamplesTTestTable$addColumnInfo(name = "BF"       , title = bf.title    , type = "number", format = "sf:4;dp:3")
  indSamplesTTestTable$addColumnInfo(name = "error"    , title = "error %"   , type = "number", format = "sf:4;dp:3")
  indSamplesTTestTable$addColumnInfo(name = "pValue"   , title = "p"         , type = "number", format = "sf:4;dp:3")
  
  jaspResults[["indSamplesTTestTable"]] <- indSamplesTTestTable
  
  # extract rows from summaryStatsBinomialResults
  indSamplesTTestTable$addRows(tableResults)
}

# Prior and Posterior plot ----
.summaryStatsIndSamplesPriorPosteriorPlot <- function(jaspResults, options, summaryStatsBinomialResults) {
  
  plotResults      <- summaryStatsIndSamplesResults[["ttestPriorPosteriorPlot"]]
  hypothesis_list  <- summaryStatsIndSamplesResults[["hypothesis_list"]]
  
  # extract parameters needed for prior and posterior plot
  
  # Prior and posterior plot
  if(options$plotPriorAndPosterior) {
    
      p <- .plotPriorPosterior(
        t                      = plotResults$t,
        n1                     = plotResults$n1,
        n2                     = plotResults$n2,
        paired                 = paired,
        oneSided               = plotResults$oneSided,
        BF                     = plotResults$BF,
        BFH1H0                 = plotResults$BFH1H0,
        rscale                 = rscale,
        delta                  = delta,
        addInformation         = options$plotPriorAndPosteriorAdditionalInfo,
        wilcoxTest             = wilcoxTest,
        options                = options,
        ...
      )
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

# # Robustness plot ----
# .summaryStatsBinomialPlot <- function(jaspResults, options, summaryStatsBinomialResults) {
#   
#   plotResults <- summaryStatsBinomialResults[["binomPlot"]]
#   hypothesis  <- summaryStatsBinomialResults[["hypothesis"]]
#   
#   if (hypothesis == "two.sided") {
#     bfSubscripts <- "BF[1][0]"
#   }
#   else if (hypothesis == "greater"){
#     bfSubscripts <- "BF['+'][0]"
#   }
#   else if (hypothesis == "less"){
#     bfSubscripts <- "BF['-'][0]"
#   }
#   
#   # extract parameters needed for prior and posterior plot
#   a         <- plotResults$a
#   b         <- plotResults$b
#   successes <- plotResults$successes
#   n         <- plotResults$n
#   theta0    <- plotResults$theta0
#   BF10      <- plotResults$BF[["BF10"]]
#   
#   # Prior and posterior plot
#   if(options$plotPriorAndPosterior) {
#     quantiles       <- .credibleIntervalPlusMedian(credibleIntervalInterval = .95, a, b, successes, n, hyp = hypothesis, theta0 = theta0)
#     medianPosterior <- quantiles$ci.median
#     CIlower         <- quantiles$ci.lower
#     CIupper         <- quantiles$ci.upper
#     ppCri           <- c(CIlower, CIupper)
#     dfLinesPP       <- .dfLinesPP( a = a, b = b, hyp = hypothesis, theta0 = theta0, counts = successes, n = n)
#     dfPointsPP      <- .dfPointsPP(a = a, b = b, hyp = hypothesis, theta0 = theta0, counts = successes, n = n)
#     xName           <- expression(paste("Population proportion ", theta))
#     
#     if(options$plotPriorAndPosteriorAdditionalInfo){
#       p <- JASPgraphs::PlotPriorAndPosterior(dfLines = dfLinesPP, dfPoints = dfPointsPP, xName = xName, BF01 = 1/BF10,
#                                              CRI = ppCri, median = medianPosterior, drawCRItxt = TRUE, bfSubscripts = bfSubscripts)
#     } 
#     else {
#       p <- JASPgraphs::PlotPriorAndPosterior(dfLines = dfLinesPP, dfPoints = dfPointsPP, xName = xName, bfSubscripts = bfSubscripts)
#     }
#     
#     # create JASP object
#     plot <- createJaspPlot(
#       title       = "Prior and Posterior",
#       width       = 530,
#       height      = 400,
#       plot        = p,
#       aspectRatio = 0.7
#     )
#     plot$position <- 2
#     plot$dependOn(optionsFromObject = jaspResults[["stateSummaryStatsBinomialResults"]], 
#                   options           = c("plotPriorAndPosterior, plotPriorAndPosteriorAdditionalInfo"))
#     jaspResults[["priorPosteriorPlot"]] <- plot
#   }
# }

# helper functions
.hypothesisType.summarystats.ttest.independent <- function(hypothesis) {
  if (hypothesis == "groupsNotEqual") {
    
    nullInterval <- c(-Inf, Inf)
    hypothesis   <- "twoSided"
    oneSided     <- FALSE
    message      <- NULL
    
  } else if (hypothesis == "groupOneGreater") {
    
    nullInterval <- c(0, Inf)
    hypothesis   <- "plusSided"
    oneSided     <- "right"
    message      <- paste("For all tests, the alternative hypothesis specifies that group 1 is greater than group 2", sep = "")
    
  } else if (hypothesis == "groupTwoGreater") {
    
    nullInterval <- c(-Inf, 0)
    hypothesis   <- "minSided"
    oneSided     <- "left"
    message      <- paste("For all tests, the alternative hypothesis specifies that group 1 is lesser than group 2", sep = "")
  }
  
  return(list(nullInterval = nullInterval,
              hypothesis   = hypothesis,
              oneSided     = oneSided,
              message      = message)
  )
}
.setBFtitle.summarystats <- function(bayesFactorType, hypothesis){
  # set title for different Bayes factor types
  if (bayesFactorType == "BF01") {
    
    if (hypothesis == "twoSided") 
      bf.title <- "BF\u2080\u2081"
    else if (hypothesis == "plusSided")
      bf.title <- "BF\u2080\u208A"
    else if (hypothesis == "minSided")
      bf.title <- "BF\u2080\u208B"
    
  } else if (bayesFactorType == "BF10") {
    
    if (hypothesis == "twoSided")
      bf.title <- "BF\u2081\u2080"
    else if (hypothesis == "plusSided")
      bf.title <- "BF\u208A\u2080"
    else if (hypothesis == "minSided")
      bf.title <- "BF\u208B\u2080"
    
  } else if (bayesFactorType == "LogBF10") {
    
    if (hypothesis == "twoSided")
      bf.title <- "Log(\u0042\u0046\u2081\u2080)"
    else if (hypothesis == "plusSided")
      bf.title <-"Log(\u0042\u0046\u208A\u2080)"
    else if (hypothesis == "minSided")
      bf.title <- "Log(\u0042\u0046\u208B\u2080)"
    
  }
  
  return(bf.title)
  
}

# citations for summary stats module
.summaryStatsCitations <- c(
  "MoreyRounder2015"  = "Morey, R. D., & Rouder, J. N. (2015). BayesFactor (Version 0.9.11-3)[Computer software].",
  "RounderEtAl2009"   = "Rouder, J. N., Speckman, P. L., Sun, D., Morey, R. D., & Iverson, G. (2009). Bayesian t tests for accepting and rejecting the null hypothesis. Psychonomic Bulletin & Review, 16, 225–237.",
  "GronauEtAl2017"    = "Gronau, Q. F., Ly, A., & Wagenmakers, E.-J. (2017). Informed Bayesian T-Tests. The American Statistician.",
  "Jeffreys1961"      = "Jeffreys, H. (1961). Theory of Probability. Oxford, Oxford University Press.",
  "OHaganForster2004" = "O'Hagan, A., & Forster, J. (2004). Kendall's advanced theory of statistics vol. 2B: Bayesian inference (2nd ed.). London: Arnold.",
  "Haldane1932"       = "Haldane, J. B. S. (1932). A note on inverse probability. Mathematical Proceedings of the Cambridge Philosophical Society, 28, 55-61."
)
