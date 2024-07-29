library(conflicted)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(readr)
library(readxl)
library(BayesFactor)
library(stringr)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("count", "dplyr")
library(MetaboDiff)
library(tidyHeatmap)
conflict_prefer("heatmap", "tidyHeatmap")
library(RColorBrewer)

# Implementation of quotient implementation
# according to Dieterle et al., 2006; https://www.ncbi.nlm.nih.gov/pubmed/16808434
#
# Jan K
# last update: 2016-03-27
#
quotNorm <- function(X, vars = 1:dim(X)[2], NAerror = F) {
        # x:        data frame to be normalized
        # vars:     index vector of variables o be used, default: all
        # NAerrors: throw error for NA's or just ignore?

        # check if there are any NAs
        if (sum(is.na(X[, vars])) > 0) {
                # throw warning or error?
                if (NAerror) {
                        stop("Data matrix contains NAs")
                } else {
                        warning("Data matrix contains NAs")
                }
        }

        # median reference sample
        ref <- apply(X[, vars], 2, function(x) median(x, na.rm = T))
        # get dilution factors
        d <- apply(X[, vars], 1, function(s) median(as.numeric(s / ref), na.rm = T))
        # apply to each sample  (for each row=sample, divide values by median dilution factor)
        Y <- t(sapply(1:dim(X)[1], function(i) unlist(X[i, ] / d[i])))

        # Y = t(apply(X,1,  function(s) s /  d) )
        rownames(Y) <- rownames(X)

        # return
        list(X = Y, dilution = d)
}

categorize_bayes_factor <- Vectorize(function(bf) {
        factor(if (bf < 0) {
                "H0"
        } else if (bf < 3) {
                "anectdotal evidence"
        } else if (bf < 10) {
                "moderate evidence"
        } else if (bf < 30) {
                "strong evidence"
        } else if (bf < 100) {
                "very\nstrong evidence"
        } else {
                "extreme evidence"
        }, levels = c("extreme evidence", "very\nstrong evidence", "strong evidence", "moderate evidence", "anectdotal evidence", "H0"))
})


detect_group <- Vectorize(function(sample_name) {
        if (str_starts(sample_name, "GF")) {
                "GF"
        } else if (str_starts(sample_name, "11mix")) {
                "11-mix"
        } else if (str_starts(sample_name, "7mix")) {
                "7-mix"
        } else if (str_starts(sample_name, "4mix")) {
                "4-mix"
        } else if (str_starts(sample_name, "10mix")) {
                "10-mix"
        } else {
                stop("invalid sample name")
        }
})


