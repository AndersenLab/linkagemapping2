% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postmapping.R
\name{annotate_lods}
\alias{annotate_lods}
\title{Annotate LOD peaks with variance explained, effect size, and confidence
interval bounds}
\usage{
annotate_lods(lods, cross, annotate_all = FALSE, bayes = FALSE)
}
\arguments{
\item{lods}{A data frame output by the mapping functions to be converted to a
\code{scanone} object}

\item{cross}{The cross object used for the original mapping}

\item{annotate_all}{Boolean whether or not to annotate all markers with
variance explained and effect size. If \code{FALSE} (default), only peak lods
will be annotated.}

\item{bayes}{Boolean whether or not to calculate confidence intervals based
on Bayes statistics (LOD drop 1.5 used to find CI by default)}
}
\value{
The annotated lods data frame with information added for peak markers
of each iteration
}
\description{
Annotate LOD peaks with variance explained, effect size, and confidence
interval bounds
}
