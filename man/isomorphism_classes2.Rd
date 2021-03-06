% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toposearch.R
\name{isomorphism_classes2}
\alias{isomorphism_classes2}
\title{Find identical graphs}
\usage{
isomorphism_classes2(igraphlist)
}
\arguments{
\item{igraphlist}{A list with admixture graphs}
}
\value{
An integer vector with isomorphism classes.
Graphs with the same number have identical topology and leaf labels (but may have different internal labels).
}
\description{
Find identical graphs
}
\seealso{
\code{\link{isomorphism_classes}}
}
