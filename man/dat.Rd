% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{dat}
\alias{dat}
\title{Overall survival data with treatment crossover and post-discontinuation
therapies}
\format{
\subsection{\code{dat}}{

A data frame with 1000 rows and 12 columns:
\describe{
\item{x1}{A continuous baseline covariate}
\item{x2}{A binary baseline covariate}
\item{a}{The initial treatment}
\item{c}{The censoring time}
\item{t.pfs}{The progression-free survival (PFS) time}
\item{t.co}{The survival time from the start to the end of the
crossover period}
\item{t.os}{The overall survival (OS) time}
\item{delta.pfs}{The event indicator for PFS}
\item{delta.dp}{The disease progression indicator, where \code{delta.dp = 1}
indicates the individual experiences disease progression, and
\code{delta.dp = 0} otherwise}
\item{delta.co}{The treatment crossover indicator}
\item{delta.pdt}{The post-discontinuation
therapy (PDT) indicator}
\item{delta.os}{The event indicator for OS}
}
}
}
\usage{
dat
}
\description{
A subset of data from the World Health Organization Global Tuberculosis
Report ...
}
\keyword{datasets}
