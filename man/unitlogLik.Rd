% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{unitlogLik}
\alias{unitlogLik}
\title{Extract log likelihood of units of panel models}
\usage{
unitlogLik(object, ...)
}
\arguments{
\item{object}{an object for which log likelihood values for units can be extracted.}

\item{...}{additional arguments.}
}
\value{
Log likelihood extracted for each unit of the panel model \code{object}.

\unitlogLikReturn
}
\description{
\code{unitlogLik()} is a generic function that extracts the log
likelihood for each unit of panel objects returned by panel modeling functions.
While the \code{numeric} value with the log likelihood for the entire panel
is useful and possible via S4 methods \code{logLik()}, the contributions to it
by panel units can be implemented via \code{unitlogLik()}.
}
\details{
This is a generic function: methods can be defined for it.
}
\examples{
# filter, which generates log likelihoods
pfrw <- pfilter(panelRandomWalk(),Np=10)
# extract log likelihood for each panel unit
unitlogLik(pfrw)
}
\seealso{
\link{pfilter}
}
\author{
Carles \Breto
}
\keyword{internal}
