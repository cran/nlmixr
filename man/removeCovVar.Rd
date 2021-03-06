% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/autocovar.R
\name{removeCovVar}
\alias{removeCovVar}
\title{Remove covariate from function string}
\usage{
removeCovVar(fitobject, varName, covariate, categorical = FALSE, isHS = FALSE)
}
\arguments{
\item{fitobject}{an nlmixr 'fit' object}

\item{varName}{a string giving the variable name to which covariate needs to be added}

\item{covariate}{a string giving the covariate name; must be present in the data used for 'fit'}

\item{categorical}{a boolean to represent if the covariate to be added is categorical}

\item{isHS}{a boolean to represent if the covariate to be added is hockey-stick normalized}
}
\value{
returns a list containing the updated model and the parameter names for the covariates added
}
\description{
Function to remove covariates from a given variable's equation in the function string text
}
\author{
Vipul Mann, Matthew Fidler
}
