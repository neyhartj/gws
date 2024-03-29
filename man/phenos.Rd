% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/package_data.R
\docType{data}
\name{phenos}
\alias{phenos}
\alias{genos}
\alias{map}
\alias{crosses}
\title{Example data}
\format{
A dataset with four objects:

\describe{
  \item{\code{phenos}}{A \code{data.frame} with 245 rows and 5 columns of phenotypic observations of
  four traits on 245 barley lines. Some lines have missing data.}

  \item{\code{genos}}{A \code{matrix} of realized molecular marker genotypes for 245 barley lines
  and 742 SNP markers. Genotypes are coded as \code{1} for homozygous first allele, \code{0} for heterozygous,
  and \code{-1} for homozygous second allele.}

  \item{\code{map}}{A \code{data.frame} with 742 rows and 3 columns containing genetic map positions
  of the 742 SNP markers. The first column contains the marker name, the second column the chromosome,
  and the third column the genetic position.}

  \item{\code{crosses}}{A \code{data.frame} with 186 rows and 2 columns containing specific combinations
  of two parents for prospective breeding crosses.}

 }

An object of class \code{matrix} (inherits from \code{array}) with 245 rows and 742 columns.

An object of class \code{data.frame} with 742 rows and 3 columns.

An object of class \code{data.frame} with 186 rows and 2 columns.
}
\source{
Data are exported from the \code{\link[PopVar]{think_barley.rda}} package.
}
\usage{
phenos

genos

map

crosses
}
\description{
Example phenotypic and molecular marker data on a population of six-row barley lines.
}
\details{
This dataset includes phenotypic observations and realized molecular marker genotypes for 245 barley
lines. Phenotypic observations were recorded for four traits, and the population was genotyped with
742 SNP markers. The data formats are described below:
}
\keyword{datasets}
