\name{panelSelectCorr}
\alias{panelSelectCorr}
\title{Selecting a Panel(s) using the Mouse}
\description{
  To visualize relevant parameters for correlation matrix, removing those panels
  related to variances.
}
\usage{
panelSelectCorr(CovPlotData, range.logical.contour = FALSE,
                range.logical.all = TRUE)
}
\arguments{
  \item{CovPlotData}{the output list returned by \code{VisCov}.}
  \item{range.logical.all}{a logical value indicating whether the ranges in the 
  panels are set as the same over different sets of samples. For the histograms,
  the ranges of frequencies are set as the same over the different sets of 
  samples.}
  \item{range.logical.contour}{a logical value indicating whether the ranges of
  the contour plots are set as the same over different sets of samples when 
  \code{range.logical.all} is TRUE.}
}
\value{
  No return value, called for side effects.
}
\details{
  After drawing a graph by \code{VisCov}, one can select and draw relevant 
  panels for correlation matrix.
}
\references{
  Tokuda, T., Goodrich, B., Van Mechelen, I., Gelman, A. and Tuerlinckx, F. 
  (submitted). Visualizing Distributions of Covariance Matrices.
}
\examples{
\donttest{
set.seed(1234)
distribution = "LKJ"
eta = 5
dim = 50
param = list(prob = 0.5, dim = dim, eta=eta, scaleCov = diag(1,dim))
CovPlotData = VisCov(distribution, param,title.logical = FALSE)
panelSelectCorr(CovPlotData)
}
}
