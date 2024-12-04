\name{panelSelect}
\alias{panelSelect}
\title{Selecting a Panel(s) using the Mouse}
\description{
  To select (using the mouse) one or more panels from a four-layered plot for a 
  single distribution. The selected panels will be drawn in a separate graphics 
  window.
}
\usage{
panelSelect(panel.no, CovPlotData)
}
\arguments{
  \item{panel.no}{an integer indicating the number of panels that one wishes to 
  select.}
  \item{CovPlotData}{the output list returned by \code{VisCov}.}
}
\value{
  No return value, called for side effects.
}
\details{
  After drawing a graph by \code{VisCov}, one can select and draw the specified 
  number of panels by clicking the mouse on the panels of interest.
}
\references{
  Tokuda, T., Goodrich, B., Van Mechelen, I., Gelman, A. and Tuerlinckx, F. 
  (submitted). Visualizing Distributions of Covariance Matrices.
}
\examples{
\donttest{
set.seed(1234)
CovPlotData = VisCov()
panelSelect(panel.no = 1, CovPlotData) # Click once on the plot you want to select
panelSelect(panel.no = 4, CovPlotData) # Click four times on the plots you want to select
}
}
