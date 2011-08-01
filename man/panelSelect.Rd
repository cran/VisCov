\name{panelSelect}
\alias{panelSelect}
\title{Selecting a panel(s) by mouse}
\description{
To select (using the mouse) one or more panels from a four-layered plot for a single distribution. The selected panels will be drawn in separate graphics window.
}
\usage{
panelSelect(panel.no, CovPlotData)
}
\arguments{
\item{panel.no}{an integer indicating the number of panels that one wishes to select}
\item{CovPlotData}{the output list returned by \code{VisCov}}
}

\details{
After drawing a graph by \code{VisCov}, one can select and draw the specified number of panels by clicking the mouse on the panels of interest.
}

\references{
Tokuda, T., Goodrich, B.,  Van Mechelen, I., Gelman, A., & Tuerlinckx, F. (submitted). Visualizing Distributions of Covariance Matrices
}

\examples{
\dontrun{
CovPlotData = VisCov()
panelSelect(panel.no = 1, CovPlotData) ## Click once on the plot you want to select
panelSelect(panel.no = 4, CovPlotData) ## Click four times on the plots you want to select}
}



