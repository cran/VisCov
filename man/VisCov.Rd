% File src/library/base/man/load.Rd
\name{VisCov}
\alias{VisCov}
\alias{panelSelectMultiple}
\alias{panelSelect}
\alias{CheckPD}
\title{Four-layer plots of distributions of covariance matrices}
\description{
To generate samples from a specific prior distribution, draw the four-layered graphs and return the sampled matrices and the related information on the distribution.
}
\usage{
VisCov(distribution = "Inverse Wishart", param = list(prob = 0.5, dim = 4, nu = 5, scaleCov = diag(1,4)), title = distribution, Ncontours = 100, Nsamples = 1000, logSD = TRUE, histogram.Variance = TRUE, histogram.Correlation = TRUE, histogram.Effective.Variance = TRUE, histogram.Effective.Dependence = TRUE, extreme.regio = "Effective Dependence", title.logical = TRUE) 
panelSelectMultiple(selected.condition, CovPlotDataMultiple, range.logical.contour = FALSE, range.logical.all = TRUE, row = FALSE)
panelSelect(panel.no, CovPlotData)
CheckPD(R) 
}
\arguments{
\item{distribution}{a distribution to be plotted among:}
\item{}{"Inverse Wishart", "Scaled Inverse Wishart","Scaled Inverse Wishart for correlation",
"Scaled with uniform on correlation", "Wishart", and "User defined distribution"}

\item{title.logical}{a logical value indicating whether a titile is given to the graph}
\item{title}{a title of the graph}
\item{Ncontours}{the number of samples for the contour plot}
\item{Nsamples}{the number of samples for the other plots}
\item{logSD}{a logical value indicating whether the standard deviation (sd) is rescaled by the natural logarithm}
\item{histogram.Variance}{a logical value indicating whether the graph of sd is a histogram (TRUE) or
a density (FALSE)}
\item{histogram.Correlation}{a logical value indicating whether the graph of correlation is a histogram (TRUE) or
a density curve (FALSE)}
\item{histogram.Effective.Variance}{a logical value indicating whether the graph of effective variance is a histogram (TRUE)
or a density curve (FALSE)}
\item{histogram.Effective.Dependence}{a logical value indicating whether the graph of effective dependence is a histogram (TRUE) or a density curve (FALSE)}
\item{extreme.regio}{a string of coloring type of samples based on "Effective Dependence" or
"Effective Variance". If neither of them are specified, there is no coloring.}
\item{panel.no}{an integer indicating the number of panels that one wishes to select}
\item{CovPlotData}{a list of covariance matrices and other related information on the panels returned by \code{VisCov}}
\item{selected.condition}{a vector of strings containing types of the graphs that one wishes to select:
For layer 1, "vari", "cor";
for layer 2, "scatter1", "scatter2", "scatter3", "scatter4", "scatter5";
for layer 3, "contour", "threeD";
for layer 4, "Effective.Variance", "Effective.Dependence", "Effective.Dependence.submatrix"}
\item{CovPlotDataMultiple}{a list of objects returned by \code{VisCov}}
\item{range.logical.contour}{a logical value indicating whether the ranges of contour plots are set as the same
over different sets of samples}
\item{range.logical.all}{a logical value indicating whether the ranges of variables are set as the same over
different sets of samples. For the histograms, the ranges of frequencies are
set as the same over the different sets of samples}
\item{row}{a logical value whether the panels belonging to the same set of samples are arranged in the same
row. If it is FALSE, they are arranged in the same column}
\item{R}{a correlation matrix}
\item{param}{a list of parameteres consisting of the following objects:
	    \code{mat}, a list of matrices, applicable only when distribution = "User defined distribution";
            \code{prob}, mass concentration in the contour plot;
	    \code{dim}, dimension of matrices;
	    \code{nu}, degrees of freedom for the inverse Wishart or Wishart type of distributions;
	    \code{scaleCov}, scale matrix for the inverse Wishart or Wishart type of distributions;
	    \code{mu0}, a mean for the folded normal distribution, applicable only when \code{distribution} is "Scaled Inverse Wishart","Scaled Inverse Wishart for correlation" or "Scaled with uniform on correlation";
	    \code{s0}, a standard deviation for the folded normal, applicable only when \code{distribution} is "Scaled Inverse Wishart","Scaled Inverse Wishart for correlation" or "Scaled with uniform on correlation"}
}

\details{
\code{VisCov} generates samples from a specific prior distribution (or user defined distribution), draws the four-layered graphs and returns the sampled matrices and the related information on the distribution.
Using the objects returned by \code{VisCov}, \code{panelSelectMultiple} draws selected panels from several distributions, and
\code{panelSelect} selects panels in a single distribution by mouse and draw these panels. \code{CheckPD} is an internal-use function to examine whether a given correlation matrix is positive-definite. 
}

\references{
Tokuda, T., Goodrich, B.,  Van Mechelen, I., Gelman, A., & Tuerlinckx, F. (submitted). Visualizing Distributions of Covariance Matrices
}

\examples{
###### Example 1: Figures in the paper (Tokuda, Goodrich, Van Mechelen, Gelman, and Tuerlinckx)
## Figure 1
distribution = "Inverse Wishart"
dim = 4
param = list(prob = 0.5, dim = dim, nu = dim+1, scaleCov = diag(1,dim))
CovPlotData = VisCov(distribution, param, title.logical = FALSE)

## Figure 2
distribution = "Inverse Wishart"
dim = 4
param = list(prob = 0.5, dim = dim, nu = dim+1, scaleCov = diag(1,dim))
CovPlotData1 = VisCov(distribution, param, title.logical = FALSE)
distribution = "Inverse Wishart"
dim = 4
param = list(prob = 0.5, dim = dim, nu = dim+50, scaleCov = diag(1,dim))
CovPlotData2 = VisCov(distribution, param, title.logical = FALSE)
selected.condition = c("scatter1", "scatter4", "contour", "Effective.Dependence")
panelSelectMultiple(selected.condition, list(CovPlotData1, CovPlotData2))

## Figure 3
distribution = "Inverse Wishart"
dim = 4
param = list(prob = 0.5, dim = dim, nu = dim+1, scaleCov = diag(1,dim))
CovPlotData1 = VisCov(distribution, param, title.logical = FALSE)
dim = 100
param = list(prob = 0.5, dim = dim, nu = dim+1, scaleCov = diag(1,dim))
CovPlotData2 = VisCov(distribution, param, title.logical = FALSE)
selected.condition = c("Effective.Dependence")
panelSelectMultiple(selected.condition, list(CovPlotData1, CovPlotData2))

## Figure 4
distribution = "Scaled Inverse Wishart for correlation"
dim = 4
param = list(mu0 = 0, s0 = 1, prob = 0.5, dim = dim, nu = dim+1, scaleCov = diag(1,dim))
CovPlotData1 = VisCov(distribution, param, title.logical = FALSE)
param = list(mu0 = 0, s0 = 1,prob = 0.5, dim = dim, nu = dim+50, scaleCov = diag(1,dim))
CovPlotData2 = VisCov(distribution, param, title.logical = FALSE)
selected.condition = c("scatter1", "scatter2", "scatter4", "contour")
panelSelectMultiple(selected.condition, list(CovPlotData1, CovPlotData2))

## Figure 5
distribution = "Scaled with uniform on correlation"
dim = 4
param = list(mu0 = 0, s0 = 1, prob = 0.5, dim = dim)
CovPlotData1 = VisCov(distribution, param, title.logical = FALSE)
dim = 50
param = list(mu0 = 0, s0 = 1, prob = 0.5, dim = dim)
CovPlotData2 = VisCov(distribution, param,title.logical = FALSE)
selected.condition = c("scatter2", "scatter4", "contour","Effective.Dependence.submatrix")
panelSelectMultiple(selected.condition, list(CovPlotData1, CovPlotData2))

## Figure 6
distribution = "Scaled Inverse Wishart"
dim = 4
param = list(mu0 = 0, s0 = 1, prob = 0.5, dim = dim, nu = dim+1, scaleCov = diag(1,dim))
CovPlotData1 = VisCov(distribution, param, title.logical = FALSE)
param = list(mu0 = 0, s0 = 1,prob = 0.5, dim = dim, nu = dim+50, scaleCov = diag(1,dim))
CovPlotData2 = VisCov(distribution, param, title.logical = FALSE)
selected.condition = c("scatter2", "scatter4", "contour","Effective.Dependence")
panelSelectMultiple(selected.condition, list(CovPlotData1,CovPlotData2))

## Figure 7
distribution = "Wishart"
dim = 4
param = list(prob = 0.5, dim = dim, nu = dim+1, scaleCov = diag(1,dim))
CovPlotData1 = VisCov(distribution, param, title.logical = FALSE)
dim = 50
param = list(prob = 0.5, dim = dim, nu = dim+1, scaleCov = diag(1,dim))
CovPlotData2 = VisCov(distribution, param,title.logical = FALSE)
selected.condition = c("scatter2", "scatter4", "contour","Effective.Dependence")
panelSelectMultiple(selected.condition, list(CovPlotData1,CovPlotData2))

## Figure 8
distribution = "User defined distribution"
param = list()
mat = list()
# Generation of covariance matrices
for (i in 1:1000){
index = 0
dim = 50
while (index == 0){
W = matrix(rnorm(dim*dim,0,1),nrow = dim)
svdW= svd(W)
Rho = svdW$u\%*\%t(svdW$v)
D = diag(rbeta(dim,0.5,5))
Sigma = diag(rnorm(dim,0,1))
Sigma = abs(Sigma)
S = Sigma\%*\%Rho\%*\%D\%*\%t(Rho)\%*\%Sigma
index = CheckPD(S)
if (index ==0){print("zero")}
}
mat[[i]] = S
}
param = list(prob = 0.5, mat = mat)
CovPlotData = VisCov(distribution, param , title.logical = FALSE)

##### Example 2: Selection of panels by mouse
\dontrun{
CovPlotData = VisCov()
panelSelect(panel.no = 1, CovPlotData) ## Click once on the plot you want to select
panelSelect(panel.no = 4, CovPlotData) ## Click four times on the plots you want to select
}
}



