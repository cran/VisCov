\name{VisCov}
\alias{VisCov}
\title{Visualizing of Distributions of Covariance Matrices}
\description{
  To generate samples from a specific distribution of covariance matrices, draw 
  the four-layered graphs and return the sampled matrices and the related 
  information on the distribution.
}
\usage{
VisCov(distribution = "Inverse Wishart", param = list(prob = 0.5, dim = 4, 
       nu = 5, eta = 1, scaleCov = diag(1,4)), title = distribution, 
       Nsamples = 1000, Ncontours = 100, logSD = TRUE, 
       histogram.Variance = TRUE, histogram.Correlation = TRUE, 
       histogram.Effective.Variance = TRUE, 
       histogram.Effective.Dependence = TRUE, 
       extreme.regio = "Effective Dependence", title.logical = TRUE)
}
\arguments{
  \item{distribution}{a distribution to be plotted among:\cr\cr 
  "Inverse Wishart", "Scaled Inverse Wishart", "Scaled Inverse Wishart for 
  correlation", "Scaled with uniform on correlation", "Wishart", "LKJ", and 
  "User defined distribution".}
  \item{title.logical}{a logical value indicating whether a title should be 
  given to the graph.}
  \item{title}{title of the graph.}
  \item{Nsamples}{the number of samples for all plots except the contour plot.}
  \item{Ncontours}{the number of samples for the contour plot.}
  \item{logSD}{a logical value indicating whether the natural logarithm of the 
  standard deviation is plotted.}
  \item{histogram.Variance}{a logical value indicating whether the univariate 
  graph of the (log) standard deviation is a histogram (TRUE) or a density 
  (FALSE).}
  \item{histogram.Correlation}{a logical value indicating whether the univariate
  graph of the correlation is a histogram (TRUE) or a density curve (FALSE).}
  \item{histogram.Effective.Variance}{a logical value indicating whether the 
  univariate graph of the effective variance is a histogram (TRUE) or a density 
  curve (FALSE).}
  \item{histogram.Effective.Dependence}{a logical value indicating whether the 
  univariate graph of the effective dependence is a histogram (TRUE) or a 
  density curve (FALSE).}
  \item{extreme.regio}{indicates whether the coloring of extreme samples (in red
  and blue) is based on the "Effective Dependence" or "Effective Variance". If 
  neither of them is specified, there is no coloring.}
  \item{param}{a list of parameters consisting of the following objects:
        \code{mat}, a list of covariance matrices; applicable only when 
        distribution = "User defined distribution"; \code{prob}, mass 
        concentration in the equiprobability contour plot; \code{dim}, dimension
        of matrices; \code{nu}, degrees of freedom for the inverse Wishart or 
        Wishart distributions; \code{scaleCov}, scale matrix for the inverse 
        Wishart or Wishart distributions; \code{mu0} and \code{s0} are the mean
        and the standard deviation for the folded normal distribution, 
        applicable when \code{distribution} is "Scaled Inverse Wishart", 
        "Scaled Inverse Wishart for correlation" or "Scaled with uniform on 
        correlation"; \code{eta}, shape parameter for the LKJ distribution.
  }
}
\value{
  A list of variables to pass for selecting panels.
}
\details{
  \code{VisCov} generates samples from a specific distribution of covariance 
  matrices (or user defined distribution), draws the four-layered graphs and 
  returns the sampled matrices and the related information on the distribution. 

  For more detail, see Tokuda et al. in the following reference.
}
\references{
  Tokuda, T., Goodrich, B., Van Mechelen, I., Gelman, A. and Tuerlinckx, F. 
  (submitted). Visualizing Distributions of Covariance Matrices.
}
\examples{
## Figures in the paper (Tokuda, Goodrich, Van Mechelen, Gelman, and Tuerlinckx)
## Example 1 (Figure 1 in the paper)
set.seed(1234)
distribution = "Inverse Wishart"
dim = 4
param = list(prob = 0.5, dim = dim, nu = dim+1, scaleCov = diag(1,dim))
CovPlotData = VisCov(distribution, param, title.logical = FALSE)

## Example 2  (Figure 2 in the paper)
set.seed(1234)
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

## Example 3  (Figure 3 in the paper)
\donttest{set.seed(1234)
distribution = "Inverse Wishart"
dim = 4
param = list(prob = 0.5, dim = dim, nu = dim+1, scaleCov = diag(1,dim))
CovPlotData1 = VisCov(distribution, param, title.logical = FALSE)
dim = 100
param = list(prob = 0.5, dim = dim, nu = dim+1, scaleCov = diag(1,dim))
CovPlotData2 = VisCov(distribution, param, title.logical = FALSE)
selected.condition = c("Effective.Dependence")
panelSelectMultiple(selected.condition, list(CovPlotData1, CovPlotData2))}

## Example 4  (Figure 4 in the paper)
set.seed(1234)
distribution = "Scaled Inverse Wishart for correlation"
dim = 4
param = list(mu0 = 0, s0 = 1, prob = 0.5, dim = dim, nu = dim+1, scaleCov = diag(1,dim))
CovPlotData1 = VisCov(distribution, param, title.logical = FALSE)
param = list(mu0 = 0, s0 = 1,prob = 0.5, dim = dim, nu = dim+50, scaleCov = diag(1,dim))
CovPlotData2 = VisCov(distribution, param, title.logical = FALSE)
selected.condition = c("scatter1", "scatter2", "scatter4", "contour")
panelSelectMultiple(selected.condition, list(CovPlotData1, CovPlotData2))

## Example 5  (Figure 5 in the paper)
\donttest{set.seed(1234)
distribution = "Scaled with uniform on correlation"
dim = 4
param = list(mu0 = 0, s0 = 1, prob = 0.5, dim = dim)
CovPlotData1 = VisCov(distribution, param, title.logical = FALSE)
dim = 50
param = list(mu0 = 0, s0 = 1, prob = 0.5, dim = dim)
CovPlotData2 = VisCov(distribution, param,title.logical = FALSE)
selected.condition = c("scatter2", "scatter4", "contour", "Effective.Dependence.submatrix")
panelSelectMultiple(selected.condition, list(CovPlotData1, CovPlotData2))}

## Example 6  (Figure 6 in the paper)
\donttest{set.seed(1234)
distribution = "LKJ"
eta = 5
dim = 100
param = list(prob = 0.5, dim = dim, eta=eta, scaleCov = diag(1,dim))
CovPlotData = VisCov(distribution, param,title.logical = FALSE)
panelSelectCorr(CovPlotData)}

## Example 7  (Figure 7 in the paper)
set.seed(1234)
distribution = "Scaled Inverse Wishart"
dim = 4
param = list(mu0 = 0, s0 = 1, prob = 0.5, dim = dim, nu = dim+1, scaleCov = diag(1,dim))
CovPlotData1 = VisCov(distribution, param, title.logical = FALSE)
param = list(mu0 = 0, s0 = 1,prob = 0.5, dim = dim, nu = dim+50, scaleCov = diag(1,dim))
CovPlotData2 = VisCov(distribution, param, title.logical = FALSE)
selected.condition = c("scatter2", "scatter4", "contour", "Effective.Dependence")
panelSelectMultiple(selected.condition, list(CovPlotData1,CovPlotData2))

## Example 8  (Figure 8 in the paper)
\donttest{distribution = "Wishart"
set.seed(1234)
dim = 4
param = list(prob = 0.5, dim = dim, nu = dim+1, scaleCov = diag(1,dim))
CovPlotData1 = VisCov(distribution, param, title.logical = FALSE)
dim = 50
param = list(prob = 0.5, dim = dim, nu = dim+1, scaleCov = diag(1,dim))
CovPlotData2 = VisCov(distribution, param,title.logical = FALSE)
selected.condition = c("scatter2", "scatter4", "contour", "Effective.Dependence")
panelSelectMultiple(selected.condition, list(CovPlotData1, CovPlotData2))}

## Example 9  (Figure 9 in the paper)
\donttest{set.seed(1234)
distribution = "User defined distribution"
param = list()
mat = list()
# Generation of covariance matrices
for (i in 1:1000){
  index = 0
  dim = 50
  while (index == 0){
    W = matrix(rnorm(dim * dim, 0, 1), nrow = dim)
    svdW= svd(W)
    Rho = svdW$u \%*\% t(svdW$v)
    D = diag(rbeta(dim, 0.5, 5))
    Sigma = diag(rnorm(dim, 0, 1))
    Sigma = abs(Sigma)
    S = Sigma \%*\% Rho \%*\% D \%*\% t(Rho) \%*\% Sigma
    # Checking positive definiteness
    eigenv = eigen(S)$values
    if (all(eigenv > 0)) index = 1
    if (index == 0){ print("zero") }
    mat[[i]] = S
  }
}
param = list(prob = 0.5, mat = mat)
CovPlotData = VisCov(distribution, param , title.logical = FALSE)
}
}
