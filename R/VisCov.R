################################### R code for "VisCov"

##### Main functions
VisCov <- function(distribution = "Inverse Wishart", param = list(prob = 0.5, dim = 4, nu = 5, eta = 1, scaleCov = diag(1, 4)), title = distribution, Nsamples = 1000, Ncontours = 100, logSD = TRUE, histogram.Variance = TRUE, histogram.Correlation = TRUE, histogram.Effective.Variance = TRUE, histogram.Effective.Dependence = TRUE, extreme.regio = "Effective Dependence", title.logical = TRUE) {
  ##### Fundamental constants used by the program
  ## List of distributions
  distribution.list <- c(
    "User defined distribution", "Wishart", "Inverse Wishart", "Scaled Inverse Wishart",
    "Scaled Inverse Wishart for correlation", "Scaled with uniform on correlation",
    "Jeffreys", "Hierarchical Inverse Wishart", "LKJ", "Scaled LKJ"
  )

  ## List of variables to pass for selecting panels
  CovPlotData.list <- c(
    "mat", "distribution", "dim", "ppcex", "Ncontours", "Nsamples", "logSD",
    "HistEffective.Variance.out", "HistEffective.Dependence.out", "extreme.regio.logical", "extreme.regio",
    "size.label", "Contour.out", "num", "Histogram.out", "size.label", "Scatter.out", "Submatrix.out", "ThreeDCor.out",
    "SelectPlot", "layout.matrix", "PlotFunction.list", "histogram.Effective.Variance", "histogram.Effective.Dependence",
    "histogram.Variance", "histogram.Correlation", "upperRegion", "lowerRegion", "numDevice", "range.yaxis.all"
  )

  ## List of functions for plots
  PlotFunction.list <- c(
    "ContourPlot", "HistogramVariance", "HistogramCor", "HistogramEffectiveVariance",
    "HistogramEffectiveDependence", "SubmatrixPlot", "ThreeDcorPlot", "ScatterPlot"
  )

  #### Generating the covariance matrices
  N <- max(Nsamples, Ncontours)
  if (extreme.regio == "Effective Variance" | extreme.regio == "Effective Dependence") {
    extreme.regio.logical <- TRUE
  } else {
    extreme.regio.logical <- FALSE
  }
  if (!(any(distribution == distribution.list))) {
    stop("Distribution must be any of the following: ", paste(distribution.list, ", ", sep = ""), "\n")
  }
  if (distribution == "Jeffreys") {
    message("This may take a while!")
  }
  if (distribution == "User defined distribution") {
    mat <- param$mat
    param$dim <- ncol(mat[[1]])
    dim <- ncol(mat[[1]])

    ## Check the number of matrices
    if (length(mat) < N) {
      stop("NOTE! The number of matrices is not sufficient (less than specified in Ncontours or Nsamples)\n")
    }
    ## Check dimension 1 (should be more than 1)
    if (is.null(dim)) {
      stop("The dimension should be more than one!\n")
    }
    ## Check dimensions 2 (should be the square and the same for all matrices)
    CheckIndex.dimension <- lapply(mat, DimensionCheck, dim = dim)
    sumCheck.dimension <- sum(unlist(CheckIndex.dimension))
    if (sumCheck.dimension != length(mat)) {
      message(unlist(CheckIndex.dimension))
      stop("NOTE! Some of matrices have different dimensions from the first matrix. Check them in the above indicators (The matrix with zero has different dimensions from the first one)\n")
    }
    ## Check positive Definitive
    CheckIndex.positiveDefinite <- lapply(mat, CheckPD)
    sumCheck.positiveDefinite <- sum(unlist(CheckIndex.positiveDefinite))
    if (sumCheck.positiveDefinite != length(mat)) {
      message(unlist(CheckIndex.positiveDefinite))
      stop("NOTE! Some of matrices are not positive definite. Check them in the above indicators (The matrix number with zero is NOT positive difinite)\n")
    }
  }
  if (distribution != "User defined distribution") {
    ## Check dimension
    if (param$dim < 2) {
      stop("The dimension should be more than one!\n")
    }
    mat <- Sample(distribution, N, param)
  }

  ## Getting correlation matrices
  cormat <- lapply(mat, cov2cor)

  ##### Keeping the current setting for the graphic
  def.par <- par(no.readonly = TRUE)
  on.exit(par(def.par))

  #### Graphic settings
  numDevice <- as.numeric(dev.cur())
  SelectPlot <- NULL
  size.label <- 0.9
  ppcex <- 0.6 ## Pointsize
  vec <- c(2, 3, 0, 0, 0, 8, 9, 10, 11, 12, 1, 7, 0, 0, 0, 4, 5, 6, 0, 0)
  layout.matrix <- t(matrix(vec, 4, 5, byrow = TRUE))
  layout(layout.matrix)
  if (!title.logical) {
    par(omi = c(0, 0, 0.1, 0), mar = c(3, 0, 0, 0), pty = "s")
  } else {
    par(omi = c(0, 0, 0.5, 0), mar = c(3, 0, 0, 0), pty = "s")
  }
  upperRegion <- 0.95 ## probability for upper regions of effective variance or dependence
  lowerRegion <- 0.05 ## probability for lower regions of effective variance or dependence



  #### Transforming the data (mat) for drawing graphs
  Histogram.out <- Histogram(param, mat[1:Nsamples])
  if (logSD) {
    Histogram.out$SD <- log(Histogram.out$SD)
    exp.lab.SD <- expression(log(sigma[1]))
  }
  MainTitle <- GetMainTitle(distribution, distribution.list, param)
  HistEffective.Variance.out <- Histogram.Effective.Variance(param, mat[1:Nsamples])
  HistEffective.Dependence.out <- Histogram.Effective.Dependence(param, mat[1:Nsamples])
  num <- 200 ## Number of lines of equiprobability
  Contour.out <- Contour(param, mat[1:Ncontours], num)
  Scatter.out <- Scatter(param, mat, Nsamples)
  # print(Scatter.out)
  Submatrix.out <- Submatrix(param, cormat, Ncontours, Nsamples)
  ThreeDCor.out <- ThreeDCor(param, cormat)
  dim <- param$dim ## Modification on June 24, 2011

  #### Keep all information for graphs in the form of list
  CovPlotData <- list()
  range.yaxis.all <- list(1) ## dummy
  endp <- length(CovPlotData.list)
  for (i in 1:endp) {
    name <- CovPlotData.list[i]
    CovPlotData[[i]] <- get(name)
  }
  names(CovPlotData) <- CovPlotData.list

  #### Drawing graphs
  for (i in 1:length(PlotFunction.list)) {
    range.yaxis.all[[i]] <- do.call(PlotFunction.list[i], list(CovPlotData))
  }

  if (title.logical) {
    mtext(MainTitle, outer = TRUE, line = 1.2, cex = 1.2)
  }
  ## Save the ranges
  CovPlotData$range.yaxis.all <- range.yaxis.all

  class(CovPlotData) <- "CovPlot"

  ## For selectioning panels
  return(CovPlotData)
}

panelSelect <- function(panel.no, CovPlotData) {
  def.par <- par(no.readonly = TRUE)
  on.exit(par(def.par))

  if (interactive()) {
    m <- panel.no
    if (m > 16) stop("The maximum 16 plots!\n")
    PlotFunction.list <- CovPlotData$PlotFunction.list
    layout.matrix <- CovPlotData$layout.matrix
    no.panel <- length(layout.matrix[layout.matrix > 0])
    ## Adjusting layout.matrix for dim = 2 and 3 (smaller number of panels)
    if (CovPlotData$dim == 2) {
      layout.matrix[2, 3] <- 0
      layout.matrix[2, 4] <- 0
      layout.matrix[2, 5] <- 0
      layout.matrix[3, 2] <- 0
      layout.matrix[4, 3] <- 0
    }
    if (CovPlotData$dim == 3) {
      layout.matrix[2, 5] <- 0
    }
    num <- dev.cur()
    dev.set(num)
    ## Identifying coordinates of the selected panels
    dim.layout <- dim(t(layout.matrix))
    devidedPointsx <- seq(0, 1, length.out = dim.layout[1] + 1)
    devidedPointsy <- seq(0, 1, length.out = dim.layout[2] + 1)
    xdiv <- grconvertX(devidedPointsx, from = "ndc", to = "user")
    ydiv <- grconvertY(devidedPointsy, from = "ndc", to = "user")
    sizex <- length(xdiv)
    sizey <- length(ydiv)
    selected.panels <- rep(NA, panel.no)
    many <- 1
    while (many < m + 1) {
      selected.position <- locator(n = 1)
      x <- selected.position$x
      y <- selected.position$y
      panel.coordinate <- c(sum(xdiv < x), sum(ydiv < y))
      pos <- panel.coordinate
      select.number <- layout.matrix[sizey - pos[2], pos[1]]
      ## Check if the selected one is an actual panel
      if (select.number > 0) {
        selected.panels[many] <- select.number
        many <- many + 1
      }
      ## end of while
    }




    ## Drawing the selected panels
    dev.new()
    side <- ceiling(sqrt(m))
    par(mfrow = c(side, side), pty = "s", omi = c(0.2, 0.2, 0.2, 0.2), mar = c(2.5, 2.5, 1, 1))
    ## Drawing
    for (i in 1:m) {
      par(xaxt = "s", yaxt = "s")
      select <- selected.panels[i]
      CovPlotData$SelectPlot <- select
      if (select > (no.panel - 4)) {
        select <- no.panel - 4
      }
      if (select > 0) {
        do.call(PlotFunction.list[select], list(CovPlotData))
      }
    }
    invisible(dev.set(num))
  }
}

panelSelectMultiple <- function(selected.condition, CovPlotDataMultiple, range.logical.contour = FALSE, range.logical.all = TRUE, row = FALSE) {
  num.datasets <- length(CovPlotDataMultiple)
  num.condition <- length(selected.condition)

  def.par <- par(no.readonly = TRUE)
  on.exit(par(def.par))

  if (row) {
    par(mfcol = c(num.datasets, num.condition), pty = "s", omi = c(0.2, 0.2, 0.2, 0.2), mar = c(2.5, 2.5, 1, 1))
  } else {
    par(mfrow = c(num.condition, num.datasets), pty = "s", omi = c(0.2, 0.2, 0.2, 0.2), mar = c(2.5, 2.5, 1, 1))
  }
  range.matrix <- matrix(nrow = 8, ncol = 2)
  for (i in 1:num.condition) {
    select2 <- getSelect(selected.condition[i])
    ## For having the same ranges across different datasets
    contour.no <- 1
    if (range.logical.contour == TRUE) {
      contour.no <- 0
    }
    if (select2 != 6 & select2 != 7 & (select2 > contour.no) & range.logical.all == TRUE) {
      if (select2 > 8) select2 <- 8
      range.get.condition <- list()
      for (k in 1:length(CovPlotDataMultiple)) {
        temp <- CovPlotDataMultiple[[k]]
        range.yaxis.all <- temp$range.yaxis.all
        range.get.condition[[k]] <- range.yaxis.all[[select2]]
      }
      range.unlist <- unlist(range.get.condition)
      range.yaxis <- c(min(range.unlist), max(range.unlist))
    } else {
      range.yaxis <- NULL
    }
    ## For drawing plots
    for (j in 1:num.datasets) {
      select <- getSelect(selected.condition[i])
      CovPlotData <- CovPlotDataMultiple[[j]]
      CovPlotData$SelectPlot <- select
      if (select > 8) select <- 8
      PlotFunction.list <- CovPlotData$PlotFunction.list
      do.call(PlotFunction.list[select], list(CovPlotData, range.yaxis = range.yaxis))
    }
  }
}


panelSelectCorr <- function(CovPlotData, range.logical.contour = FALSE, range.logical.all = TRUE) {
  def.par <- par(no.readonly = TRUE)
  on.exit(par(def.par))

  if (!is(CovPlotData, "CovPlot")) {
    stop("CovPlotData is not a VisCov object")
  }

  selected.condition <- c("cor", "scatter4", "scatter5", "contour", "threeD", "Effective.Dependence", "Effective.Dependence.submatrix")
  num.condition <- length(selected.condition)

  dev.new(width = 8, height = 5, unit = "in")
  par(pty = "s", omi = c(0.2, 0.2, 0.2, 0.2), mar = c(2.5, 2.5, 1, 1))
  layout(matrix(c(1, 2, 4, 6, 0, 3, 5, 7), 2, 4, byrow = T))

  for (i in 1:num.condition) {
    select <- getSelect(selected.condition[i])
    CovPlotData$SelectPlot <- select
    if (select > 8) select <- 8

    ## For having the same ranges across different datasets
    contour.no <- 1
    if (range.logical.contour == TRUE) {
      contour.no <- 0
    }
    if (select != 6 & select != 7 & (select > contour.no) & range.logical.all == TRUE) {
      range.get.condition <- CovPlotData$range.yaxis.all[[select]]
      range.yaxis <- c(min(range.get.condition), max(range.get.condition))
    } else {
      range.yaxis <- NULL
    }
    ## For drawing plots
    PlotFunction.list <- CovPlotData$PlotFunction.list
    do.call(PlotFunction.list[select], list(CovPlotData, range.yaxis = range.yaxis))
  }
}

######### The other functions that are used in the main functions
DimensionCheck <- function(matrix, dim) {
  matrixDim <- dim(matrix)
  index <- 0
  if (matrixDim[1] == dim & matrixDim[2] == dim) {
    index <- 1
  }
  return(index)
}

foo <- function(Sigma) ## Calculates 1-determinants of submatrices of correlation matrices
{
  sapply(1:(ncol(Sigma)), FUN = function(i) 1 - det(Sigma[1:i, 1:i, drop = FALSE])^(1 / i))
}

exuptriang <- function(Sigma) ## Extracts and vectorizes upper triangular part
{
  m <- lower.tri(Sigma)
  y <- Sigma[m]
}

bar <- function(simulations, size.label) ## Kernel density plot
{
  L <- nrow(simulations)
  k <- ncol(simulations)
  smoothScatter(rep(2, L), simulations[, 2],
    nbin = c(2, 64), xlim = c(1, k),
    ylim = c(0, 1), xlab = "", ylab = "", axes = FALSE
  )
  axis(1, padj = -0.8)
  axis(2, padj = 0.5)
  op <- par(no.readonly = TRUE)
  dimension.plot <- op$mfrow
  if (dimension.plot[1] + dimension.plot[2] < 5) {
    size.label <- 1
  } else {
    size.label <- 0.65 * size.label
  }
  mtext(side = 1, text = "p (size of leading subR)", line = 1.4, cex = size.label)
  mtext(side = 2, text = expression(1 - abs(subR)^{
    ~ frac(1, p)
  }), line = 1.6, cex = size.label)
  sapply(3:k, FUN = function(i) {
    smoothScatter(rep(i, L), simulations[, i],
      nbin = c(3, 64), add = TRUE
    )
  })
  return(invisible(NULL))
}

ThreeDCor <- function(param, cormat) {
  if (param$dim == 2) {
    ThreeDcor.out <- 1
  } ## dummy
  else {
    cors <- lapply(cormat, exuptriang)
    dim <- param$dim
    CC <- do.call(rbind, cors)
    ThreeDcor.out <- CC[, c(1, 2, dim)]
  }
  return(ThreeDcor.out)
}

ThreeDcorPlot <- function(CovPlotData, range.yaxis = NULL) {
  if (CovPlotData$dim == 2) {
    plot(x = 1, y = 1, type = "n", axes = FALSE, xlab = "", ylab = "")
  } else {
    upperRegion <- CovPlotData$upperRegion
    lowerRegion <- CovPlotData$lowerRegion
    CC <- CovPlotData$ThreeDCor.out
    op <- par(no.readonly = TRUE)
    x <- CC[, 1]
    y <- CC[, 2]
    z <- CC[, 3]
    scatter3 <- scatterplot3d(
      x = x, y = y, z = z,
      highlight.3d = FALSE, color = "grey35", angle = 24, pch = 20, cex.symbols = .5,
      mar = op$mar, cex.axis = .5, x.ticklabs = c(-1, 0, 1), y.ticklabs = c(-1, 0, 1),
      z.ticklabs = c(-1, 0, 1), xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
      lab = c(2, 2, 3), lab.z = 2, xlab = "", ylab = "", zlab = "", y.margin.add = 0.2
    )
    mtext(side = 1, text = expression(rho[1:2]), line = 0.5, adj = 0.3, cex = 0.7)
    mtext(side = 2, text = expression(rho[1:3]), line = 0.5, adj = 0.4, cex = 0.7)

    # Calculate the angle of the Y-axis in 2D (in degrees)
    y_axis_diff <- scatter3$xyz.convert(c(0, 0), c(-1, 1), c(0, 0))

    dx <- y_axis_diff$x[2] - y_axis_diff$x[1]
    dy <- y_axis_diff$y[2] - y_axis_diff$y[1]
    y_axis_angle <- atan2(dy, dx) * 180 / pi

    dims <- par("usr")
    x <- dims[1] + 0.80 * diff(dims[1:2])
    y <- dims[3] + 0.05 * diff(dims[3:4])

    text(x, y, labels = expression(rho[2:3]), srt = y_axis_angle, xpd = TRUE)

    extreme.regio.logical <- CovPlotData$extreme.regio.logical
    if (CovPlotData$extreme.regio == "Effective Variance") {
      base <- CovPlotData$HistEffective.Variance.out
    }
    if (CovPlotData$extreme.regio == "Effective Dependence") {
      base <- CovPlotData$HistEffective.Dependence.out
    }
    if (extreme.regio.logical) {
      low <- subset(CC, base < quantile(base, probs = lowerRegion))
      scatter3$points3d(
        x = low[, 1], y = low[, 2], z = low[, 3],
        cex = .8, pch = 20, col = "blue"
      )
      high <- subset(CC, base > quantile(base, probs = upperRegion))
      scatter3$points3d(
        x = high[, 1], y = high[, 2], z = high[, 3],
        cex = .8, pch = 20, col = rgb(1, 0.35, 0.35)
      )
    }
  }
}

Submatrix <- function(param, cormat, Ncontours, Nsamples) {
  dim <- param$dim
  N <- max(Nsamples, Ncontours)
  subcormat <- lapply(cormat, foo)
  subcormat <- do.call(rbind, subcormat)
  Submatrix.out <- subcormat[1:N, ]
  return(Submatrix.out)
}

SubmatrixPlot <- function(CovPlotData, range.yaxis = NULL) {
  if (CovPlotData$dim == 2) {
    plot(x = 1, y = 1, type = "n", axes = FALSE, xlab = "", ylab = "")
  } else {
    dim <- CovPlotData$dim
    size.label <- CovPlotData$size.label
    upperRegion <- CovPlotData$upperRegion
    lowerRegion <- CovPlotData$lowerRegion
    subcormat <- CovPlotData$Submatrix.out
    bar(subcormat, size.label)
    extreme.regio.logical <- CovPlotData$extreme.regio.logical
    if (CovPlotData$extreme.regio == "Effective Variance") {
      base <- CovPlotData$HistEffective.Variance.out
    }
    if (CovPlotData$extreme.regio == "Effective Dependence") {
      base <- CovPlotData$HistEffective.Dependence.out
    }
    if (extreme.regio.logical) {
      extr.low <- subset(subcormat, base < quantile(base, probs = lowerRegion))
      dl <- dim(extr.low)[1]
      for (j in 1:dl) {
        lines(1:dim, extr.low[j, 1:dim], lwd = .5, pch = 20, col = "blue")
      }
      extr.high <- subset(subcormat, base > quantile(base, probs = upperRegion))
      dl <- dim(extr.high)[1]
      for (j in 1:dl) {
        lines(1:dim, extr.high[j, 1:dim], lwd = .5, pch = 20, col = rgb(1, 0.35, 0.35))
      }
    } else {
      extr.selected <- subcormat[1:CovPlotData$Ncontours, ]
      dl <- dim(extr.selected)[1]
      for (j in 1:dl) {
        lines(1:dim, extr.selected[j, 1:dim], lwd = .5, pch = 20)
      }
    }
  }
}

Scatter <- function(param, mat, Nsamples) {
  k <- param$dim
  N <- Nsamples
  x <- rep(NA, N)
  y <- rep(NA, N)
  xlist <- list(c(2, 2), c(1, 2), c(2, 3), c(2, 3), c(3, 4))
  ylist <- list(c(1, 1), c(1, 1), c(1, 1), c(1, 2), c(1, 2))
  NumPlot <- c(0, 2, 4, 5)
  if (k > 4) {
    k <- 4
  }
  corMat <- unlist(lapply(mat, CovCor, k = k))
  scatter.out <- list()
  for (j in 1:NumPlot[k]) {
    index1 <- xlist[[j]]
    index2 <- ylist[[j]]
    elementx <- k * (index1[1] - 1) + index1[2]
    elementy <- k * (index2[1] - 1) + index2[2]
    x <- corMat[seq(elementx, length.out = N, by = k * k)]
    y <- corMat[seq(elementy, length.out = N, by = k * k)]
    ## save the data
    namex <- sprintf("x%d", j)
    scatter.out[[namex]] <- x
    namey <- sprintf("y%d", j)
    scatter.out[[namey]] <- y
  }
  return(scatter.out)
}

ScatterPlot <- function(CovPlotData, range.yaxis = NULL) {
  layout.matrix <- CovPlotData$layout.matrix
  upperRegion <- CovPlotData$upperRegion
  lowerRegion <- CovPlotData$lowerRegion
  size.label <- CovPlotData$size.label
  extreme.regio <- CovPlotData$extreme.regio
  ppcex <- CovPlotData$ppcex
  extreme.regio.logical <- CovPlotData$extreme.regio.logical
  if (extreme.regio == "Effective Variance") {
    base <- CovPlotData$HistEffective.Variance.out
  }
  if (extreme.regio == "Effective Dependence") {
    base <- CovPlotData$HistEffective.Dependence.out
  }
  Scatter.out <- CovPlotData$Scatter.out
  switch <- rep(1, 5)
  no.panel <- length(layout.matrix[layout.matrix > 0])
  if (!is.null(CovPlotData$SelectPlot)) {
    switch <- rep(0, 5)
    select <- CovPlotData$SelectPlot
    switch[select - (no.panel - 5)] <- 1
  }
  k <- CovPlotData$dim
  N <- CovPlotData$Nsamples
  x <- rep(NA, N)
  y <- rep(NA, N)
  xlist <- list(c(2, 2), c(1, 2), c(2, 3), c(2, 3), c(3, 4))
  ylist <- list(c(1, 1), c(1, 1), c(1, 1), c(1, 2), c(1, 2))
  NumPlot <- c(0, 2, 4, 5)
  if (k > 4) {
    k <- 4
  }

  ## Find the range of log(sigma2)
  if (k == 2) {
    data.sigma2 <- c(
      CovPlotData$Scatter.out[["x1"]], CovPlotData$Scatter.out[["y1"]],
      CovPlotData$Scatter.out[["y2"]]
    )
  }
  if (k > 2) {
    data.sigma2 <- c(
      CovPlotData$Scatter.out[["x1"]], CovPlotData$Scatter.out[["y1"]],
      CovPlotData$Scatter.out[["y2"]], CovPlotData$Scatter.out[["y3"]]
    )
  }
  range.sigma2 <- range(data.sigma2)
  if (!is.null(range.yaxis)) {
    range.sigma2 <- range.yaxis
  }
  range.yaxis <- range.sigma2

  for (j in 1:NumPlot[k]) {
    index1 <- xlist[[j]]
    index2 <- ylist[[j]]
    elementx <- k * (index1[1] - 1) + index1[2]
    elementy <- k * (index2[1] - 1) + index2[2]
    namex <- sprintf("x%d", j)
    x <- Scatter.out[[namex]]
    namey <- sprintf("y%d", j)
    y <- Scatter.out[[namey]]
    xlim <- c(-1, 1)
    ylim <- c(-1, 1)
    xlab <- GetLabel(index1[1], index1[2])
    ylab <- GetLabel(index2[1], index2[2])
    if (index1[1] == index1[2]) {
      xlim <- range.sigma2
      if (index1[1] == 1) {
        xlab <- expression(log(sigma[1]))
      }
      if (index1[1] == 2) {
        xlab <- expression(log(sigma[2]))
      }
    }
    if (index2[1] == index2[2]) {
      ylim <- range.sigma2
      if (index2[1] == 1) {
        ylab <- expression(log(sigma[1]))
      }
      if (index2[1] == 2) {
        ylab <- expression(log(sigma[2]))
      }
    }
    if (switch[j] == 1) {
      if (j == 1) {
        plot(x, y,
          axes = FALSE, xlab = NA, ylab = NA, xlim = xlim,
          ylim = ylim, type = "p", pch = 20, cex = ppcex, asp = 1, col = "grey35"
        )
        if (extreme.regio.logical) {
          xy <- cbind(x, y, base)
          xy2 <- subset(xy, base < quantile(base, probs = lowerRegion))
          points(xy2[, 1:2], col = "blue", pch = 20, cex = ppcex)
          xy2 <- subset(xy, base > quantile(base, probs = upperRegion))
          points(xy2[, 1:2], col = rgb(1, 0.35, 0.35), pch = 20, cex = ppcex)
        }
      }
      if (j > 1) {
        plot(x, y,
          axes = FALSE, xlab = NA, ylab = NA, xlim = xlim,
          ylim = ylim, type = "p", pch = 20, col = "grey35", cex = ppcex, lab = c(2, 2, 4)
        )
        if (extreme.regio.logical) {
          xy <- cbind(x, y, base)
          xy2 <- subset(xy, base < quantile(base, probs = lowerRegion))
          points(xy2[, 1:2], pch = 20, col = "blue", cex = ppcex)
          xy2 <- subset(xy, base > quantile(base, probs = upperRegion))
          points(xy2[, 1:2], pch = 20, col = rgb(1, 0.35, 0.35), cex = ppcex)
        }
      }
      title(xlab = xlab, ylab = ylab, cex.lab = size.label, line = 1.3)
      axis(1, padj = -0.8, c(roundoff(xlim[1], 1, 1), round(mean(xlim),
        digits = 1
      ), roundoff(xlim[2], 1, -1)))
      axis(2, padj = 1, c(roundoff(ylim[1], 1, 1), round(mean(ylim),
        digits = 1
      ), roundoff(ylim[2], 1, -1)))
      box()
      ## end of if (switch[j]==1){
    }
  }
  return(range.yaxis)
}

ScatterPlotCor <- function(CovPlotData, range.yaxis = NULL) {
  layout.matrix <- CovPlotData$layout.matrix
  upperRegion <- CovPlotData$upperRegion
  lowerRegion <- CovPlotData$lowerRegion
  size.label <- CovPlotData$size.label
  extreme.regio <- CovPlotData$extreme.regio
  ppcex <- CovPlotData$ppcex
  extreme.regio.logical <- CovPlotData$extreme.regio.logical
  if (extreme.regio == "Effective Variance") {
    base <- CovPlotData$HistEffective.Variance.out
  }
  if (extreme.regio == "Effective Dependence") {
    base <- CovPlotData$HistEffective.Dependence.out
  }
  Scatter.out <- CovPlotData$Scatter.out
  switch <- rep(1, 5)
  no.panel <- length(layout.matrix[layout.matrix > 0])
  if (!is.null(CovPlotData$SelectPlot)) {
    switch <- rep(0, 5)
    select <- CovPlotData$SelectPlot
    switch[select - (no.panel - 5)] <- 1
  }
  k <- CovPlotData$dim
  N <- CovPlotData$Nsamples
  x <- rep(NA, N)
  y <- rep(NA, N)
  xlist <- list(c(2, 2), c(1, 2), c(2, 3), c(2, 3), c(3, 4))
  ylist <- list(c(1, 1), c(1, 1), c(1, 1), c(1, 2), c(1, 2))
  NumPlot <- c(0, 2, 4, 5)
  if (k > 4) {
    k <- 4
  }

  ## Find the range of log(sigma2)
  if (k == 2) {
    data.sigma2 <- c(
      CovPlotData$Scatter.out[["x1"]], CovPlotData$Scatter.out[["y1"]],
      CovPlotData$Scatter.out[["y2"]]
    )
  }
  if (k > 2) {
    data.sigma2 <- c(
      CovPlotData$Scatter.out[["x1"]], CovPlotData$Scatter.out[["y1"]],
      CovPlotData$Scatter.out[["y2"]], CovPlotData$Scatter.out[["y3"]]
    )
  }
  range.sigma2 <- range(data.sigma2)
  if (!is.null(range.yaxis)) {
    range.sigma2 <- range.yaxis
  }
  range.yaxis <- range.sigma2

  for (j in 1:NumPlot[k]) {
    index1 <- xlist[[j]]
    index2 <- ylist[[j]]
    elementx <- k * (index1[1] - 1) + index1[2]
    elementy <- k * (index2[1] - 1) + index2[2]
    namex <- sprintf("x%d", j)
    x <- Scatter.out[[namex]]
    namey <- sprintf("y%d", j)
    y <- Scatter.out[[namey]]
    xlim <- c(-1, 1)
    ylim <- c(-1, 1)
    xlab <- GetLabel(index1[1], index1[2])
    ylab <- GetLabel(index2[1], index2[2])
    if (index1[1] == index1[2]) {
      xlim <- range.sigma2
      if (index1[1] == 1) {
        xlab <- expression(log(sigma[1]))
      }
      if (index1[1] == 2) {
        xlab <- expression(log(sigma[2]))
      }
    }
    if (index2[1] == index2[2]) {
      ylim <- range.sigma2
      if (index2[1] == 1) {
        ylab <- expression(log(sigma[1]))
      }
      if (index2[1] == 2) {
        ylab <- expression(log(sigma[2]))
      }
    }
    if (switch[j] == 1) {
      if (j == 1) {
        plot(x, y,
          axes = FALSE, xlab = NA, ylab = NA, xlim = xlim,
          ylim = ylim, type = "p", pch = 20, cex = ppcex, asp = 1, col = "grey35"
        )
        if (extreme.regio.logical) {
          xy <- cbind(x, y, base)
          xy2 <- subset(xy, base < quantile(base, probs = lowerRegion))
          points(xy2[, 1:2], col = "blue", pch = 20, cex = ppcex)
          xy2 <- subset(xy, base > quantile(base, probs = upperRegion))
          points(xy2[, 1:2], col = rgb(1, 0.35, 0.35), pch = 20, cex = ppcex)
        }
      }
      if (j > 1) {
        plot(x, y,
          axes = FALSE, xlab = NA, ylab = NA, xlim = xlim,
          ylim = ylim, type = "p", pch = 20, col = "grey35", cex = ppcex, lab = c(2, 2, 4)
        )
        if (extreme.regio.logical) {
          xy <- cbind(x, y, base)
          xy2 <- subset(xy, base < quantile(base, probs = lowerRegion))
          points(xy2[, 1:2], pch = 20, col = "blue", cex = ppcex)
          xy2 <- subset(xy, base > quantile(base, probs = upperRegion))
          points(xy2[, 1:2], pch = 20, col = rgb(1, 0.35, 0.35), cex = ppcex)
        }
      }
      title(xlab = xlab, ylab = ylab, cex.lab = size.label, line = 1.3)
      axis(1, padj = -0.8, c(roundoff(xlim[1], 1, 1), round(mean(xlim),
        digits = 1
      ), roundoff(xlim[2], 1, -1)))
      axis(2, padj = 1, c(roundoff(ylim[1], 1, 1), round(mean(ylim),
        digits = 1
      ), roundoff(ylim[2], 1, -1)))
      box()
      ## end of if (switch[j]==1){
    }
  }
  return(range.yaxis)
}


HistogramEffectiveDependence <- function(CovPlotData, range.yaxis = NULL) {
  upperRegion <- CovPlotData$upperRegion
  lowerRegion <- CovPlotData$lowerRegion
  histogram.Effective.Dependence <- CovPlotData$histogram.Effective.Dependence
  HistEffective.Dependence.out <- CovPlotData$HistEffective.Dependence.out
  HistEffective.Variance.out <- CovPlotData$HistEffective.Variance.out
  size.label <- CovPlotData$size.label
  extreme.regio <- CovPlotData$extreme.regio
  extreme.regio.logical <- CovPlotData$extreme.regio.logical

  if (extreme.regio == "Effective Variance") {
    base <- CovPlotData$HistEffective.Variance.out
  }
  if (extreme.regio == "Effective Dependence") {
    base <- CovPlotData$HistEffective.Dependence.out
  }

  if (histogram.Effective.Dependence) {
    tmp <- hist(HistEffective.Dependence.out, xlab = "", ylab = "", main = "", xlim = c(0, 1), axes = FALSE)
    range.freq <- c(0, max(tmp$counts))
    if (!is.null(range.yaxis)) {
      range.freq <- range.yaxis
    }
    range.yaxis <- range.freq
    ylab <- "Frequency"
    ## Coloring both tails
    if (extreme.regio == "Effective Dependence") {
      tu <- par("usr")
      clip(tu[1], quantile(HistEffective.Dependence.out, lowerRegion), tu[3], tu[4])
      plot(tmp, col = "blue", add = TRUE)
      clip(quantile(HistEffective.Dependence.out, upperRegion), tu[2], tu[3], tu[4])
      plot(tmp, col = rgb(1, 0.35, 0.35), add = TRUE)
    }
  } else {
    dens.depen <- density(HistEffective.Dependence.out)
    plot(dens.depen, xlab = "", ylab = "", main = "", xlim = c(0, 1), axes = FALSE)
    ylab <- "Density"
  }
  xlab <- "Effective Dependence"
  title(xlab = xlab, ylab = ylab, cex.lab = size.label, line = 1.5)
  axis(1, padj = -0.8, c(0, 0.5, 1))
  axis(2, padj = 0.5)
  if (extreme.regio.logical) {
    AddTailrug(data = HistEffective.Dependence.out, base = base, lowerRegion, upperRegion)
  }
  return(range.yaxis)
}

HistogramEffectiveVariance <- function(CovPlotData, range.yaxis = NULL) {
  upperRegion <- CovPlotData$upperRegion
  lowerRegion <- CovPlotData$lowerRegion
  size.label <- CovPlotData$size.label
  HistEffective.Dependence.out <- CovPlotData$HistEffective.Dependence.out
  HistEffective.Variance.out <- CovPlotData$HistEffective.Variance.out
  extreme.regio <- CovPlotData$extreme.regio
  extreme.regio.logical <- CovPlotData$extreme.regio.logical
  histogram.Effective.Variance <- CovPlotData$histogram.Effective.Variance
  extreme.regio <- CovPlotData$extreme.regio
  if (extreme.regio == "Effective Variance") {
    base <- CovPlotData$HistEffective.Variance.out
  }
  if (extreme.regio == "Effective Dependence") {
    base <- CovPlotData$HistEffective.Dependence.out
  }
  if (histogram.Effective.Variance) {
    tmpnoplot <- hist(HistEffective.Variance.out, plot = FALSE)
    range.freq <- c(0, max(tmpnoplot$counts))
    if (!is.null(range.yaxis)) {
      range.freq <- range.yaxis
    }
    range.yaxis <- range.freq
    tmp <- hist(HistEffective.Variance.out, ylim = range.freq, xlab = "", ylab = "", main = "", axes = FALSE)
    ylab <- "Frequency"
    ## Coloring both tails
    if (extreme.regio == "Effective Variance") {
      tu <- par("usr")
      clip(tu[1], quantile(HistEffective.Variance.out, lowerRegion), tu[3], tu[4])
      plot(tmp, col = "blue", add = TRUE)
      LowerVariance <- HistEffective.Variance.out[HistEffective.Variance.out < quantile(HistEffective.Variance.out, lowerRegion)]
      rug(LowerVariance, col = "blue")
      clip(quantile(HistEffective.Variance.out, upperRegion), tu[2], tu[3], tu[4])
      plot(tmp, col = "red", add = TRUE)
      HigerVariance <- HistEffective.Variance.out[HistEffective.Variance.out > quantile(HistEffective.Variance.out, upperRegion)]
      rug(HigerVariance, col = "red")
    }
  } else {
    dens.effect <- density(HistEffective.Variance.out)
    plot(dens.effect, ylab = "", xlab = "", main = "", axes = FALSE)
    ylab <- "Density"
    range.freq <- c(0, max(dens.effect$y))
  }
  xlab <- "Effective Variance"
  title(xlab = xlab, ylab = ylab, cex.lab = size.label, line = 1.5)
  axis(1, padj = -0.8)
  axis(2, padj = 0.5)
  if (extreme.regio.logical) {
    AddTailrug(data = HistEffective.Variance.out, base = base, lowerRegion, upperRegion)
  }
  return(range.yaxis)
}

HistogramVariance <- function(CovPlotData, range.yaxis = NULL) {
  lowerRegion <- CovPlotData$lowerRegion
  upperRegion <- CovPlotData$upperRegion
  histogram.Variance <- CovPlotData$histogram.Variance
  size.label <- CovPlotData$size.label
  ## For xlabel
  if (CovPlotData$logSD == TRUE) {
    exp.lab.SD <- expression(log(sigma[1]))
  } else {
    exp.lab.SD <- expression(sigma[1])
  }
  Histogram.out <- CovPlotData$Histogram.out
  extreme.regio <- CovPlotData$extreme.regio
  extreme.regio.logical <- CovPlotData$extreme.regio.logical
  if (extreme.regio == "Effective Variance") {
    base <- CovPlotData$HistEffective.Variance.out
  }
  if (extreme.regio == "Effective Dependence") {
    base <- CovPlotData$HistEffective.Dependence.out
  }
  if (histogram.Variance) {
    tmpnoplot <- hist(Histogram.out$SD, plot = FALSE)
    range.freq <- c(0, max(tmpnoplot$counts))
    if (!is.null(range.yaxis)) {
      range.freq <- range.yaxis
    }
    range.yaxis <- range.freq
    hist(Histogram.out$SD, ylim = range.freq, xlab = "", ylab = "", main = "", axes = FALSE)
    ylab <- "Frequency"
  } else {
    dens.SD <- density(Histogram.out$SD)
    plot(dens.SD, xlab = exp.lab.SD, main = "", axes = FALSE)
    ylab <- "Density"
  }
  xlab <- exp.lab.SD
  title(xlab = xlab, ylab = ylab, cex.lab = size.label, line = 1.5)
  axis(1, padj = -0.8)
  axis(2, padj = 0.5)
  if (extreme.regio.logical) {
    AddTailrug(data = Histogram.out$SD, base = base, lowerRegion, upperRegion)
  }
  return(range.yaxis)
}

AddTailrug <- function(data, base, lower, upper) {
  rug(data, col = "grey35")
  Lower <- data[base < quantile(base, lower)]
  rug(Lower, col = "blue")
  Higher <- data[base > quantile(base, upper)]
  rug(Higher, col = rgb(1, 0.35, 0.35))
}

HistogramCor <- function(CovPlotData, range.yaxis = NULL) {
  upperRegion <- CovPlotData$upperRegion
  lowerRegion <- CovPlotData$lowerRegion
  histogram.Correlation <- CovPlotData$histogram.Correlation
  Histogram.out <- CovPlotData$Histogram.out
  size.label <- CovPlotData$size.label
  Nsamples <- CovPlotData$Nsamples
  extreme.regio <- CovPlotData$extreme.regio
  extreme.regio.logical <- CovPlotData$extreme.regio.logical
  if (extreme.regio == "Effective Variance") {
    base <- CovPlotData$HistEffective.Variance.out
  }
  if (extreme.regio == "Effective Dependence") {
    base <- CovPlotData$HistEffective.Dependence.out
  }

  ## Do not draw the histogram. Just to get information
  result.hist <- hist(Histogram.out$Correlation, breaks = seq(-1, 1, length = 9), plot = FALSE)
  d <- CovPlotData$dim

  x <- seq(-1, 1, .001)
  if (histogram.Correlation) {
    scale <- length(result.hist$mids) / 2
    den.beta <- Nsamples * 0.5 * dbeta(.5 * (x + 1), shape1 = d / 2, shape2 = d / 2) / scale
    maxfreq <- max(result.hist$counts)
    ylim <- c(0, max(den.beta, maxfreq))
    if (!is.null(range.yaxis)) {
      ylim <- range.yaxis
    }
    range.yaxis <- ylim
    hist(Histogram.out$Correlation, ylim = ylim, breaks = seq(-1, 1, length = 9), xlab = "", ylab = "", main = "", axes = FALSE)
    ylab <- "Frequency"
    ## Draw a reference density curve
    lines(x, den.beta, col = "green")
  } else {
    den.beta <- 0.5 * dbeta(.5 * (x + 1), shape1 = d / 2, shape2 = d / 2)
    dens.rho <- density(Histogram.out$Correlation, from = -1, to = 1)
    minden <- min(min(den.beta), min(dens.rho$y))
    maxden <- max(max(den.beta), max(dens.rho$y))
    plot(dens.rho, xlab = expression(rho[12]), ylab = "", main = "", ylim = c(minden, maxden), axes = FALSE)
    ylab <- "Density"
    lines(x, den.beta, col = "green")
  }
  xlab <- expression(rho[12])
  title(xlab = xlab, ylab = ylab, cex.lab = size.label, line = 1.5)
  axis(1, padj = -0.8, c(-1, 0, 1))
  axis(2, padj = 0.5)
  if (extreme.regio.logical) {
    AddTailrug(data = Histogram.out$Correlation, base = base, lowerRegion, upperRegion)
  }
  return(range.yaxis)
}

Contour <- function(param, mat, num) {
  p <- param$prob
  N <- length(mat)
  radsq <- qchisq(p, 2)
  rad <- sqrt(radsq)
  angle <- (2 * pi / num) * (1:num)
  points <- rad * cbind(cos(angle), sin(angle))
  pointsAll <- lapply(mat, GetPoints, points = points, index = 1)
  pointsAll2 <- unlist(lapply(mat, GetPoints, points = points, index = 2))
  Contour.out <- list()
  Contour.out$pointsAll <- pointsAll
  Contour.out$pointsAll2 <- pointsAll2
  return(Contour.out)
}

ContourPlot <- function(CovPlotData, range.yaxis = NULL) {
  lowerRegion <- CovPlotData$lowerRegion
  upperRegion <- CovPlotData$upperRegion
  size.label <- CovPlotData$size.label
  extreme.regio.logical <- CovPlotData$extreme.regio.logical
  extreme.regio <- CovPlotData$extreme.regio
  if (extreme.regio == "Effective Variance") {
    base <- CovPlotData$HistEffective.Variance.out
  }
  if (extreme.regio == "Effective Dependence") {
    base <- CovPlotData$HistEffective.Dependence.out
  }
  N <- CovPlotData$Ncontours
  num <- CovPlotData$num
  Contour.out <- CovPlotData$Contour.out
  pointsAll <- Contour.out$pointsAll
  pointsAll2 <- Contour.out$pointsAll2
  ind.sel1 <- seq(1, length.out = N * num, by = 2)
  ind.sel2 <- seq(2, length.out = N * num, by = 2)
  xlim <- range(pointsAll2[ind.sel1])
  ylim <- range(pointsAll2[ind.sel2])
  xlim <- c(min(xlim[1], ylim[1]), max(xlim[2], ylim[2]))
  if (!is.null(range.yaxis)) {
    xlim <- range.yaxis
  }
  range.yaxis <- xlim
  ylim <- xlim
  plot(
    axes = FALSE, x = 0, y = 0, type = "n",
    xlim = xlim, ylim = ylim, asp = 1, xlab = NA, ylab = NA, pty = "s"
  )
  axis(1, padj = -0.8, c(
    roundoff(xlim[1], 1, 1), round(mean(xlim), digits = 1),
    roundoff(xlim[2], 1, -1)
  ))
  axis(2, padj = 0.5, c(
    roundoff(ylim[1], 1, 1), round(mean(ylim), digits = 1),
    roundoff(ylim[2], 1, -1)
  ))
  box()
  xlab <- expression(X[1])
  ylab <- expression(X[2])
  title(xlab = xlab, ylab = ylab, cex.lab = size.label, line = 1.5)
  lapply(pointsAll, points, col = "grey35", type = "l", lwd = 0.5)
  if (extreme.regio.logical) {
    pointsAllbase <- list()
    many <- 1
    for (s in 1:length(pointsAll)) {
      if (base[s] > quantile(base, probs = upperRegion)) {
        pointsAllbase[[many]] <- pointsAll[[s]]
        many <- many + 1
      }
    }
    lapply(pointsAllbase, points, col = rgb(1, 0.35, 0.35), type = "l", lwd = 0.5)
    pointsAllbase <- list()
    many <- 1
    for (s in 1:length(pointsAll)) {
      if (base[s] < quantile(base, probs = lowerRegion)) {
        pointsAllbase[[many]] <- pointsAll[[s]]
        many <- many + 1
      }
    }
    lapply(pointsAllbase, points, col = "blue", type = "l", lwd = 0.5)
  }
  return(range.yaxis)
}

Histogram <- function(param, mat) {
  k <- nrow(mat[[1]])
  N <- length(mat)
  corMat <- unlist(lapply(mat, CovCor, k = k))
  SD <- exp(corMat[seq(1, N * k^2, by = k^2)])
  Correlation <- corMat[seq(2, N * k^2, by = k^2)]
  Histogram.out <- list(SD = SD, Correlation = Correlation)
  return(Histogram.out)
}

Histogram.Effective.Variance <- function(param, mat) {
  k <- nrow(mat[[1]])
  detMat <- unlist(lapply(mat, det))
  return(detMat^(1 / k))
}

Histogram.Effective.Dependence <- function(param, mat) {
  k <- nrow(mat[[1]])
  mat <- lapply(mat, cov2cor)
  detMat <- unlist(lapply(mat, det))
  return(1 - detMat^(1 / k))
}

GetMainTitle <- function(distribution, distribution.list, param) {
  distribution.nr <- which(distribution == distribution.list)
  if ((distribution.nr < 6) & (distribution.nr > 1)) {
    k <- nrow(param$scaleCov)
    name <- sprintf("%s (with k=%d, degrees of freedom ", distribution, k)
    c <- param$nu - k
    MainTitle <- bquote(paste(.(name), nu == k + .(c), ")", sep = ""))
  } else {
    name <- paste(distribution, " (with dimension k = ", sep = "")
    dim <- param$dim
    MainTitle <- bquote(paste(.(name), .(dim), ")", sep = ""))
  }
  return(MainTitle)
}

CovCor <- function(mat, k) {
  cov <- mat[1:k, 1:k]
  cor <- cov2cor(cov)
  diag(cor) <- log(diag(cov)) / 2
  return(cor)
}

roundoff <- function(a, digits, up) {
  if (up == -1) {
    add <- 10^{
      -digits
    }
    b <- round(a, digits = digits)
    if (a >= b) {
      x <- b
    }
    if (a < b) {
      x <- b - add
    }
  }
  if (up == 1) {
    add <- 10^{
      -digits
    }
    b <- round(a, digits = digits)
    if (a <= b) {
      x <- b
    }
    if (a > b) {
      x <- b + add
    }
  }
  return(x)
}

GetPoints <- function(mat, points, index) {
  cD <- chol(mat[1:2, 1:2])
  cDlow <- t(cD)
  if (index == 1) {
    points2 <- t(cDlow %*% t(points))
  } else {
    points2 <- cDlow %*% t(points)
  }
  return(points2)
}

GetLabel <- function(i, j) {
  label <- NA
  if (i == 1 & j == 2) {
    label <- expression(rho[12])
  }
  if (i == 2 & j == 3) {
    label <- expression(rho[23])
  }
  if (i == 3 & j == 4) {
    label <- expression(rho[34])
  }
  return(label)
}

dJeffreys <- function(D) {
  m <- nrow(D)
  det <- det(D)
  return(det^(-0.5 * (m + 1)))
}

Sample <- function(distribution, N, param) {
  name <- paste("r", gsub(" ", "", distribution), sep = "")
  dummy <- as.list(rep(NA, N))
  mat <- lapply(dummy, name, param = param)
  return(mat)
}

rInverseWishart <- function(dummy, param) {
  nu <- param$nu
  Q <- param$scaleCov
  mat <- rwishart(nu, solve(Q))$IW
  return(mat)
}

rWishart <- function(dummy, param) {
  nu <- param$nu
  Q <- param$scaleCov
  mat <- rwishart(nu, Q)$W
  return(mat)
}

rHierarchicalInverseWishart <- function(dummy, param) {
  k <- param$dim
  scale <- param$scale
  upper <- param$upper
  nu <- k + 2 + rexp(1, 1 / scale)
  kappa <- (nu - k - 1) * upper * runif(1, 0, 1)
  Q <- diag(kappa, k)
  mat <- rwishart(nu, solve(Q))$IW
  return(mat)
}

rJeffreys <- function(dummy, param) {
  range.r <- param$range.r
  range.v <- param$range.v
  m <- param$dim
  rmax <- max(abs(range.r))
  vmin <- min(range.v)
  R2 <- matrix(rep(rmax, m * m), m, m)
  diag(R2) <- rep(1, m)
  S2 <- sqrt(diag(vmin, m))
  mat2 <- S2 %*% R2 %*% S2
  max.dens.Jeffreys <- dJeffreys(mat2)
  indexJef <- 0
  while (indexJef == 0) {
    R <- rcorrmatrix(m)
    S <- sqrt(diag(runif(m, range.v[1], range.v[2])))
    mat <- S %*% R %*% S
    ratio <- dJeffreys(mat) / max.dens.Jeffreys
    ran <- runif(1, 0, 1)
    if (ran < ratio) {
      indexJef <- 1
      message("Yes, another Jeffreys realization! \n")
    }
  }
  return(mat)
}

CheckPD <- function(R) {
  m <- nrow(R)
  eigenv <- eigen(R)$values
  indexPosi <- 1
  if (any(eigenv <= 0)) indexPosi <- 0
  return(indexPosi)
}

rScaledInverseWishart <- function(dummy, param) {
  nu <- param$nu
  Q <- param$scaleCov
  mu0 <- param$mu0
  s0 <- param$s0
  m <- nrow(Q)
  D <- rInverseWishart(dummy, param)
  qq <- pnorm(0, mu0, s0)
  d.u <- runif(m, qq, 1)
  d <- qnorm(d.u, mu0, s0)
  matDelta <- diag(d)
  mat <- matDelta %*% D %*% matDelta
  return(mat)
}

rScaledInverseWishartforcorrelation <- function(dummy, param) {
  nu <- param$nu
  Q <- param$scaleCov
  mu0 <- param$mu0
  s0 <- param$s0
  m <- nrow(Q)
  D <- rInverseWishart(dummy, param)
  R <- cov2cor(D)
  qq <- pnorm(0, mu0, s0)
  d.u <- runif(m, qq, 1)
  d <- qnorm(d.u, mu0, s0)
  matDelta <- diag(d)
  mat <- matDelta %*% R %*% matDelta
  return(mat)
}

rScaledwithuniformoncorrelation <- function(dummy, param) {
  nu <- param$nu
  m <- param$dim
  mu0 <- param$mu0
  s0 <- param$s0
  R <- rcorrmatrix(m)
  qq <- pnorm(0, mu0, s0)
  d.u <- runif(m, qq, 1)
  d <- qnorm(d.u, mu0, s0)
  matDelta <- diag(d)
  mat <- matDelta %*% R %*% matDelta
  return(mat)
}

rLKJ <- function(dummy, param) {
  eta <- param$eta
  dim <- param$dim
  mat <- rlkjcorr(1, dim, eta = eta)
  return(mat)
}

rScaledLKJ <- function(dummy, param) {
  eta <- param$eta
  dim <- param$dim
  mu0 <- param$mu0
  s0 <- param$s0
  R <- rlkjcorr(1, dim, eta = eta)
  qq <- pnorm(0, mu0, s0)
  d.u <- runif(dim, qq, 1)
  d <- qnorm(d.u, mu0, s0)
  matDelta <- diag(d)
  mat <- matDelta %*% R %*% matDelta
  return(mat)
}


getSelect <- function(name) {
  name.list <- c("contour", "vari", "cor", "Effective.Variance", "Effective.Dependence", "Effective.Dependence.submatrix", "threeD", "scatter1", "scatter2", "scatter3", "scatter4", "scatter5")
  selected.number <- which(name.list == name)
  if (length(selected.number) == 0) {
    stop("Wrong specification of conditions!\n")
  }
  return(selected.number)
}

######### Record of modifications
## First edition: Feb 10, 2011
## Second edition: June 20, 2011. Correction of errors:
##   1. Adding dim = param$dim in the function of VisCov
##   2. For default, the setting is dim = 4, param = list(prob = 0.5, dim = 4, nu = 4+1, scaleCov = diag(1,4)), instead of
##      param = list(prob = 0.5, dim = 4, nu = dim+1, scaleCov = diag(1,dim)).
