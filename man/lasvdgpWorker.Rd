\name{lasvdgpWorker}
\alias{lasvdgpWorker}
\alias{lasvdgpms}
\alias{lasvdgpParallel}
\alias{lasvdgpmsParal}
\alias{lasvdgpOmp}
\alias{lasvdgpmsOmp}

\title{
  Local Approximate SVD-baed GP Models
}

\description{
  Fit a local approximate SVD-based GP model with test set \code{X0},
  design set \code{design} and response matrix \code{resp}. The local design
  consists of \code{nn} out of which \code{n0} points are selected by
  the Euclidean distance. If \code{n0} = \code{nn}, it performs the
  naive approach knnsvdGP. The functions \code{lasvdgpWorker} and
  \code{lasvdgpms} are executed in sequential way. The functions
  \code{lasvdgpParallel} and \code{lasvdgpmsParal} parallelize the
  execution via the R package "parallel". The functions
  \code{lasvdgpOmp} and \code{lasvdgpmsOmp} parallelize the execution by
  openMP.}

\usage{
lasvdgpWorker(X0, design, resp, n0, nn,
              nfea = min(1000,nrow(design)),
              nsvd = nn, nadd = 1,
              frac = .9, gstart = 0.001,
              resvdThres = min(5, nn-n0),
              every = min(5,nn-n0),centralize=FALSE,
              maxit=100, verb=0)
lasvdgpms(X0, design, resp, n0, nn,
          nfea = min(1000,nrow(design)),
          nsvd = nn, nadd = 1,
          frac = .9, gstart = 0.001,
          resvdThres = min(5, nn-n0),
          every = min(5,nn-n0),
          nstarts = 5,centralize=FALSE,
          maxit=100, verb=0)
lasvdgpParallel(X0, design, resp, n0, nn,
                nfea = min(1000,nrow(design)),
                nsvd = nn, nadd = 1,
                frac = .9, gstart = 0.001,
                resvdThres = min(5, nn-n0),
                every = min(5,nn-n0),centralize=FALSE,
                maxit=100, verb=0, nthread = 4, clutype="FORK")
lasvdgpmsParal(X0, design, resp, n0, nn,
               nfea = min(1000,nrow(design)),
               nsvd = nn, nadd = 1,
               frac = .9, gstart = 0.001,
               resvdThres = min(2, nn-n0),
               every = min(5,nn-n0),
               nstarts = 5,centralize=FALSE,
               maxit=100, verb=0,
               nthread = 4, clutype="FORK")
lasvdgpOmp(X0, design, resp, n0, nn,
           nfea = min(1000,nrow(design)),
           nsvd = nn, nadd = 1,
           frac = .9, gstart = 0.001,
           resvdThres = min(5, nn-n0),
           every = min(5,nn-n0),centralize=FALSE,
           maxit=100, verb=0,nthread=4)
lasvdgpmsOmp(X0, design, resp, n0, nn,
             nfea = min(1000,nrow(design)),
             nsvd = nn, nadd = 1,
             frac = .9, gstart = 0.001,
             resvdThres = min(5, nn-n0),
             every = min(5,nn-n0),nstarts=5,
             centralize=FALSE, maxit=100, verb=0,
             nthread=4)
}
\arguments{
  \item{X0}{
    An \eqn{M} by \eqn{d} matrix of test inputs, where \eqn{M} is the number of test points
    and \eqn{d} is the dimension of input. The neighborhood will be search and
    the SVD-based GP models will be fitted on every point (row) of
    \code{X0}.
  }
  \item{design}{
    An \eqn{N} by \eqn{d} matrix of design inputs, where \eqn{N} is the number of design
    points. The neighborhood points will be selected out of the points
    (rows) in \code{design}.
  }
  \item{resp}{
    An \eqn{L} by \eqn{N} response matrix of \code{design}, where \eqn{L} is the length of
    the time series outputs, \eqn{N} is the number of design points.
  }
  \item{n0}{
    The number of points in the initial neighborhood set. The initial
    neighborhood set is selected by the Euclidean distance.
  }
  \item{nn}{
    The total number of neighborhood points. The \code{nn-n0} points are
    selected sequentially by the proposed algorithm.
  }
  \item{nfea}{
    The number of feasible points within which to select the
    neighborhood points. This function will only consider the
    \code{nfea} design points closest to the test point in terms of
    Euclidean distance when selecting neighborhood points.
  }
  \item{nsvd}{
    The number of design points closet to the test points on whose
    response matrix to perform the initial singular value decomposition.
  }
  \item{nadd}{
    The number of neighborhood points selected at one iteration.
  }
  \item{frac}{
    The threshold in the cumulative percentage criterion to select the
    number of SVD bases.
  }
  \item{gstart}{
    The starting number and upper bound of for estimating the nugget
    parameter. If \code{gstart = sqrt(.Machine$double.eps)}, the nugget
    will be fixed at \code{sqrt(.Machine$double.eps)}, since the lower
    bound is prefixed at \code{sqrt(.Machine$double.eps)}.
  }
  \item{resvdThres}{
    The threshold to re-perform SVD. After every \code{resvdThres}
    points have been included into the neighborhood set, the SVD of the
    response matrix will be re-performed and the SVD-based GP model will
    be refitted.
  }
  \item{every}{
    The threshold to refit GP models without re-perform SVD. After every
    \code{every} points have been included into the neighborhood set,
    the GP models will be refitted. But the SVD will not be
    re-performed. It is suggested \code{every} <= \code{resvdThres}.
  }
  \item{nstart}{
    The number of starting points used in the numerical maximization of
    the posterior density function. The larger \code{nstart} will
    typically lead to more accurate prediction but longer computational
    time.
  }
  \item{centralize}{
    If \code{centralize=TRUE} the response matrix will be centralized
    (subtract the mean) before the start of the algorithm. The mean will
    be added to the predictive mean at the finish of the algorithm.
  }
  \item{maxit}{
    Maximum number of iterations in the numerical optimization algorithm
    for maximizing the posterior density function.
  }
  \item{verb}{
    A nonnegative integer indicates the level of printing on the
    screen. If \code{verb=0} the function is executed in silence.
  }
  \item{nthread}{
    The number of threads (processes) used in parallel execution of this
    function.
  }
  \item{clutype}{
    The type of cluster in the R package "parallel" to perform
    parallelization.
  }
}
\value{
  \item{pmean}{
    An \eqn{L} by \eqn{M} matrix of posterior predictive mean for the response at
    the test set \code{X0}
  }
  \item{ps2}{
    An \eqn{L} by \eqn{M} matrix of posterior predictive variance for the response at
    the test set \code{X0}
  }
}
\author{Ru Zhang \email{14rz14@queensu.ca}}

\examples{
## load the R package of simulation functions
library("simfuncs")
library("lhs")

timepoints <- seq(0,1,len=200)
design <- lhs::randomLHS(10000,3)
test <- lhs::randomLHS(2000,3)
## evaluate the response matrix on the design matrix
resp <- apply(design,forretal,1,timepoints)

n0 <- 15
nn <- 30
gs <- sqrt(.Machine$double.eps)
## knnsvdGP approach with neighborhood size nn, openMP parallelization
retnn <- lasvdgpOmp(test,design,resp,nn,nn,frac=.95,gstart=gs,
                    centralize=TRUE,nthread=4)
## lasvdGP approach with neighborhood size nn and initial neighborhood size n0,
## openMP parallelization
retgg <- lasvdgpOmp(test,design,resp,n0,nn,frac=.95,gstart=gs,
                    centralize=TRUE,nthread=4)
## knnsvdGP using R package "parallel" for parallelization
retnnp <- lasvdgpParallel(test,design,resp,nn,nn,frac=.95,gstart=gs,
                          centralize=TRUE,nthread=4,clutype="PSOCK")
## lasvdGP using R package "parallel" for parallelization
retggp <- lasvdgpParallel(test,design,resp,n0,nn,frac=.95,gstart=gs,
                          centralize=TRUE,nthread=4,clutype="PSOCK")
}
\keyword{GP model}
\keyword{SVD}
\keyword{neighborhood}
\keyword{prediction}
