\name{lasvdgp-package}
\alias{lasvdgp-package}
\docType{package}

\title{Local Approximate SVD-baed GP Models}

\description{
  Local approximate SVD-based GP model (lasvdGP) for large-scale dynamic
  computer experiments. This package intends to address the issue that
  the input set of the dynamic computer experiments is too large (as
  large as more than tens of thousands) for the full GP models to
  emulate. As is well known, the time complexity of fitting a GP model
  is \eqn{O(N^3)} where \eqn{N} is the number of design points. To resolve the
  big \eqn{N} issue, we fit local SVD-based GP models on a small
  neighborhood set of every test inputs. The neighborhood sets are
  selected by the algorithm proposed in Zhang et al. (2017).

  The neighborhood selection and SVD-based GP model fitting algorithm is
  convenient for parallelization. In this package, we provide two ways
  to parallelize the algorithm,
  \itemize{
    \item parallelization via the R package "parallel". It is a
    process-level parallelism.
    \item parallelization via openMP. It is a thread-level parallelism.
  }
  The parallelization can achieve nearly linear speed since the
  procedure on each test point is independent and identical.
}

\author{Ru Zhang \email{14rz14@queensu.ca}}

\references{
  Zhang, R., Lin, C.D. and Ranjan, P. (2017) \emph{Local Gaussian
    Process Model for Large-scale Dynamic Computer Experiments},
  arXiv preprint arXiv:1611.09488.
}

\keyword{GP model}
\keyword{package}
\keyword{parallelization}
\keyword{SVD}
