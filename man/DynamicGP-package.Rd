\name{DynamicGP-package}
\alias{DynamicGP-package}
\docType{package}

\title{GP Models for Large-Scale Dynamic Computer Experiments}

\description{This R package provides three functions for emulating
  dynamic computer experiments. The function \code{svdGP} fits full
  SVD-based GP model which is computationally demanding for large-scale
  dyanmic computer experiments. As is well known, the time complexity of
  fitting a GP model is \eqn{O(N^3)} where \eqn{N} is the number of
  training/design points. Since fitting a common GP model for really
  large \eqn{N} would be computationally burdensome, we fit local
  SVD-based GP models on a sequentially selected small neighborhood set
  for every test inputs. The function \code{knnsvdGP} fits K-nearest
  neighbor SVD-based GP models which selects neighborhood sets based on
  the Euclidean distance with repect to the test points. The function
  \code{lasvdGP} fits local approximate SVD-based GP model using the
  new algorithm proposed by Zhang et al. (2018).

  The lasvdGP is an extension of the local approximate GP (laGP) model
  developed by Gramacy and Lee (2015) for the emulation of large-scale
  scalar valued computer experiments. The neighborhood selection and
  SVD-based GP model fitting algorithm is suitable for
  parallelization. We use the R package "parallel" for this task. The
  parallelization can achieve nearly linear speed since the procedure on
  each test point is independent and identical.}

\author{Ru Zhang \email{heavenmarshal@gmail.com},

  C. Devon Lin \email{devon.lin@queensu.ca},

  Pritam Ranjan \email{pritam.ranjan@gmail.com}}

\references{Gramacy, R. B. and Apley, D. W. (2015) \emph{Local Gaussian
    process approximation for large computer experiments}, Journal of
    Computational and Graphical Statistics 24(2), 561-578.

  Zhang, R., Lin, C. D. and Ranjan, P. (2018) \emph{Local Gaussian
    Process Model for Large-scale Dynamic Computer Experiments},
    arXiv:1611.09488.}

\keyword{GP model}
\keyword{package}
\keyword{parallelization}
\keyword{SVD}
