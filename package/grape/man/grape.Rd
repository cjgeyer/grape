\name{grape}
\alias{grape}
\title{MCMC for Markov Graph Model}
\usage{
grape(nnode, nbatch, state, theta, nspac = 1, blen = 1)
}
\description{
   Markov chain Monte Carlo (MCMC) simulation for th Markov graph
   model that is an exponential family with canonical statistics that are
   the number of nodes having degree \eqn{k}
   for \eqn{k = 0, 1, 2, \ldots, n - 1}, where \eqn{n} is \code{nnode}.
}
\arguments{
  \item{nnode}{the number of nodes of the graph.}
  \item{nbatch}{the number of batches for the method of batch means.}
  \item{blen}{the length of the batches (number of Markov chain iterations).}
  \item{nspac}{the spacing of iterations (number of elementary updates per
      Markov chain iteration).}
  \item{state}{an \code{nnode} by \code{nnode} symmetric zero-or-one-valued
      matrix with zeros on the diagonal.  The initial state of the Markov
      chain.  The connectivity matrix of the graph.}
  \item{theta}{the canonical parameter vector.}
}
\details{
  The canonical statistic vector can be computed by
\preformatted{
  fred <- apply(state, 1, sum)
  stat <- tabulate(fred + 1, nnode)
}
}
\value{
  a list containing the input of this function plus the additional
  component \code{batch}, which is an \code{nbatch} by \code{nnode}
  matrix, each row of which is the average of the canonical statistic
  vector over one batch of iterations.
}
\keyword{misc}
