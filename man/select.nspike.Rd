\name{select.nspike}
\alias{select.nspike}
\title{Finding Distant Spikes}
\description{
Estimates the number of distant spikes in the population based on the Generalized Spiked Population model. A finite upper bound (\code{n.spikes.max}) of the number of distant spikes must be provided.}
\usage{
select.nspike(samp.eval, p, n, n.spikes.max, evals.out = FALSE, smooth = TRUE)
}
\arguments{
  \item{samp.eval}{
Numeric vector containing the sample eigenvalues. The vector must have dimension \code{n} or \code{n-1}, it may be unordered.
}
  \item{p}{
The number of features.
}
  \item{n}{
The number of samples.
}
  \item{n.spikes.max}{
Upper bound of the number of distant spikes in the population.
}
  \item{evals.out}{
Logical. If \code{TRUE}, the estimated spikes and non-spikes are returned.
}
  \item{smooth}{
Logical. If \code{TRUE}, kernel smoothing will be performed on the estimated population eigenvalue spectrum. Default is \code{TRUE}.
}
}
\details{
The function searches between \eqn{0} and \code{n.spikes.max} to find out the number of distant spikes in the population. It also estimates both non-spiked and spiked eigenvalues based on the \eqn{\lambda}-estimation method.

The argument \code{smooth} is useful when the user assumes the population spectral distribution to be continuous.
}
\value{
\item{n.spikes}{Estimated number of distant spikes.}
\item{spikes}{If \code{evals.out=TRUE}, estimated distant spikes are returned.}
\item{nonspikes}{If \code{evals.out=TRUE}, estimated non-spikes are returned.}
\item{loss}{If \code{evals.out=TRUE}, L-infinity loss function for the spectrum estimation is returned.}
}
\references{
Dey, R. and Lee, S. (2019).
\emph{Asymptotic properties of principal component analysis and shrinkage-bias adjustment under the generalized spiked population model.}
Journal of Multivariate Analysis, Vol 173, 145-164.
}
\author{
Rounak Dey, \email{deyrnk@umich.edu}
}
\seealso{
\code{\link{hdpc_est},\link{pc_adjust}}
}
\examples{
data(hapmap)
#n = 198, p = 75435 for this data

####################################################
\dontrun{
#If you just want the estimated number of spikes
train.eval<-hapmap$train.eval
n<-hapmap$nSamp
p<-hapmap$nSNP

select.nspike(train.eval,p,n,n.spikes.max=10,evals.out=FALSE)

#If you want the estimated spikes and non-spikes
out<-select.nspike(train.eval,p,n,n.spikes.max=10,evals.out=TRUE)
}
}
\keyword{multivariate}
\keyword{optimize}
