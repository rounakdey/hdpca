\name{hdpc_est}
\alias{hdpc_est}
\title{
High-dimensional PCA estimation
}
\description{
Estimates the population eigenvalues, angles between the sample and population eigenvectors, correlations between the sample and population PC scores, and the asymptotic shrinkage factors. Three different estimation methods can be used.
}
\usage{
hdpc_est(samp.eval, p, n, method = c("d.gsp", "l.gsp", "osp"), 
n.spikes, n.spikes.max, n.spikes.out, nonspikes.out = FALSE, smooth = TRUE)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{method}{
String specifying the estimation method. Possible values are "\code{d.gsp}" (default),"\code{l.gsp}" and "\code{osp}".
}
  \item{n.spikes}{
Number of distant spikes in the population (Optional).
}
  \item{n.spikes.max}{
Upper bound of the number of distant spikes in the population. Optional, but needed if \code{n.spikes} is not specified. Ignored if \code{n.spikes} is specified.
}
  \item{n.spikes.out}{
Number of distant spikes to be returned in the output (Optional). If not specified, all the estimated distant spikes are returned.
}
  \item{nonspikes.out}{
Logical. If \code{TRUE} and \code{method="l.gsp"}, the estimated set of non-spikes are returned. If \code{TRUE} and \code{method="osp"}, the estimated value of the non-spike is returned.
}
  \item{smooth}{
Logical. If \code{TRUE} and \code{method="l.gsp"}, kernel smoothing will be performed on the estimated population eigenvalue spectrum. Default is \code{TRUE}.
}
}
\details{
The different choices for \code{method} are:
\itemize{
\item "\code{d.gsp}": \eqn{d}-estimation method based on the Generalized Spiked Population (GSP) model.
\item "\code{l.gsp}": \eqn{\lambda}-estimation method based on the GSP model.
\item "\code{osp}": Estimation method based on the Ordinary Spiked Population (OSP) model.
}
At least one of \code{n.spikes} and \code{n.spikes.max} must be provided. If \code{n.spikes} is provided then \code{n.spikes.max} is ignored, else \code{n.spikes.max} is used to find out the number of distant spikes using \code{\link{select.nspike}}.

The argument \code{nonspikes.out} is ignored if \code{method="d.gsp"}.

The argument \code{smooth} is useful when the user assumes the population spectral distribution to be continuous.
}
\value{
\item{spikes}{An array of estimated distant spikes. If \code{n.spikes.out} is specified, only largest \code{n.spikes.out} many eigenvalues are returned.}
\item{n.spikes}{Number of distant spikes. If \code{n.spikes} is not provided, then the estimated value is returned.}
\item{angles}{An array of estimated cosines of angles between the sample and population eigenvectors corresponding to the distant spikes. The \eqn{k^{th}} element of the array is the estimated cosine of the angle between \eqn{k^{th}} sample and population eigenvectors. If \code{n.spikes.out} is specified, only first \code{n.spikes.out} many \eqn{\cos}(angle)-s are returned.}
\item{correlations}{An array of estimated correlations between the sample and population PC scores corresponding to the distant spikes. The \eqn{k^{th}} element of the array is the estimated correlation between \eqn{k^{th}} sample and population PC scores. If \code{n.spikes.out} is specified, only first \code{n.spikes.out} many correlations are returned.}
\item{shrinkage}{An array of estimated asymptotic shrinkage factors corresponding to the distant spikes. If \code{n.spikes.out} is specified, only first \code{n.spikes.out} many shrinkage factors are returned.}
\item{loss}{If \code{method="l.gsp"}, L-infinity loss function for the spectrum estimation is returned.}
\item{nonspikes}{If \code{nonspikes.out=TRUE}, estimated non-spikes are returned. If \eqn{\lambda}-estimation method is used then this is a numeric vector of length \code{p-n.spikes}. If OSP model based method is used then this is a scalar number.}
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
\code{\link{select.nspike},\link{pc_adjust}}
}
\examples{
data(hapmap)
#n = 198, p = 75435 for this data

####################################################
\dontrun{
train.eval<-hapmap$train.eval
n<-hapmap$nSamp
p<-hapmap$nSNP

m<-select.nspike(train.eval,p,n,n.spikes.max=10,evals.out=FALSE)$n.spikes
out<-hdpc_est(train.eval, p, n, method = "d.gsp", 
n.spikes=m, n.spikes.out=2, nonspikes.out = FALSE)	#Output 2 spikes, no non-spike

out<-hdpc_est(train.eval, p, n, method = "l.gsp", 
n.spikes=m, nonspikes.out = FALSE)	#Output m many spikes, no non-spike

out<-hdpc_est(train.eval, p, n, method = "l.gsp", 
n.spikes.max=10, nonspikes.out = TRUE)	#Output all eigenvalues

out<-hdpc_est(train.eval, p, n, method = "osp", 
n.spikes=m, n.spikes.out=2, nonspikes.out = TRUE)	#Output m many spikes, no non-spike
}
}
\keyword{multivariate}
\keyword{optimize}