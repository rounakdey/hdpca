\name{pc_adjust}
\alias{pc_adjust}
\title{
Adjusting shrinkage in PC scores
}
\description{
Adjusts the shrinkage bias in the predicted PC scores based on the estimated shrinkage factors.
}
\usage{
pc_adjust(train.eval, p, n, test.scores, method = c("d.gsp", "l.gsp", "osp"),
n.spikes, n.spikes.max, smooth = TRUE)
}
\arguments{
  \item{train.eval}{
Numeric vector containing the sample eigenvalues. The vector must have dimension \code{n} or \code{n-1}, it may be unordered.
}
  \item{p}{
The number of features.
}
  \item{n}{
The number of training samples.
}
  \item{test.scores}{
An \eqn{m\times k} matrix or data frame containing the first \eqn{k} predicted PC scores of \eqn{m} many test samples.
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
The \eqn{(i,j)^{th}} element of \code{test.scores} should denote the \eqn{j^{th}} predicted PC score for the \eqn{i^{th}} subject in the test sample.

At least one of \code{n.spikes} and \code{n.spikes.max} must be provided. If \code{n.spikes} is provided then \code{n.spikes.max} is ignored, else \code{n.spikes.max} is used to find out the number of distant spikes using \code{\link{select.nspike}}.

The argument \code{nonspikes.out} is ignored if \code{method="d.gsp"} or "\code{osp}".

The argument \code{smooth} is useful when the user assumes the population spectral distribution to be continuous.
}
\value{
A matrix containing the bias-adjusted PC scores. The dimension of the matrix is the same as the dimension of \code{test.scores}.

A printed message shows the number of top PCs that were adjusted for shrinkage bias.
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
\code{\link{hdpc_est},\link{select.nspike}}
}
\examples{
data(hapmap)
#n = 198, p = 75435 for this data

####################################################
\dontrun{
#First estimate the number of spikes and then adjust test scores based on that
train.eval<-hapmap$train.eval
n<-hapmap$nSamp
p<-hapmap$nSNP
trainscore<-hapmap$trainscore
testscore<-hapmap$testscore

m<-select.nspike(train.eval,p,n,n.spikes.max=10,evals.out=FALSE)$n.spikes
score.adj.o1<-pc_adjust(train.eval,p,n,testscore,method="osp",n.spikes=m)
score.adj.d1<-pc_adjust(train.eval,p,n,testscore,method="d.gsp",n.spikes=m)
score.adj.l1<-pc_adjust(train.eval,p,n,testscore,method="l.gsp",n.spikes=m)

#Or you can provide an upper bound n.spikes.max
score.adj.o2<-pc_adjust(train.eval,p,n,testscore,method="osp",n.spikes.max=10)
score.adj.d2<-pc_adjust(train.eval,p,n,testscore,method="d.gsp",n.spikes.max=10)
score.adj.l2<-pc_adjust(train.eval,p,n,testscore,method="l.gsp",n.spikes.max=10)

#Plot the training score, test score, and adjusted scores
plot(trainscore,pch=19)
points(testscore,col='blue',pch=19)
points(score.adj.o1,col='red',pch=19)
points(score.adj.d2,col='green',pch=19)
}
}
\keyword{multivariate}