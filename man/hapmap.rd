\name{hapmap}
\alias{hapmap}
\title{Example dataset - Hapmap Phase III}
\description{
The example dataset is from the Hapmap Phase III project (\url{https://www.ncbi.nlm.nih.gov/variation/news/NCBI_retiring_HapMap/}). Our training sample
consisted of unrelated individuals from two different populations: a) Utah residents
with Northern and Western European ancestry (CEU), and b) Toscans in Italy (TSI).
We present the eigenvalues and PC scores obtained from performing PCA on the SNPs on chromosome 7.
}
\format{
  This example dataset is a list containing the following elements:
  \describe{
    \item{train.eval}{Sample eigenvalues of the training sample.}
    \item{trainscore}{PC scores of the training sample. This has PC1 and PC2 scores
for 198 observations.}
    \item{testscore}{We obtained the predicted scores by leaving one observation 
out at a time, applying PCA to the rest of the data and then predicting the 
PC score of the left out observation. This has PC1 and PC2 scores of 198 observations.}
    \item{nSamp}{Number of observations in the training set = 198.}
    \item{nSNP}{Number of SNPs on chromosome 7.}
     }
}


