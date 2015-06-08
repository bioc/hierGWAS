\name{simGWAS}
\docType{data}
\alias{simGWAS}
\title{Simulated GWAS data}
\description{
  This data set was simulated using \code{PLINK}. Please refer to the vignette
  for more details.
}
\format{
  The dataset contains the following components:
  \describe{
    \item{\code{SNP.1}}{The first SNP, of dimension \code{500 x 1}. Each
    row represents a subject.}
    \item{}{...}
    \item{\code{SNP.1000}}{The last SNP, of dimension \code{500 x 1}. Each
    row represents a subject.}
    \item{\code{y}}{The response vector. It can be continuous or discrete.}
    \item{\code{sex}}{The first covariate, represeting the sex of the subjects:
    0 for men and 1 for women.}
    \item{\code{age}}{The second covariate, represeting the age of the subjects.}
  }
}
\value{data.frame}
\usage{simGWAS}
\examples{
data(simGWAS)
}
\keyword{datasets}
