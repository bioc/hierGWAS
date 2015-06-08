\docType{package}
\name{hierGWAS}
\alias{hierGWAS-package}
\alias{hierGWAS}
\title{Asessing statistical significance in predictive GWA studies}

\description{
Testing individual SNPs, as well as arbitrarily large groups of SNPs in GWA
studies, using a joint model of all SNPs. The method controls the FWER, and
provides an automatic, data-driven refinement of the SNP clusters to smaller
groups or single markers.
}
\details{
hierGWAS is a package designed to assess statistical significance in GWA studies,
using a hierarchical approach.

There are 4 functions provided: \code{\link[hierGWAS]{cluster.snp}},
\code{\link[hierGWAS]{multisplit}}, \code{\link[hierGWAS]{test.hierarchy}} and
\code{\link[hierGWAS]{compute.r2}}. \code{\link[hierGWAS]{cluster.snp}} performs
the hierarchical clustering of the SNPs, while \code{\link[hierGWAS]{multisplit}}
runs multiple penalized regressions on repeated random subsamples. These 2 functions
need to be executed before \code{\link[hierGWAS]{test.hierarchy}}, which does the
hierarchical testing, though the order in which the 2 functions are executed does
not matter. \code{\link[hierGWAS]{test.hierarchy}} provides the final output of
the method: a list of SNP groups or individual SNPs, along with their corresponding
p-values. Finally, \code{\link[hierGWAS]{compute.r2}} computes the explained variance
of an arbitrary group of SNPs, of any size. This group can encompass all SNPs, SNPs
belonging to a certain chromosome, or an individual SNP.
}
\author{
Laura Buzdugan laura.buzdugan@stat.math.ethz.ch
}
\references{
Buzdugan, L. et al. (2015), Assessing statistical significance in
predictive genome-wide association studies (unpublished)
}

