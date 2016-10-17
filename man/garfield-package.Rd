\name{garfield-package}
\alias{garfield-package}
\alias{garfield}
\docType{package}
\title{
GARFIELD - GWAS Analysis of Regulatory or Functional Information Enrichment
with LD correction
}
\description{
GARFIELD leverages GWAS findings with regulatory or functional annotations
to find features relevant to a phenotype of interest. It performs greedy
pruning of GWAS SNPs (LD r2>0.1) and then annotates them based on functional
information overlap. Next, it quantifies Fold Enrichment (FE) at various
GWAS signficance cut-offs and assesses them by permutation testing, while
matching for minor allele frequency, distance to nearest transcription start
site and number of LD proxies (r2>0.8). Finally, it includes visualization
tools.
}
\details{
\tabular{ll}{
Package: \tab garfield\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-12-14\cr
License: \tab GPL-3\cr
}
See garfield.run for example analysis usage and garfield.plot for plotting
examples.
}
\author{
Sandro Morganella <email: sm22@sanger.ac.uk>

Maintainer: Valentina Iotchkova <email: vi1@sanger.ac.uk>
}
\references{
Valentina Iotchkova, Graham Ritchie, Matthias Geihs, Sandro Morganella,
Josine Min, Klaudia Walter, Nicholas Timpson, UK10K Consortium, Ian Dunham,
Ewan Birney and Nicole Soranzo. GARFIELD - GWAS Analysis of Regulatory or
Functional Information Enrichment with LD correction. In preparation
}
\keyword{ package }
\seealso{
\code{garfield.run}, \code{garfield.plot}
}

