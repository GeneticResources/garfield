\name{garfield.plot}
\alias{garfield.plot}
\title{
    Garfield plotting function
}
\description{
    \code{garfield.plot} is used for visualization of the enrichment
        analysis results obtained by the \code{garfield.run} permutation 
        step. Internally, it uses \code{garfield.plot.fnx} which has been
        adapted from the 'radial.plot' function from the 'plotrix' package.
}
\usage{
    garfield.plot(input_file, num_perm = 100000, output_prefix = "plot", 
        plot_title = "", filter = 10, tr = Inf)
}
\arguments{
    \item{input_file}{
        Input file name as produced by garfield.run permutation step.
    }
    \item{num_perm}{
        Number of permutations used in the garfield.run permutation step.
    }
    \item{output_prefix}{
        Figure file prefix. This would create the following files 
            output_prefix.Chromatin_States.pdf, output_prefix.Footprints.pdf, 
            output_prefix.Histone_Modifications.pdf, output_prefix.Peaks.pdf, 
            output_prefix.FAIRE.pdf, output_prefix.Genic.pdf, 
            output_prefix.Hotspots.pdf and output_prefix.TFBS.pdf
    }
    \item{plot_title}{
        Optional figure title
    }
    \item{filter}{
        Optional filter for the minimum number of variants at a given
        threshold. Minimum of 1 should be used, but advisable to set to 
        a larger value (e.g. 10).
    }
    \item{tr}{
        Threshold for denoting significance of an observed enrichment
        on a -log10 scale. A value of Inf denotes using default threshold 
        of -log10(0.05/498).
    }
}
\value{
    No value is produced, instead output files are generated. 
    See Details and 'output_prefix' for more information.
}
\details{
    This function is used for visualization of the enrichment 
    analysis results and produces pdf figures for each class of 
    annotations. Each figure shows the FE values (radial values) at 
    different GWAS thresholds (bottom legend) for each annotation 
    (outer circle and right legend). It further shows the significance 
    at the top 4 GWAS thresholds (if present) as dots on the outer 
    circle, with the most stringent threshold being shown at the 
    inner most side.
}
\references{
    Valentina Iotchkova, Graham Ritchie, Matthias Geihs, Sandro Morganella, 
    Josine Min, Klaudia Walter, Nicholas Timpson, UK10K Consortium, Ian Dunham, 
    Ewan Birney and Nicole Soranzo. GARFIELD - GWAS Analysis of Regulatory or 
    Functional Information Enrichment with LD correction. In preparation
}
\author{
    Sandro Morganella <email: sm22@sanger.ac.uk>

    Maintainer: Valentina Iotchkova <email: vi1@sanger.ac.uk>
}
\seealso{
    \code{\link{garfield.run}}, \code{\link{garfield}},
    \code{\link{garfield.plot.fnx}}
}
\examples{

garfield.run("tmp", data.dir=system.file("extdata",package = "garfield"),
    trait="trait",run.option = "prep", chrs = c(22), 
    exclude = c(895, 975, 976, 977, 978, 979, 98))

garfield.run("tmp", data.dir=system.file("extdata",package = "garfield"),
    trait="", run.option = "perm", nperm = 1000, 
    thresh = c(0.001, 1e-04, 1e-05), pt_thresh = c(1e-04, 1e-05),
    maf.bins = 2, tags.bins = 3, tss.bins = 3, prep.file = "tmp.prep",
    optim_mode = TRUE, minit = 100, thresh_perm = 0.05)

garfield.plot("tmp.perm", num_perm = 1000, output_prefix = "tmp",
    plot_title = "Sample run", filter = 1, tr = -log10(0.05))
#system("ls -lh tmp.*.pdf")

}
