\name{garfield.run}
\alias{garfield.run}
\alias{garfield_prep}
\alias{garfield_perm}
\title{
    GARFIELD enrichment analysis function
}
\description{
    \code{garfield.run} is used to perform greedy pruning of variants
    from a genome-wide association study, calculate fold enrichment and
    test its significance at a given genome-wide significance threshold.
}
\usage{
garfield.run(out.file, data.dir, trait, run.option = "complete",
    chrs = c(1:22, "X"), exclude = c(895, 975, 976, 977, 978, 979, 98),
    nperm = 100000, thresh = c(0.1, 0.01, 0.001, 1e-04, 1e-05, 1e-06,
    1e-07, 1e-08), pt_thresh = c(1e-05, 1e-06, 1e-07, 1e-08),
    maf.bins = 5, tags.bins = 5, tss.bins = 5, prep.file = "",
    optim_mode = 1, minit = 100, thresh_perm = 1e-04)
}
\arguments{
    \item{out.file}{
        Prefix for output file. Full garfield analysis creates out.file.prep
        and out.file.perm files from prunning and annotation step and
        permutation for significance testing step, respectively. These
        steps can additionally be run separately - for more information
        see run.option flag.
    }
    \item{data.dir}{
        Path to annotation and p-value files. The directory must contain
        "annotation","maftssd","pval" and "tags" subdirectories with per 
        chromosome files of input data. See \code{details} section for 
        further information.
    }
    \item{trait}{
        GWAS phenotype name. This must match a folder name in the
        data.dir folder. See details for more information.
    }
    \item{run.option}{
        an object specifying which part of the analysis to run. Valid 
        options are \code{complete}, \code{prep} and \code{perm}, where
        \code{prep} denotes the preparation step (pruning and annotation
        of variants), \code{perm} denotes the permutation step 
        (calculating fold enrichment and its significance) and 
        \code{complete} executes both the preparation and permutation steps.
    }
    \item{chrs}{
        A vector of the chromosomes for which to run the enrichment analysis.
        \code{chrs} can have all or subsets of values from \code{c(1:22,'X')}.
    }
    \item{exclude}{
        A numeric vector of indeces of annotations for which LD tags 
        should not be used for annotation of pruned variants. Value of -1 
        denotes using LD for annotation for all features.
    }
    \item{nperm}{
        A numeric value of the number of permutations to be performed.
    }
    \item{thresh}{
        A numeric vector of genome-wide significance thresholds to be used
        for fold enrichment calculation.
    }
    \item{pt_thresh}{
        A numeric vector of genome-wide significance thresholds to be used 
        for calculating the significance of the observed fold enrichment.
        All values must be contained in the \code{thresh} vector.
    }
    \item{maf.bins}{
        A numeric value denoting the number of bins for the minor allele 
        frequency matching during permutation testing. Must be greater 
        or equal to 1.
    }
    \item{tags.bins}{
        A numeric value denoting the number of bins for the number of LD
        tags (r2>=0.8) matching during permutation testing. Must be 
        greater or equal to 1.
    }
    \item{tss.bins}{
        A numeric value denoting the number of bins for the distance 
        to nearest transcription start site matching during permutation
        testing. Must be greater or equal to 1.
    }
    \item{prep.file}{
        File from pruning stage of algorithm. Only required if using 
        the 'run.option'=perm flag.
    }
    \item{optim_mode}{
        A binary flag denoting whether to run fast version of method 
        (1) or general version (0), where the fast version checks if 
        significance of a given enrichment would still be possible to 
        be obtained after 'minit' number of iterations and terminated 
        permutations if not. 
    }
    \item{minit}{
        An integer value for the minimum number of permutations to 
        be performed before checking if \code{thresh_perm} condition 
        can still be met. Only used if 'optim_mode'=1.
    }
    \item{thresh_perm}{
        After 'minit' number of permutations, at each iteration check
        if EmpPval can still reach 'thresh_perm' value. If not 
        terminate permutations and output obtained results at that 
        stage. Only used if 'optim_mode'=1.
    }
}
\value{
    No value is produced, instead output files are generated. 
    See Details and 'out.file' for more information.
}
\details{
    Output files:
    out.file.prep contains the genomic positions of pruned variants,
    p-values for association with the trait of interest, number of
    LD tags (r2>0.8), MAF, distance to the nearest TSS and binary 
    representation of annotation information (with LD-tagging r2>0.8).
    out.file.perm contains enrichment analysis results for each 
    annotation, where PThresh is the GWAS p-value threshold used 
    for analysis, FE denotes the fold enrichment statistic (equals -1
    if no sufficient data was available for the FE calculation),
    EmpPval shows the empirical p-value of enrichment (equals -1 if FE
    is calculated but significance of enrichment analysis is not run 
    at that threshold), NAnnotThresh - the number of variants at the
    threshold which are annotated with the given feature, NAnnot -
    the total number of annotated variants, NThresh - the total 
    number of variants at that threshold and N - the total number
    of pruned variants. The remaining columns show additional 
    information on the annotations used for analysis. 

    Data directory: data.dir, should point to a location 
    containing "annotation","maftssd","pval" and "tags" 
    subdirectories, where (i) the "pval" folder contains subfolders
    with trait names, which in turn contain per chromosome space
    separated files with genomic position in the first column and
    p-value in the second. They should be named chr1, chr2, etc, 
    and be numerically sorted with respect to genomic position; 
    (ii) "annotation" should contain per chromosome space separated 
    files with position in the first column and annotations in 
    stacked binary format in the second. The files should be 
    named chr1, chr2, etc and be sorted numerically according to 
    position. Additionally the directory should contain a 
    link_file.txt file that links the annotations to relevant 
    information about them; (iii) "maftssd" should contain per 
    chromosome space separated files with position in the first 
    column, minor allele frequency in the second and distance to 
    the nearest TSS in the third. The files should be named chr1,
    chr2, etc and be sorted numerically according to position; 
    (iv) "tags" should contain two subfolders named r01 and r08,
    which in turn contain per chromosome space separated files 
    with variant position in the first column and comma separated
    positions of all variants with r2>0.1 or 0.8, respectively 
    with the variant in the first column. The files should again 
    be named chr1, chr2, etc and be sorted numerically according 
    to position. For further information on the data.dir structure,
    please see 
    http://www.ebi.ac.uk/birney-srv/GARFIELD/documentation/GARFIELD.pdf

    Pre-computed LD (European samples - UK10K sequence data), MAF,
    TSS distance, p-value files for two example traits (Crohn's Disease 
    from the IBD Consortium and Height from the GIANT consortium) and 
    annotation files for 1005 GENCODE, ENCODE and Roadmap Epigenomics 
    annotations can be downloaded from 
    http://www.ebi.ac.uk/birney-srv/GARFIELD/package/garfield-data.tar.gz. 
    Note the data is 5.9Gb in compressed format and needs to be uncompressed 
    prior to analysis (83Gb). Variant genomic position (build 37) is 
    used as an identifier in all data files.

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
    \code{\link{garfield.plot}}, \code{\link{garfield}}
}
\examples{

garfield.run("tmp", data.dir=system.file("extdata",package = "garfield"),
    trait="trait",run.option = "prep", chrs = c(22),
    exclude = c(895, 975, 976, 977, 978, 979, 98))

garfield.run("tmp", data.dir=system.file("extdata",package = "garfield"),
    run.option = "perm", nperm = 1000, thresh = c(0.001, 1e-04, 1e-05),
    pt_thresh = c(1e-04, 1e-05), maf.bins = 2, tags.bins = 3, tss.bins = 3,
    prep.file = "tmp.prep", optim_mode = TRUE, minit = 100, thresh_perm = 0.05)

if (file.exists("tmp.perm")){ 
    perm = read.table("tmp.perm", header=TRUE)
    head(perm)
} else { print("Error: tmp.perm does not exist!") }

###### To get the sample data for enrichment analysis in European samples
###### execute the following - note this can take a long time to run and 
###### needs a substantial disk space (see Details)
#
# download data and decompress
# system("wget http://www.ebi.ac.uk/birney-srv/GARFIELD/package/
# garfield-data.tar.gz")
# system("tar -zxvf garfield-data.tar.gz")
#
# if downloaded in current working directory use the following to execute 
# garfield, otherwise please change data.dir location
# garfield.run("cd-meta.output", data.dir="garfield-data", trait="cd-meta",
# run.option = "prep", chrs = c(1:22), exclude = c(895, 975, 976, 977, 978,
# 979, 980))
#
# garfield.run("cd-meta.output", data.dir="garfield-data", run.option = "perm",
# nperm = 100000, thresh = c(0.1,0.01,0.001, 1e-04, 1e-05, 1e-06, 1e-07, 1e-08),
# pt_thresh = c(1e-05, 1e-06, 1e-07, 1e-08), maf.bins = 5, tags.bins = 5,
# tss.bins = 5, prep.file = "cd-meta.output.prep", optim_mode = TRUE,
# minit = 100, thresh_perm = 0.0001)
#
# garfield.plot("cd-meta.output.perm", num_perm = 100000,
# output_prefix = "cd-meta.output", plot_title = "Crohn's Disease",
# filter = 10, tr = -log10(0.05/498))

}
