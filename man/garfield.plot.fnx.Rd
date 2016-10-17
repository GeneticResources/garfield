\name{garfield.plot.fnx}
\alias{garfield.plot.fnx}
\title{
    Internal radial plotting function for the garfield.plot function
}
\description{
    This function has been adapted from the 'radial.plot' function from 
        the 'plotrix' package.
}
\usage{
    garfield.plot.fnx(lengths, radial.pos = NULL, labels = NA, breaks = NA,
        label.pos = NULL, radlab = FALSE, start = 0, clockwise = FALSE,
        rp.type = "r", label.prop = 1.05, main = "", xlab = "", ylab = "",
        line.col = par("fg"), lty = par("lty"), lwd = par("lwd"),
        mar = c(2, 2, 3, 2), show.grid = TRUE, show.grid.labels = 4,
        show.radial.grid = TRUE, grid.col = "grey", grid.bg = "transparent",
        grid.left = FALSE, grid.unit = NULL, point.symbols = 1,
        point.col = par("fg"), show.centroid = FALSE, radial.lim = NULL,
        radial.labels = NULL, poly.col = NA, add = FALSE, ann.col = 1,
        ann.pch = 15, ann.col.mx = 1, compact = TRUE, ...)
}
\arguments{
    \item{lengths}{
        A numeric data vector or matrix. If 'lengths' is a matrix,
            the rows will be considered separate data vectors.
    }
    \item{radial.pos}{
        A numeric vector or matrix of positions in radians.  These
            are interpreted as beginning at the right (0 radians) and
            moving counterclockwise. If 'radial.pos' is a matrix, the
            rows must correspond to rows of 'lengths'.
    }
    \item{labels}{
        Character strings to be placed at the outer ends of the
            lines. If set to NA, will suppress printing of labels, but if
            missing, the radial positions will be used.
    }
    \item{breaks}{
        A vector of (potentially different) labels to 'labels' according to
            which to draw radial lines.
    }
    \item{label.pos}{
        The positions of the labels around the plot in radians.
    }
    \item{radlab}{
        Whether to rotate the outer labels to a radial orientation.
    }
    \item{start}{
        Where to place the starting (zero) point. Defaults to the 3
            o'clock position.
    }
    \item{clockwise}{
        Whether to interpret positive positions as clockwise from
            the starting point. The default is counterclockwise.
    }
    \item{rp.type}{
        Whether to draw (r)adial lines, a (p)olygon, (s)ymbols or
           some combination of these. If 'lengths' is a matrix and
           rp.type is a vector, each row of 'lengths' can be displayed
           differently.
    }
    \item{label.prop}{
        The label position radius as a proportion of the maximum
            line length.
    }
    \item{main}{
        The title for the plot.
    }
    \item{xlab,ylab}{
        Normally x and y axis labels are suppressed.
    }
    \item{line.col}{
        The color of the radial lines or polygons drawn.
    }
    \item{lty}{
        The line type(s) to be used for polygons or radial lines.
    }
    \item{lwd}{
        The line width(s) to be used for polygons or radial lines.
    }
    \item{mar}{
        Margins for the plot. Allows the user to leave space for
            legends, long labels, etc.
    }
    \item{show.grid}{
        Logical - whether to draw a circular grid.
    }
    \item{show.grid.labels}{
        Whether and where to display labels for the grid -
            see Details.
    }
    \item{show.radial.grid}{
        Whether to draw radial lines to the plot labels.
    }
    \item{grid.col}{
        Color of the circular grid.
    }
    \item{grid.bg}{
        Fill color of above.
    }
    \item{grid.left}{
        Whether to place the radial grid labels on the left side.
    }
    \item{grid.unit}{
        Optional unit description for the grid.
    }
    \item{point.symbols}{
        The symbols for plotting (as in pch).
    }
    \item{point.col}{
        Colors for the symbols.
    }
    \item{show.centroid}{
        Whether to display a centroid.
    }
    \item{radial.lim}{
        The range of the grid circle. Defaults to
            'pretty(range(lengths))', but if more than two values are
            passed, the exact values will be displayed.
    }
    \item{radial.labels}{
        Optional labels for the radial grid. The default is the
            values of radial.lim.
    }
    \item{poly.col}{
        Fill color if polygons are drawn. Use NA for no fill.
    }
    \item{add}{
        Whether to add one or more series to an existing plot.
    }
    \item{ann.col}{
        A vector of colours for symbols to be drawn on the outer circle.
    }
    \item{ann.pch}{
        A pch value for symbols to be drawn on the outer circle.
    }
    \item{ann.col.mx}{
        A matrix of colours for dots to be drawn on the outer circle, 
            inside the symbols with 'ann.col' colours.
    }
    \item{compact}{
       A flag to specify if labels should be presented in a compact form or not.
    }
    \item{\dots}{
        Additional arguments are passed to 'plot'.
    }
}
\value{
    No value, but a plot is produced. See Details for further information.
}
\details{
    This function is only used internally as part of the 
        \code{\link{garfield.plot}} function.
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
    \code{\link{garfield.plot}}, \code{\link{garfield.run}},
    \code{\link{garfield}}
}
\examples{

DATA = rbind(rnorm(10,5,0.5),rnorm(10,3.5,0.5),rnorm(10,2,0.5))
garfield.plot.fnx(DATA,ann.col.mx=DATA!=0, ann.col=rep(1:2,each=5),
    ann.pch=15, rp.type="p",line.col=1:3,show.grid=TRUE, show.radial.grid=TRUE,
    labels=paste("label ",c(1:10)," ",sep=""),breaks=(1:10), radlab=TRUE, 
    poly.col=1:3)
}

