# Copyright (C) 2014 Genome Research Ltd / EMBL - 
# European Bioinformatics Institute 
#
# The function garfield.plot.fnx in this file is based on the 'radial.plot'
# function from the R package 'plotrix'. The original code is released under
# the GPL v2 |GPL v3 license, which has been modified here and all changes
# to the original code are released under the GPL v3 license.
#
# The original license information can be found at 
# http://cran.r-project.org/web/packages/plotrix/index.html

garfield.plot <- function(input_file, num_perm=100000, output_prefix="plot",
                          plot_title="", filter=10, tr=Inf){ 
    if (tr==Inf) tr=-log10(0.05/498)
    if (file.exists(input_file)){
        input <- read.table(input_file, header=TRUE)
    } else {
        print(paste("Error! Input file" , input_file, "does not exist!",
            sep=""))
        return(1);
    }
    input = input[which(input[,7]>=filter),]
    input[which(input[,3]==-1),3] = NA # set -1 FE values to NA
    input[which(input[,4]==0),4] = 1/num_perm ## use num_perm to transform
        ## 0 pvalues to a finite limit of the -log10 scale
    input[which(input[,4]==(-1)),4] = NA # set -1 pvalues to NA 
    input[which(is.na(input[,3])),3] = 0
    input$Category = as.character(input$Category)
    if (length(which(is.na(input$Category)))>1)
        input$Category[which(is.na(input$Category))] = "Custom"
    input$Category = as.factor(input$Category)

    for (category in levels(input$Category)){
        ids = which(input$Category==category)
        thresholds = sort(unique(input[ids,2]))
        thresholdsP = sort(unique(input[ids,2][which(!is.na(input[ids,4]))]))
        annotations = unique(as.character(input[ids,1]))
        if (category %in% c("Genic","Histone_Modifications",
            "Chromatin_States")){
            tissues = as.character(input$Type[ids][match(annotations,
                input[ids,1])])
            nms = as.character(input$Celltype[ids][match(annotations,
                input[ids,1])])
            compact=FALSE
            tissue_label = "Feature"
            if (category %in% c("Genic")) nms=tissues
        } else if (category %in% c("TFBS","FAIRE","Hotspots","Peaks",
            "Footprints")){
            nms = as.character(input$Celltype[ids][match(annotations,
                input[ids,1])])
            tissues = as.character(input$Tissue[ids][match(annotations,
                input[ids,1])])
            compact=TRUE
            tissue_label = "Tissue"
            if (category %in% c("Hotspots","Peaks","Footprints")){
                nms=tissues
            }
        } else {
            nms = as.character(input$Celltype[ids][match(annotations,
                input[ids,1])])
            tissues = as.character(input$Tissue[ids][match(annotations,
                input[ids,1])])
        }
        DATA = matrix(NA, nrow=length(thresholds)+1, ncol=length(annotations))
        DATA_p = matrix(NA, nrow=length(thresholdsP)+1, 
            ncol=length(annotations))
        for (j in 1:length(annotations)){
            for (i in 1:length(thresholds)){
                DATA[i,j] = input[which(input[,1]==annotations[j] & 
                    input[,2]==thresholds[i]),3]
            }
            for (i in 1:length(thresholdsP)){
                DATA_p[i,j] = input[which(input[,1]==annotations[j] & 
                    input[,2]==thresholdsP[i]),4]
            }
        }
        DATA[length(thresholds)+1,] = DATA_p[length(thresholdsP)+1,] = 1
        DATA_p = -log10(DATA_p)
        ann.col = colorRampPalette(c("tomato","skyblue3","yellow","brown2",
            "lightgreen","lightgoldenrod3","purple","pink","darkblue","gray",
            "darkgreen"))( length(unique(tissues)) )[as.numeric(
            as.factor(tissues))]
        ann.col.mx = matrix(ann.col,nrow=length(thresholdsP),
            ncol=length(ann.col),byrow=TRUE)
        for (i in 1:(length(thresholdsP))){
            ann.col.mx[i,which(DATA_p[i,]<tr)]=0
        }
        ord = order(tissues)
        col.thresh = colorRampPalette(c("black","firebrick3","tomato",
            "RosyBrown2","dodgerblue3","skyblue2","lightskyblue1","gray70",
            "blanchedalmond"))( length(unique(thresholds))+1 )
        if (length(thresholdsP)<4){
            rws = length(thresholdsP):1
        } else {
            rws = 4:1
        }
        tissues = gsub("_", " ", as.character(tissues))
        nms = gsub("_", " ", as.character(nms))
        pdf(paste(output_prefix,".",category,".pdf",sep=""),12,10)
        layout(matrix(c(1,2,3,3),nrow=2), widths = c(9,2),heights = c(9, 1), respect = FALSE)
        par(oma = c(0, 0, 0, 0))
        garfield.plot.fnx(DATA[,ord],ann.col.mx=matrix(ann.col.mx[rws,ord],
            nrow=length(rws)), ann.col=ann.col[ord], ann.pch=15, rp.type="p",
            line.col=col.thresh,show.grid=TRUE, show.radial.grid=TRUE,
            labels=nms[ord],breaks=tissues[ord], radlab=TRUE,cex.axis=0.1,
            cex.lab=0.1, mar = c(5.5, 2, 6, 1), label.prop=1.05,
            poly.col=col.thresh, compact=compact)
        if (plot_title!=""){
            title(main=paste(plot_title," ",category,sep=""),line=7,cex=2) 
        }
        par(mar = c(0, 0, 0, 0))
        plot(1:2,1:2, type="n", axes=FALSE, xlab="", ylab="")
        legend("bottom", c(thresholds,"1"),col=col.thresh,lty=1,lwd=6,
            title="GWAS P-value Threshold",horiz=TRUE,cex=1,bty="n")
        plot(1:2,1:2, type="n", axes=FALSE, xlab="", ylab="")
        legend("right",unique(tissues[ord]),col=unique(ann.col[ord]),
            lty=1,lwd=5,title=tissue_label,cex=1,bty="n")
        dev.off()
    }
}


garfield.plot.fnx<-function (lengths, radial.pos = NULL, labels = NA, 
    breaks= NA, label.pos = NULL, radlab = FALSE, start = 0, 
    clockwise = FALSE, rp.type = "r", label.prop = 1.05, main = "", xlab = "", 
    ylab = "", line.col = par("fg"), lty = par("lty"), lwd = par("lwd"), 
    mar = c(2, 2, 3, 2), show.grid = TRUE, show.grid.labels = 4, 
    show.radial.grid = TRUE, grid.col = "grey", grid.bg = "transparent", 
    grid.left = FALSE, grid.unit = NULL, point.symbols = 1, 
    point.col = par("fg"), show.centroid = FALSE, radial.lim = NULL, 
    radial.labels = NULL, poly.col = NA, add = FALSE, ann.col=1,ann.pch=15,
    ann.col.mx=1,compact=TRUE,...){
    
    lengths.ann = rep(max(lengths,na.rm=TRUE),ncol(lengths))
    radial.lim=c(0,max(lengths,na.rm=TRUE))
    
    length.dim <- dim(lengths)
    if (is.null(length.dim)) {
        npoints <- length(lengths)
        nsets <- 1
        lengths <- matrix(lengths, nrow = 1)
    }
    else {
        npoints <- length.dim[2]
        nsets <- length.dim[1]
        lengths <- as.matrix(lengths)
    }
    lengths <- lengths - radial.lim[1]
    lengths[lengths < 0] <- NA
    if (is.null(radial.pos[1])) 
        radial.pos <- seq(0, pi * (2 - 2 * (rp.type != "l")/npoints), 
            length.out = npoints)
    radial.pos.dim <- dim(radial.pos)
    if (is.null(radial.pos.dim)) 
        radial.pos <- matrix(rep(radial.pos, nsets), nrow = nsets, 
            byrow = TRUE)
    else radial.pos <- as.matrix(radial.pos)
    if (rp.type == "l") {
        clockwise <- TRUE
        start <- pi/2
    }
    if (clockwise) 
        radial.pos <- -radial.pos
    if (start) {
        radial.pos <- radial.pos + start
    }
    if (show.grid) {
        if (length(radial.lim) < 3) 
            grid.pos <- pretty(radial.lim)
        else grid.pos <- radial.lim
        if (grid.pos[1] < radial.lim[1]) 
            grid.pos <- grid.pos[-1]
        maxlength <- max(grid.pos - radial.lim[1])
        angles <- seq(0, 1.96 * pi, by = 0.04 * pi)
    }
    else {
        grid.pos <- NA
        maxlength <- diff(radial.lim)
    }
    oldpar <- par("xpd", "mar", "pty")
    if (!add) {
        par(mar = mar, pty = "s")
        plot(c(-maxlength, maxlength), c(-maxlength, maxlength), 
            type = "n", axes = FALSE, main = main, xlab = xlab, 
            ylab = ylab)
        if (show.grid) {
            for (i in seq(length(grid.pos), 1, by = -1)) {
                xpos <- cos(angles) * (grid.pos[i] - radial.lim[1])
                ypos <- sin(angles) * (grid.pos[i] - radial.lim[1])
                polygon(xpos, ypos, border = grid.col, col = grid.bg)
            }
        }
    }
    par(xpd = TRUE)
    if (length(line.col) < nsets) 
        line.col <- 1:nsets
    if (length(rp.type) < nsets) 
        rp.type <- rep(rp.type, length.out = nsets)
    if (length(point.symbols) < nsets) 
        point.symbols <- rep(point.symbols, length.out = nsets)
    if (length(point.col) < nsets) 
        point.col <- rep(point.col, length.out = nsets)
    if (length(poly.col) < nsets) 
        poly.col <- rep(poly.col, length.out = nsets)
    if (length(lty) < nsets) 
        lty <- rep(lty, length.out = nsets)
    if (length(lwd) < nsets) 
        lwd <- rep(lwd, length.out = nsets)

    for (i in 1:nsets) {
        if (nsets > 1) {
            linecol <- line.col[i]
            polycol <- poly.col[i]
            pointcol <- point.col[i]
            pointsymbols <- point.symbols[i]
            ltype <- lty[i]
            lwidth <- lwd[i]
        }
        else {
            linecol <- line.col
            polycol <- poly.col
            pointcol <- point.col
            pointsymbols <- point.symbols
            ltype <- lty
            lwidth <- lwd
        }
        rptype <- unlist(strsplit(rp.type[i], ""))
        if (match("s", rptype, 0)) {
            if (is.null(pointsymbols)) 
                pointsymbols <- i
            if (is.null(pointcol)) 
                pointcol <- i
        }
        xpos <- cos(radial.pos[i, ]) * lengths[i, ]
        ypos <- sin(radial.pos[i, ]) * lengths[i, ]
        if (match("r", rptype, 0)) 
            segments(0, 0, xpos, ypos, col = linecol, lty = ltype, 
                lwd = lwidth, ...)
        if (match("p", rptype, 0)) 
            polygon(xpos, ypos, border = linecol, col = polycol, 
                lty = ltype, lwd = lwidth, ...)
        if (match("s", rptype, 0)) 
            points(xpos, ypos, pch = pointsymbols, col = pointcol, 
                ...)
        if (match("l", rptype, 0)) 
            lines(xpos, ypos, lty = ltype, lwd = lwidth, col = linecol, 
                ...)
        if (show.centroid) 
            if (match("p", rptype, 0)) {
                nvertices <- length(xpos)
                polygonarea <- xpos[nvertices] * ypos[1] - xpos[1] * 
                  ypos[nvertices]
                for (vertex in 1:(nvertices - 1)) polygonarea <- polygonarea + 
                  xpos[vertex] * ypos[vertex + 1] - xpos[vertex + 
                  1] * ypos[vertex]
                polygonarea <- polygonarea/2
                centroidx <- (xpos[nvertices] + xpos[1]) * (xpos[nvertices] * 
                  ypos[1] - xpos[1] * ypos[nvertices])
                centroidy <- (ypos[nvertices] + ypos[1]) * (xpos[nvertices] * 
                  ypos[1] - xpos[1] * ypos[nvertices])
                for (vertex in 1:(nvertices - 1)) {
                  centroidx <- centroidx + (xpos[vertex] + xpos[vertex + 
                    1]) * (xpos[vertex] * ypos[vertex + 1] - 
                    xpos[vertex + 1] * ypos[vertex])
                  centroidy <- centroidy + (ypos[vertex] + ypos[vertex + 
                    1]) * (xpos[vertex] * ypos[vertex + 1] - 
                    xpos[vertex + 1] * ypos[vertex])
                }
                points(centroidx/(6 * polygonarea), centroidy/(6 * 
                  polygonarea), col = point.col[i], pch = point.symbols[i], 
                  cex = 2, ...)
            }
            else points(mean(xpos), mean(ypos), col = pointcol, 
                pch = pointsymbols, cex = 2, ...)

    }
    xpos.ann <- cos(radial.pos[1, ]) * maxlength * 1.02
    ypos.ann <- sin(radial.pos[1, ]) * maxlength * 1.02

    points(xpos.ann, ypos.ann, pch = ann.pch, col = ann.col,...)
    for (ij in 1:nrow(ann.col.mx)){
        xpos.ann2 <- cos(radial.pos[1, ]) * (maxlength * (0.992 - 0.015 * 
            (ij - 1)))
        ypos.ann2 <- sin(radial.pos[1, ]) * (maxlength * (0.992 - 0.015 * 
            (ij - 1)))
        points(xpos.ann2, ypos.ann2, pch = 19, col = ann.col.mx[ij,], 
            cex=0.5,...) 
    }
    if (!add) {
        if (is.na(labels[1])) {
            label.pos <- seq(0, 1.8 * pi, length = 9)
            labels <- as.character(round(label.pos, 2))
            lablen <- length(labels)
            breaks.pos <- seq(0, pi * (2 - 2/length(breaks)), 
                length.out = length(breaks))
        }
        if (is.null(label.pos[1])) {
            lablen <- length(labels)
            label.pos <- seq(0, pi * (2 - 2/lablen), length.out = lablen)
            breaks.pos <- seq(0, pi * (2 - 2/length(breaks)), 
                length.out = length(breaks))
        }
        if (clockwise) 
            label.pos <- -label.pos
        if (start) 
            label.pos <- label.pos + start
        xpos <- cos(breaks.pos-pi/lablen) * maxlength
        ypos <- sin(breaks.pos-pi/lablen) * maxlength
            if (show.radial.grid) 
                segments(0, 0, xpos[which(!duplicated(breaks))], 
                    ypos[which(!duplicated(breaks))],col = "gray",lty=4)        
            xpos <- cos(label.pos) * maxlength * label.prop
            ypos <- sin(label.pos) * maxlength * label.prop
            label.adj <- round(abs(1 - cos(label.pos))/2-10^{-10})
            if (radlab) {
                labn = names(table(labels))
                nlab = as.numeric(table(labels))/sum(as.numeric(table(labels
                    )))*length(as.numeric(table(labels)))
                labelsn = nlab[match(labels,labn)]
                for (label in 1:length(labels)) {
                    if (!is.na(labels[label])){
                        flr=ceiling(median(which(labels == labels[label])))
                    } else {
                        flr=ceiling(median(which(is.na(labels))))
                    }
                    if ((flr==label & compact==TRUE) | compact==FALSE){
                        labelsrt <- (180 * label.pos[label]/pi) + 180 * 
                            (label.pos[label] > pi/2 && label.pos[label] < 
                            3 * pi/2)
                        text(xpos[label], ypos[label], labels[label],
                            cex = (0.45+0.3*labelsn[label]), srt = labelsrt, 
                            adj = label.adj[label],
                            #adj = label.adj[label]+lengths.ann[1]*1.07/1500,
                            col=1)
                        #text(xpos[label], ypos[label], labels[label],
                        #    cex = (0.35+0.3*labelsn[label]), srt = labelsrt, 
                        #    adj = label.adj[label],col=ann.col[label])
                    }
               }
          } else {
              for (label in 1:length(labels)) {
                  text(xpos[label], ypos[label], labels[label], 
                      cex = par("cex.axis"), adj = label.adj[label])
              }
          }
          if (show.grid.labels) {
              if (show.grid.labels%%2) {
                  ypos <- grid.pos - radial.lim[1]
                  xpos <- rep(0, length(grid.pos))
                  if (show.grid.labels == 1) 
                      ypos <- -ypos
                  } else {
                      xpos <- grid.pos - radial.lim[1]
                      ypos <- rep(0, length(grid.pos))
                      if (show.grid.labels == 2) 
                         xpos <- -xpos
                  }
                  if (is.null(radial.labels)) 
                      radial.labels = as.character(grid.pos)
                  if (!is.null(grid.unit)) 
                      radial.labels[length(grid.pos)] <- paste(radial.labels
                          [length(grid.pos)],grid.unit)
                  text(xpos+0.003*max(xpos), ypos+0.003*max(ypos), 
                  radial.labels, cex = par("cex.lab"),col="white")
                  text(xpos, ypos, radial.labels, cex = par("cex.lab"))
         }
    }
    invisible(oldpar)
}

