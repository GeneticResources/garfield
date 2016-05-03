# Copyright (C) 2014 Genome Research Ltd / EMBL - 
# European Bioinformatics Institute 

garfield.run <- function(out.file, data.dir, trait, run.option="complete", 
    chrs=c(1:22,"X"), exclude=c(895,975,976,977,978,979,98), nperm=100000,
    thresh=c(1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8), 
    pt_thresh=c(1e-5,1e-6,1e-7,1e-8), maf.bins=5, tags.bins=5, tss.bins=5,

    prep.file="", optim_mode=1, minit=100, thresh_perm=0.0001){
        if( !( run.option %in% c("complete", "prep", "perm")) ){
            message("ERROR: Invalid run.option")
            message("Valid Options are: \"complete\", \"prep\", \"perm\"")
        }
        message("\n--- START COMPUTATION ---\n")
        if( run.option %in% c("complete","prep")){
            ## Run preparation script
            message("\t[RUN] PREPARING DATA")
            prep.file <- paste(out.file, ".prep",sep="")
            if(file.exists(prep.file)){
                file.remove(prep.file)
                file.create(prep.file)
            }
            for(i in 1:length(chrs)){
                message("\t\tAnalyizing Chromosome ",chrs[i])
                prune_dir <- file.path(data.dir,"tags","r01", 
                    paste("chr",chrs[i],sep=""))
                clump_dir <- file.path(data.dir,"tags","r08",
                    paste("chr",chrs[i], sep=""))
                maf_dir <- file.path(data.dir,"maftssd", 
                    paste("chr",chrs[i], sep=""))
                pv_dir <- file.path(data.dir,"pval",trait,
                    paste("chr",chrs[i], sep=""))
                annot_dir <- file.path(data.dir,"annotation",
                    paste("chr",chrs[i], sep=""))
                excl <- paste(exclude, collapse=",")
                res <- .C("garfield_prep", prune_dir=as.character(prune_dir),
                    clump_dir=as.character(clump_dir),
                    maf_dir=as.character(maf_dir),
                    pv_dir=as.character(pv_dir),
                    annot_dir=as.character(annot_dir),
                    excl=as.character(excl),
                    prep.file=as.character(prep.file))
                }
            }
            if( run.option %in% c("complete", "perm") ){
                ## Run permutation script
                message("\n\n\t[RUN] PERMUTATION STEP")
                num_annot <- length(readLines(file.path(data.dir,"annotation",
                    "link_file.txt")))-1
                pvthresh <- paste(thresh, collapse=",")
                ptthresh <- paste(pt_thresh, collapse=",")
                perm.file <- paste(out.file, ".perm",sep="")
                if(file.exists(perm.file)){
                    file.remove(perm.file)
                    file.create(perm.file)
                }
                res <- .C("garfield_perm",input=as.character(prep.file),
                    link=as.character(file.path(data.dir,"annotation",
                    "link_file.txt")),out_file=as.character(perm.file),
                    p_thresh=as.character(pvthresh),
                    py_thresh=as.character(ptthresh),
                    npermut_par=as.character(nperm),
                    nannot_par=as.character(num_annot),
                    nqmaf_par=as.character(maf.bins),
                    nqntag_par=as.character(tags.bins),
                    nqtssd_par=as.character(tss.bins),
                    optim_mode_par=as.character(optim_mode),
                    minit_par=as.character(minit),
                    thresh=as.character(thresh_perm))
           }
        message("\n--- END COMPUTATION ---\n")
}



