#' @title Epimutations analysis based on outlier detection methods
#' @description The function identifies  differentially methylated regions
#' in a case sample by comparing it against a control panel. 
#' @param case_samples a GenomicRatioSet object containing the case samples.
#' See the constructor function \link[minfi]{GenomicRatioSet}, 
#' \link[minfi]{makeGenomicRatioSetFromMatrix}. 
#' @param control_panel a GenomicRatioSet object containing the 
#' control panel (control panel).
#' @param method a character string naming the 
#' outlier detection method to be used. 
#' This can be set as: \code{"manova"}, 
#' \code{"mlm"}, \code{"iForest"}, \code{"mahdist"}, 
#' \code{"quantile"} and \code{"beta"}. 
#' The default is \code{"manova"}. 
#' For more information see \strong{Details}. 
#' @param chr a character string containing the sequence 
#' names to be analysed. The default value is \code{NULL}. 
#' @param start an integer specifying the start position. 
#' The default value is \code{NULL}.
#' @param end an integer specifying the end position. 
#' The default value is \code{NULL}.
#' @param epi_params the parameters for each method. 
#' See the function \link[epimutacions]{epi_parameters}.  
#' @param maxGap the maximum location gap used in 
#' \link[bumphunter]{bumphunter} method. 
#' @param bump_cutoff a numeric value of the 
#' estimate of the genomic profile above the 
#' cutoff or below the negative of the 
#' cutoff will be used as candidate regions. 
#' @param min_cpg an integer specifying the minimum CpGs number in a DMR.  
# #' @param pca_correction logical. If TRUE methylation PCA correction is 
# #' applied to compensate batch effect. The default value if FALSE
#' @param verbose logical. If TRUE additional details about 
#' the procedure will provide to the user. 
#' The default is TRUE. 
#' @details The function compares a case sample against 
#' a control panel to identify epimutations in the given sample. 
#' First, the DMRs are identified using the 
#' \link[bumphunter]{bumphunter} approach. 
#' After that, CpGs in those DMRs are tested in order to detect regions
#' with CpGs being outliers.  
#' For that, different outlier detection methods can be selected:  
#'  * Multivariate Analysis of Variance (\code{"manova"}). \link[stats]{manova}
#'  * Multivariate Linear Model (\code{"mlm"})
#'  * Isolation Forest (\code{"iForest"}) \link[isotree]{isolation.forest}
#'  * Robust Mahalanobis Distance (\code{"mahdist"}) 
#'  \link[robustbase]{covMcd}
#'  * Quantile distribution (\code{"quantile"})
#'  * Beta (\code{"beta"})
#'  
#' We defined candidate epimutation regions (found in candRegsGR) 
#' based on the 450K array design. As CpGs are not equally distributed 
#' along the genome, only CpGs closer
#' to other CpGs can form an epimutation. 
#' More information can be found in candRegsGR documentation.
#' 
#' @return The function returns an object of class tibble 
#' containing the outliers regions.  
#' The results are composed by the following columns: 
#' * \code{epi_id}: systematic name for each epimutation identified. 
#' It provides the name of the used anomaly detection method. 
#' * \code{sample}: the name of the sample containing the epimutation. 
#' * \code{chromosome}, \code{start} and \code{end}: 
#' indicate the location of the epimutation.
#' * \code{sz}: the window's size of the event.
#' * \code{cpg_n}: the number of CpGs in the epimutation.
#' * \code{cpg_ids}: the names of CpGs in the epimutation.
#' * \code{outlier_score}: 
#'    * For method \code{manova} it provides the approximation 
#'    to F-test and the Pillai score, separated by \code{/}.
#'    * For method \code{mlm} it provides the approximation to 
#'    F-test and the R2 of the model, separated by \code{/}.
#'    * For method \code{iForest} it provides 
#'    the magnitude of the outlier score.
#'    * For method \code{beta} it provides the mean outlier p-value.
#'    * For methods \code{quantile} and 
#'    \code{mahdist} it is filled with NA.
#' * \code{outlier_direction}: indicates the direction 
#' of the outlier with \code{"hypomethylation"} and \code{"hypermethylation"}
#'    * For \code{manova}, \code{mlm}, \code{iForest}, and \code{mahdist} 
#'    it is computed from the values obtained from bumphunter.
#'    * For \code{quantile} it is computed from the location 
#'    of the sample in the reference distribution (left vs. right outlier).
#'    * For method \code{beta} it return a NA.
#' * \code{pvalue}: 
#'    * For methods \code{manova}, \code{mlm}, and \code{iForest} 
#'    it provides the p-value obtained from the model.
#'    * For method \code{quantile}, \code{mahdist} and \code{beta} 
#'    is filled with NA.    
#' * \code{adj_pvalue}: for methods with p-value (\code{manova} and 
#' \code{mlm} adjusted p-value with Benjamini-Hochberg based on the total 
#' number of regions detected by Bumphunter.
#' * \code{epi_region_id}: Name of the epimutation region as defined 
#' in \code{candRegsGR}.
#' * \code{CRE}: cREs (cis-Regulatory Elements) as defined by ENCODE 
#' overlapping the epimutation region. Different cREs are separated by ;.
#' * \code{CRE_type}: Type of cREs (cis-Regulatory Elements) as defined 
#' by ENCODE. Different type are separeted by, 
#' and different cREs are separated by ;.
#' @examples
#' data(GRset)
#' 
#' #Find epimutations in GSM2562701 sample of GRset dataset
#' 
#' case_samples <- GRset[,11]
#' control_panel <- GRset[,1:10]
#' epimutations(case_samples, control_panel, method = "manova")

#' @importFrom minfi annotation getBeta
#' @importFrom GenomicRanges granges makeGRangesFromDataFrame findOverlaps
#' @importFrom stats model.matrix qchisq
#' @importFrom bumphunter bumphunter
#' @importFrom S4Vectors to from
#' @importFrom matrixStats rowQuantiles
#' @import ensembldb

#' 
#' @export
epimutations <- function(case_samples, control_panel,
                        method = "manova", 
                        chr = NULL, start = NULL, end = NULL, 
                        epi_params = epi_parameters(), 
                        maxGap = 1000, bump_cutoff =  0.1, 
                        min_cpg = 3, 
                        # pca_correction = FALSE, 
                        verbose = TRUE)
{
    
    # 1. Inputs check and data extraction
    ## Inputs check
    if (is.null(case_samples)) {
        stop("The argument 'case_samples' must be introduced")
        
    }
    if (is.null(control_panel)) {
        stop("The argument 'case_samples' must be introduced")
    }
    
    if (!is(case_samples, "GenomicRatioSet")) {
        stop("'case_samples' must be of class 'GenomicRatioSet'")
    }
    if (!is(control_panel, "GenomicRatioSet")) {
        stop( "'control_panel' must be of class 'GenomicRatioSet'.
                To create a 'GenomicRatioSet' object use
                'makeGenomicRatioSetFromMatrix' function from minfi package")
    }
    if (minfi::annotation(case_samples)[1] != 
        minfi::annotation(control_panel)[1] &
        minfi::annotation(case_samples)[2] !=
        minfi::annotation(control_panel)[2]) {
        stop("'case_samples' and 'control_panel' annotations must be the same")
    }
    
    if (!is.null(start) & !is.null(end)) {
        if (is.null(chr)) {
            stop("Argument 'chr' must be inroduced with 
                    'start' and 'end' parameters")
        }
        if (length(start) != length(end) & length(chr) != length(start)) {
            stop("'start' and 'end' length must be same")
        }
        
        if (isTRUE(any(start > end))) {
            stop("'start' cannot be higher than 'end'")
        }
        
    }
    if (!is.null(start) & is.null(end) | is.null(start) & !is.null(end)) {
        stop("'start' and 'end' arguments must be introduced together")
        
    }
    
    avail <- c("manova", "mlm", "iForest", "mahdist", "quantile", "beta")
    method <- charmatch(method, avail)
    method <- avail[method]
    if (is.na(method))
        stop("Invalid method was selected'")
    
    if (verbose)
        message("Selected epimutation detection method '", method, "'")
    
    pck <- c("methods", "ensembldb")
    lapply(pck, function(x)
        if (!requireNamespace(x))
            stop("'", x, "'", " package not avaibale"))
    
    # # Apply PCA correction
    # if( pca_correction ) {
    #     if (verbose)
    #         message("Applying PCA correction")
    #     pccorr <- PCA_correction(case_samples, control_samples)
    #     case_samples <- pccorr$cases
    #     control_panel <- pccorr$controls
    #     rm(pccorr)
    # }
    
    
    ## Extract required data:
    #* feature annotation
    #* betas
    #* sample's classification
    
    ### Feature annotation
    fd <- as.data.frame(GenomicRanges::granges(case_samples))
    rownames(fd) <- rownames(case_samples)
    ### Betas
    betas_case <- minfi::getBeta(case_samples)
    betas_case <- betas_case[rownames(fd), , drop = FALSE]
    betas_control <- minfi::getBeta(control_panel)
    betas_control <- betas_control[rownames(fd), ]
    
    ### Select CpGs in the specified in the arguments 'chr', 'start' and 'end'
    
    if (!is.null(chr)) {
        if (!is.null(start) & !is.null(end)) {
            fd <- fd[fd$seqnames %in% chr & fd$start>=start & fd$end <= end, ]
        } else{
            fd <- fd[fd$seqnames %in% chr, ]
        }
        if (nrow(fd) == 0) {
            stop("No CpG was found in the specified region")
        }
        
        betas_case <- betas_case[rownames(fd), , drop = FALSE]
        betas_control <- betas_control[rownames(fd), ]
    }
    
    ### Identify case and control samples
    cas_sam <- colnames(betas_case)
    ctr_sam <- colnames(betas_control)
    
    
    # 2. Epimutations definition (using different methods)
    ##Methods that need bumphunter
    ##("manova", "mlm", "mahdist" and "iForest")
    if (method %in% c("manova", "mlm", "mahdist", "iForest")) {
        if (verbose)
            message("Selected method '", method, "' required of 'bumphunter'")
        # Prepare model to be evaluated
        rst <- do.call(rbind, lapply(cas_sam, function(case) {
            samples_names <- c(ctr_sam, case)
            status <- samples_names == case
            status <- as.data.frame(status, row.names = samples_names)
            betas <- cbind(betas_control, betas_case[, case, drop = FALSE])
            model <- stats::model.matrix( ~ status, status)
            
            # Run bumphunter for region partitioning
            bumps <- bumphunter::bumphunter( object = betas,
                                                design = model,
                                                pos = fd$start,
                                                chr = fd$seqnames,
                                                maxGap = maxGap,
                                                cutoff = bump_cutoff
                                            )$table
            
            if (all(!is.na(bumps))) {
                ## Homogeneize output of bumphunter to epimutacions naming
                bumps$chromosome <- bumps$chr
                bumps$cpg_n <- bumps$L
                bumps$sz <- bumps$end - bumps$start
                
                bumps <- bumps[bumps$L >= min_cpg,]
                if (verbose)
                    message(nrow(bumps),
                            " candidate regions were found for case sample '",
                            case, "'")
                if (nrow(bumps) != 0) {
                    # Identify outliers according to selected method
                    bumps  <- do.call(rbind, lapply(seq_len(nrow(bumps)), 
                        function(ii) {
                            bump <- bumps[ii,]
                            beta_bump <- betas_from_bump(bump, fd, betas)
                            # Add sample name and cpg_ids
                            bump$cpg_ids <- paste(rownames(beta_bump),
                                                    collapse = ",", sep = "")
                            bump$sample <- case
                            bump$delta_beta <- abs(mean(beta_bump[, ctr_sam]) -
                                                    mean(beta_bump[, case]))
                            if (method == "mahdist") {
                                dst <- try(epi_mahdist(beta_bump,
                                    epi_params$mahdist$nsamp), silent = TRUE)
                                if (!is(dst, "try-error")) {
                                    threshold <- sqrt(stats::qchisq(
                                        p = 0.999975,
                                        df = ncol(beta_bump) ))
                                    outliers <- which(dst$statistic>=threshold)
                                    outliers <- dst$ID[outliers]
                                    x <- res_mahdist(case, bump, outliers)
                                }
                            } else if (method == "mlm") {
                                sts <- try(epi_mlm(beta_bump, model), 
                                                            silent = TRUE)
                                x <- res_mlm(bump, sts)
                            } else if (method == "manova") {
                                sts <- epi_manova(beta_bump, model, case)
                                x <- res_manova(bump, sts)
                            } else if (method == "iForest") {
                                sts <- epi_iForest(beta_bump, case,
                                                    epi_params$iForest$ntrees)
                                x <- res_iForest(bump, sts,
                                    epi_params$iForest$outlier_score_cutoff)
                            }
                        }))
                    ## Filter using the adjusted p-value calculated from regions
                    ##identified in each sample ("manova and "mlm")
                    if (method == "manova" & !is.null(bumps)) {
                        bumps <- filter(bumps, epi_params$manova$pvalue_cutoff)
                    }
                    if (method == "mlm" & !is.null(bumps)) {
                        bumps <- filter(bumps, epi_params$mlm$pvalue_cutoff)
                    }
                    
                }
            }
            #Add a row filled by NAs for the samples with any epimutations
            if (is.null(nrow(bumps)) || nrow(bumps) == 0) {
                bumps <- data.frame( chromosome = 0, start = 0, end = 0,
                                sz = NA, cpg_n = NA, cpg_ids = NA,
                                outlier_score = NA, outlier_direction = NA,
                                pvalue = NA, adj_pvalue = NA, delta_beta = NA,
                                sample = case)
            }
            bumps
        }))
        ## Methods that do not need bumphunter ("quantile" and "beta")
    } else if (method == "quantile") {
        # Compute reference statistics
        if (verbose)
            message("Calculating statistics from 'quantile' method")
        if (verbose)
            message("Using quantiles ", epi_params$quantile$qinf, " and ",
                        epi_params$quantile$qsup)
        bctr_prc <- matrixStats::rowQuantiles( betas_control,
                            probs = c(epi_params$quantile$qinf, 
                                        epi_params$quantile$qsup),
                            na.rm = TRUE )
        bctr_pmin <- bctr_prc[, 1]
        bctr_pmax <- bctr_prc[, 2]
        rm(bctr_prc)
        # Run region detection
        rst <- do.call(rbind, lapply(cas_sam, function(case) {
            betas <- cbind(betas_control, betas_case[, case, drop = FALSE])
            x <- epi_quantile( betas_case[, case, drop = FALSE],
                                fd, bctr_pmin, bctr_pmax, ctr_sam,
                                betas, epi_params$quantile$window_sz,
                                min_cpg, epi_params$quantile$offset_abs )
            if (is.null(x) || nrow(x) == 0) {
                x <- data.frame( chromosome = 0, start = 0, end = 0, sz = NA,
                                cpg_n = NA, cpg_ids = NA, outlier_score = NA,
                                outlier_direction = NA, pvalue = NA,
                                adj_pvalue = NA, delta_beta = NA, 
                                sample = case )
            } else {
                x$sample <- case
            }
            x
        }))
    } else if (method == "beta") {
        # Get Beta distribution params
        message("Computing beta distribution parameters")
        beta_params <- getBetaParams(t(betas_control))
        beta_mean <- rowMeans(betas_control, na.rm = TRUE)
        
        message("Defining Regions")
        rst <- do.call(rbind, lapply(cas_sam, function(case) {
            betas <- cbind(betas_control, betas_case[, case, drop = FALSE])
            fd_sub <- fd[rownames(betas_control), ]
            x <- epi_beta( beta_params, beta_mean, 
                            betas_case[, case, drop = FALSE], case, ctr_sam,
                            betas, 
                            GenomicRanges::makeGRangesFromDataFrame(fd_sub),
                            epi_params$beta$pvalue_cutoff,
                            epi_params$beta$diff_threshold, min_cpg, maxGap )
            rm(fd_sub)
            if (nrow(x) == 0) {
                x <- data.frame( chromosome = 0, start = 0, end = 0, sz = NA,
                                cpg_n = NA, cpg_ids = NA, outlier_score = NA,
                                outlier_direction = NA, pvalue = NA,
                                adj_pvalue = NA, delta_beta = NA, 
                                sample = case )
            } else {
                x$sample <- case
            }
            x
        }))
    }
    # 3. Prepare the output and addition of CREs
    ## Prepare the output
    rst$epi_id <- vapply(seq_len(nrow(rst)), function(ii)
                            paste0("epi_", method, "_", ii), character(1))
    rownames(rst) <- seq_len(nrow(rst))
    rst <- rst[, c(13, 12, seq_len(11))]
    
    ## Add CREs and epi_region_id
    rst$CRE_type <- rst$CRE <- rst$epi_region_id <- NA
    rst_c <- rst
    rst_c <- tryCatch({
        rstGR <- GenomicRanges::makeGRangesFromDataFrame(rst)
        ## Ensure chromosomes have the same format
        seqlevelsStyle(rstGR) <- "UCSC"
        #Get candidate regions
        candRegsGR <- get_candRegsGR()
        over <- GenomicRanges::findOverlaps(rstGR, candRegsGR)
        #variables (avoid long code)
        ids <- names(candRegsGR[S4Vectors::to(over)])
        cre <- candRegsGR[S4Vectors::to(over)]$CRE
        cre_type <- candRegsGR[S4Vectors::to(over)]$CRE_type
        
        rst$epi_region_id[S4Vectors::from(over)] <- ids
        rst$CRE[S4Vectors::from(over)] <- cre
        rst$CRE_type[S4Vectors::from(over)] <- cre_type
        rm(c("ids", "cre", "cre_type"))
        rst
    }, error = function(e) {
        rst
    })
    
    ## Convert rst into a tibble class
    if (requireNamespace("tibble", quietly = TRUE)) {
        rst <- tibble::as_tibble(rst_c)
    } else {
        stop("'tibble' package not avaibale")
    }
    return(rst)
}
