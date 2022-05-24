#' @title Identifies epimutations using quantile distribution
#' @description  Identifies CpGs with 
#' outlier methylation values
#' using a sliding window
#' approach to compare individual 
#' methylation profiles of a single
#' case sample against all other 
#' samples from reference panel (controls)
#' @param case beta values for a 
#' single case (data.frame). The samples as 
#' single column and CpGs in rows (named).
#' @param fd feature description as 
#' data.frame having at least chromosome and
#' position as columns and and CpGs in rows (named).
#' @param bctr_pmin Beta value observed 
#' at 0.01 quantile in controls. A beta
#' values has to be lower or equal 
#' to this value to be considered an 
#' epimutation.
#' @param bctr_pmax Beta value observed 
#' at 0.99 quantile in controls. A beta
#' values has to be higher or equal 
#' to this value to be considered an 
#' epimutation.
#' @param window_sz Maximum distance 
#' between a pair of CpGs to defined an
#' region of CpGs as epimutation (default: 1000).
#' @param betas a matrix containing the 
#' beta values for all samples.
#' @param controls control samples names.
#' @param N Minimum number of CpGs, 
#' separated in a maximum of window_sz bass,
#' to defined an epimutation (default: 3).
#' @param offset_abs Extra enforcement 
#' defining an epimutation based on 
#' beta values at 0.005 and 0.995 quantiles (default: 0.15).
#' @return The function returns a data 
#' frame with the regions candidates to be
#' epimutations.
epi_quantile <- function(case, 
                         fd, 
                         bctr_pmin, 
                         bctr_pmax, 
                         controls, betas,
                         window_sz = 1000, 
                         N = 3, offset_abs = 0.15) {
  # Check that there is a single proband
  if(ncol(case) != 1) {
    stop("Epimutation detection with 'quantile'
         approach can only works with a singe proband")
  }
  # Check that "N" is not smaller than 3
  if(N < 3){
    stop("The minimum number of CpGs allowed is 3")
  }
  
  if (!requireNamespace("methods")) 
    stop("'methods' package not available")
  
  flag_result <- data.frame(
    chr = as.character(fd[rownames(case), "seqnames"]),
    pos = fd[rownames(case), "start"],
    flag_qm_sup = case[, 1] >= bctr_pmax + offset_abs,
    flag_qm_inf = case[, 1] <= bctr_pmin - offset_abs,
    stringsAsFactors = FALSE
  )
  
  flag_result <- flag_result[!is.na(flag_result$flag_qm_sup) & 
                             !is.na(flag_result$flag_qm_inf), ]
 

  
  # Function used to detect regions of N CpGs closer than window size
  get_regions <- function(flag_df, window_sz = 1000, N = 3, pref = "R") {
   
    if(nrow(flag_df) < N) {
      return(data.frame(chr=NA, pos=NA, region=NA))
    }
    
    # Order the input first by chromosome and then by position
    flag_df <- flag_df[with(flag_df, order(chr, pos)), ]
    
    x <- do.call(rbind, lapply(unique(flag_df$chr), function(chr) {
      # Subset by chromosome
      red_df <- flag_df[flag_df$chr == chr, ]
      if(nrow(red_df) < N) {
        return(data.frame(chr=NA, pos=NA, region=NA))
      }
      
      # Get the position of the previous and 
      #next probe for each probe in the
      # data.frame. The first and last position get its 
      #own position minus/plus
      # the window size to be sure to include 
      #them in the resulting data.frame.
      red_df$pos_next <- c(red_df$pos[seq(2, nrow(red_df))], 
                           red_df$pos[nrow(red_df)] + window_sz + 1)
      red_df$pos_prev <- c(red_df$pos[1] - window_sz - 1, 
                           red_df$pos[seq(1, nrow(red_df) - 1)])
      
      # We add two columns indicating if 
      # a probe is within the window size
      # range with its previous and with its next probe
      red_df$in_prev <- red_df$pos - red_df$pos_prev <= window_sz
      red_df$in_next <- red_df$pos_next - red_df$pos <= window_sz
      
      # We drop all the probes that do not 
      #have a previous nor a next
      # probe within the range of interest
      red_df <- red_df[red_df$in_prev | red_df$in_next, ]
      if(nrow(red_df) < N) {
        return(data.frame(chr=NA, pos=NA, region=NA))
      }
      
      # Using the cumsum function and by 
      # negating the content of the "in_next"
      # column we can define the regions 
      # of CpGs within the range since they
      # will be tagged with the same number
      red_df$cum <- cumsum(!red_df$in_next)
      
      # Correct the base position of 
      # the change in the region
      red_df$cum2 <- c(red_df$cum[1], 
                       vapply(seq(2, nrow(red_df)), function(ii) {
        if(red_df$cum[ii] != red_df$cum[ii - 1] & 
           red_df$in_prev[ii] & !red_df$in_next[ii]) {
          return (red_df$cum[ii] - 1)
        } else {
          return (red_df$cum[ii])
        }
      }, numeric(1)))

      # Computing the frequency of each 
      # "number" assign to the region we can 
      # know how may probes are in it. 
      # We can use this frequency to filter out
      # those regions with less probes than given by N.
      # We also give to the regions a proper name.
      fr <- data.frame(table(red_df$cum2), 
                       stringsAsFactors = FALSE)
      fr <- as.numeric(as.character(fr$Var1[fr$Freq >= N]))
      if(length(fr) > 0) {
        fr <- data.frame(current = fr, 
                         new = paste0(pref,
                                      "_", chr, 
                                      "_", seq_len(length(fr))))
        red_df <- red_df[red_df$cum2 %in% fr$current, ]
        rownames(fr) <- paste0("O", fr$current)
        red_df$region <- fr[paste0("O", red_df$cum2), "new"]
        
        
        # Since the first and last probe in 
        # a chromosome will have TRUE in
        # prev or next distance we need to 
        # be sure to drop them if they
        # are not in the window
        red_df$dist_next <- red_df$pos_next - red_df$pos
        red_df$dist_next[length(red_df$dist_next)] <- 0
        red_df <- red_df[red_df$dist_next <= window_sz, ]
        
        # We drop the columns with the flags 
        # used for the outlier and region detection
        red_df <- red_df[ , c("chr", "pos", "region")]
        return(red_df)
      } else {
        return(data.frame(chr=NA, pos=NA, region=NA))
      }
    }))
    
    return(x[!is.na(x$chr), ])
  }
  
  # We select the CpGs according to the direction of the outlier
  flag_sup <- flag_result[flag_result$flag_qm_sup, ]
  flag_inf <- flag_result[flag_result$flag_qm_inf, ]
  
  # We identify the regions taking into account the direction
  reg_sup <- get_regions(flag_sup, pref = "Rs")
  reg_inf <- get_regions(flag_inf, pref = "Ri")
  
  # We add a column indicating the direction of the regions/outliers
  if(nrow(reg_sup) != 0){
    reg_sup$outlier_direction <- "hypermethylation"
    reg_sup$CpG_ids <- rownames(reg_sup) 
  }
  if(nrow(reg_inf) != 0){
    reg_inf$outlier_direction <- "hypomethylation"
    reg_inf$CpG_ids <- rownames(reg_inf)
  }


  collapse_regions <- function(flag_df, 
                               case_name, 
                               controls, 
                               betas) {
    empty <- data.frame(
      chromosome = character(),
      start = numeric(),
      end = numeric(),
      sz = numeric(),
      cpg_n = numeric(),
      cpg_ids = character(),
      outlier_score = numeric(),
      outlier_direction = character(),
      pvalue = numeric(),
      adj_pvalue = numeric(),
      delta_beta = numeric()
    )
    if(nrow(flag_df) == 0) {
      return(empty)
    }
    if(all(is.na(flag_df))){
      return(empty)
    }
    do.call(rbind, lapply(unique(flag_df$region), 
                          function(reg) {
      x <- flag_df[flag_df$region == reg, ]
      if(nrow(x) > 0) {
        data.frame(
          chromosome = x$chr[1],
          start = min(x$pos),
          end = max(x$pos),
          sz = max(x$pos) - min(x$pos),
          cpg_n = nrow(x),
          cpg_ids = paste(x$CpG_ids, 
                          collapse = ",", 
                          sep = ""),
          outlier_score = NA,
          outlier_direction = x$outlier_direction[1],
          pvalue = NA,
          adj_pvalue = NA,
          delta_beta = abs(mean(betas[x$CpG_ids, controls]) -
                           mean(betas[x$CpG_ids, colnames(case)]))
        )
      } else {
        empty
      }
    }))
  }
  
  # We collapse the CpGs in regions and format the output
  clean_sup <- collapse_regions(reg_sup, 
                                colnames(case), 
                                controls, 
                                betas)
  clean_inf <- collapse_regions(reg_inf, 
                                colnames(case), 
                                controls, 
                                betas)
  
  rst <- rbind(clean_inf, clean_sup)
  rst <- rst[!is.na(rst$chromosome), ]
  return(rst)
}
