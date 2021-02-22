# #!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)
# TadsFile = args[1] # Tab-separated file containing the list of TADs. Each line (no header) should represent a TAD, with genomic coordinates (chr, start, end)
# resolution_kb = as.numeric(args[2]) # Resolution in kb (e.g. '10')
# OutFolder = args[3] # Folder where results should be saved
# ChrSizes_file = args[4] # Tab-separated file containing the chromosome sizes. Each line should correspond to a chromosome (e.g. chr1   249250621)
# ProteinPeaks_Folder = args[5] # Folder containing the three lists of peaks, one for each protein, in .bed format. Each list should be named as 'protein'_peaks.bed (e.g. CTCF_peaks.bed)

#### StructProt_EnrichBoundaries_script.R
# Script to assess the enrichment of CTCF, RAD21 and SMC3 ChIP-seq peaks at TAD boundaries for a given TAD partition.
# Input: see command line arguments.
# Output:
# - Enrich_StructProt_*.pdf: the plot of the average number of peaks every 5 kb around TAD boundaries for the three proteins
# - Enrich_StructProt_*.txt: the table used to generate the plot, with all the values for each protein
# - StructProteins_*.txt: a table containing, for each protein, the ratio of boundaries tagged, the fold change of the peak and its empirical p-value
#
# Refer to Zufferey & Tavernari et al. for details.
# 
# Contact daniele.tavernari@unil.ch for enquiries

####### START ######
start_time <- Sys.time()
####################


require(ggpubr)
require(patchwork)


###### INPUT ######

# OutFolder = paste0(OutFolder,"/")
color1 = "darkgoldenrod2" 
color2 = "royalblue4"
color3 = "orangered3"

mycols <- c("CTCF" = color1, "RAD21" = color2, "SMC3" = color3)

# TadsFile = "data/struct_prot/chr6_TopDom_final_domains.txt"
# resolution_kb = 25
# OutFolder = args[3] # Folder where results should be saved
# ChrSizes_file = args[4] # Tab-separated file containing the chromosome sizes. Each line should correspond to a chromosome (e.g. chr1   249250621)
# 


# res_step = 5
# proteins = c("CTCF", "RAD21", "SMC3")

# Tads_chr <- read.table(TadsFile, quote = '', stringsAsFactors=FALSE, sep = "\t", header = FALSE)
# chr = as.numeric(substr(unique(Tads_chr$V1),4,nchar(unique(Tads_chr$V1))))

# ChrSizes = read.table(file = ChrSizes_file, quote = '', stringsAsFactors=FALSE, sep = "\t", header = FALSE)
# ChrSize = as.numeric(ChrSizes[ChrSizes[,1]==paste0('chr',as.character(chr)),2])

# ChrSize=171075000

###################

############## ENRICHMENT AT BOUNDARY SITES ##############
# Final_Ratios_df <- data.frame(TadsFile = TadsFile, resolution_kb = NA, protein = NA, domains_ratio = NA)
# Final_EnrichPeaks_df <- as.data.frame(matrix(nrow = 0, ncol = (as.integer(1000/res_step)+3) ))
# colnames(Final_EnrichPeaks_df) <- c("TadsFile", "resolution_kb", "protein", as.character(c(1:as.integer(1000/res_step))))
# 
# 
# Allprots_enrichs <- data.frame(row.names = c(1:as.integer(1000/res_step)))

# default_prot <- c("CTCF", "RAD21", "SMC3")

# tad_dt = Tads_chr

get_structEnrichmentPlot <- function(tad_dt, resolution_kb, 
                                     proteins,
                                     # proteins = default_prot, 
                                     # ProteinPeaks_Folder = "data/struct_prot/",
                                     # ProteinPeaks_Folder,
                                     ProteinPeaks_Files,  # named list
                                     res_step = 5, 
                                     plot_tit=NULL) {
   
   
   stopifnot(proteins %in% names(ProteinPeaks_Files))
   
   print(paste0("plot_tit = ", plot_tit))
   
   # Final_Ratios_df <- data.frame(TadsFile = TadsFile, resolution_kb = NA, protein = NA, domains_ratio = NA)
   Final_Ratios_df <- data.frame(TadsFile = "foo", resolution_kb = NA, protein = NA, domains_ratio = NA)
   Final_EnrichPeaks_df <- as.data.frame(matrix(nrow = 0, ncol = (as.integer(1000/res_step)+3) ))
   colnames(Final_EnrichPeaks_df) <- c("TadsFile", "resolution_kb", "protein", as.character(c(1:as.integer(1000/res_step))))
   
   
   Allprots_enrichs <- data.frame(row.names = c(1:as.integer(1000/res_step)))
   
   
   Tads_chr <- tad_dt
   
   head(Tads_chr)
   
   # chr = as.numeric(substr(unique(Tads_chr$V1),4,nchar(unique(Tads_chr$V1))))
   # change MZ: don't use as.numeric because won't work for chrX
   chr = as.character(substr(unique(Tads_chr$V1),4,nchar(unique(Tads_chr$V1))))
   stopifnot(!is.na(chr))
   stopifnot(length(chr) == 1)
   plot_chr <- paste0("chr", chr)
   
   ChrSize <- max(as.numeric(as.character(Tads_chr$V3)))
   stopifnot(ChrSize >= as.numeric(as.character(Tads_chr$V2)))
   
   head(Tads_chr)
   # plot(x=Tads_chr$V3, y =Tads_chr$V2, main=paste0(ChrSize))
   # return()
   
   stopifnot(proteins %in% default_prot)
   for (protein in proteins)
   {
      
      print(paste0("protein = ", protein, "\n"))
      # Peaks <- read.table(paste0(ProteinPeaks_Folder, "chr6_", protein,'_peaks.bed'), quote = '', stringsAsFactors=FALSE, sep = "\t", header = FALSE)
      peakFile <- ProteinPeaks_Files[[paste0(protein)]]
      stopifnot(length(peakFile) > 0 & !is.null(peakFile))
      print(paste0("... reading: ",peakFile))
      # Peaks <- read.table(paste0(ProteinPeaks_Folder, "chr6_", protein,'_peaks.bed'), quote = '', stringsAsFactors=FALSE, sep = "\t", header = FALSE)
      Peaks <- read.table(peakFile, quote = '', stringsAsFactors=FALSE, sep = "\t", header = FALSE)
      
      Peaks_chr <- Peaks[Peaks[,1]==paste0("chr",as.character(chr)),c(2,3)]
      
      Enrich_df = data.frame(matrix(vector(), 0, as.integer(1000/res_step)))
      
      Tad_edges <- unique(c(Tads_chr[,2], (Tads_chr[,3]+1)))
      
      Tad_edges=Tad_edges[1:10] ################# !!!!!!!!!!!!!!!!!!!!!!!!!!!!! COMMENT HERE FOR FINAL VERSION <<<<<<<<<<<<<<<<<
      
      for (tad_edge in Tad_edges)
      {
         if (((tad_edge-500*1000)<0) | ((tad_edge+500*1000)>ChrSize)) next
         Tad_surroundings <- seq(from = tad_edge-500*1000, to = tad_edge+500*1000, by = res_step*1000)
         Peaks_tad <- Peaks_chr[Peaks_chr[,2]>=(tad_edge-500*1000) & Peaks_chr[,1]<=(tad_edge+500*1000),]
         Vec_tad <- vector(mode = "numeric", length = as.integer(1000/res_step))
         for (pt in 1:dim(Peaks_tad)[1])
         {
            query <- as.numeric(Peaks_tad[pt,])
            query <- query[(query>Tad_surroundings[1]) & (query<Tad_surroundings[length(Tad_surroundings)]) ]
            a <- hist(query, breaks = Tad_surroundings, plot = FALSE)
            Vec_tad <- Vec_tad + as.numeric(a$counts>0)
         }
         Enrich_df <- rbind(Enrich_df, Vec_tad)
      }
      Ratio_tagged <- sum(as.numeric(rowSums(Enrich_df[,c((length(Vec_tad)/2-resolution_kb/res_step+1):(length(Vec_tad)/2+resolution_kb/res_step))])>0))/dim(Enrich_df)[1]
      Allprots_enrichs[,protein] <- colMeans(Enrich_df)
      # Ratios_temp_df <- data.frame(TadsFile = TadsFile, resolution_kb = resolution_kb, protein = protein, domains_ratio = Ratio_tagged)
      Ratios_temp_df <- data.frame(TadsFile = "foo", resolution_kb = resolution_kb, protein = protein, domains_ratio = Ratio_tagged)
      Final_Ratios_df <- rbind(Final_Ratios_df, Ratios_temp_df)
      
   }
   
   # pdf(paste0(OutFolder, "Enrich_StructProt_chr",chr,"_res",resolution_kb,"kb.pdf"))
   ymax <- max(Allprots_enrichs) 
   ymax = max(ymax, 0.14)
   ymin <- min(Allprots_enrichs) 
   # plot_df = data.frame(row.names = c(1:as.integer(1000/res_step)), x = c(1:as.integer(1000/res_step)),
   #                      CTCF = Allprots_enrichs[,"CTCF"], 
   #                      RAD21 = Allprots_enrichs[,"RAD21"], 
   #                      SMC3 = Allprots_enrichs[,"SMC3"] )
   
   tmp_dt2 <- Allprots_enrichs[,proteins, drop=FALSE]
   colnames(tmp_dt2) <- proteins
   tmp_dt1 = data.frame(row.names = c(1:as.integer(1000/res_step)), x = c(1:as.integer(1000/res_step)))
   plot_dt <- cbind(tmp_dt1, tmp_dt2)
   
   main_tit <- ifelse(is.null(plot_tit),paste0(ds_name, " - chr",chr,", resolution ",resolution_kb,"kb") , plot_tit)
   
   print(paste0("main_tit = ", main_tit))
                        
   # write.table(plot_df, file = paste0(OutFolder, "Enrich_StructProt_chr",chr,"_res",resolution_kb,"kb.txt"), row.names = F, col.names = T, quote = F )
   plot(x=c(1:as.integer(1000/res_step)), 
        # y=Allprots_enrichs[,"CTCF"],
        y=Allprots_enrichs[,proteins[1]],
        type = "l", col = color1, main = main_tit, 
        xlab = "Distance from TAD boundary",
        ylab = paste0("Average # of peaks / ", as.character(res_step)," kb"), 
        xaxt='n', cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
   axis(1, at=seq(from = 0, to = as.integer(1000/res_step), by = as.integer(1000/res_step)/4), 
        labels=c("-500 kb","-250 kb", "0", "+250 kb", "+500 kb"), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
   if(length(proteins) > 1) {
      lines(x=c(1:as.integer(1000/res_step)),
            # y=Allprots_enrichs[,"RAD21"],
            y=Allprots_enrichs[,proteins[2]],
            col = color2)
      all_cols <- c(color1, color2)
   } else {
      all_cols <- color1
   }
   if(length(proteins) == 3){
      lines(x=c(1:as.integer(1000/res_step)), 
            y=Allprots_enrichs[,proteins[3]], col = color3)
      all_cols <- c(color1, color2, color3)
   }
   # legend(x="topright", legend = proteins, col = c(color1, color2, color3), lty = 1, cex = 1.5)
   legend(x="topright", legend = proteins, col = all_cols, lty = 1, cex = 1.5)
   mtext(side=3, text=paste0(plot_chr,"; Hi-C resolution (kb) = ", resolution_kb,"; agg. resol. step = ",res_step))
}

#############################################################################################
############################################################################################# FC AND PVAL COMPARISON
#############################################################################################
mean_fc_over_bg <- function( profile, na.rm = TRUE, ratio_from_peak = 0.01, ratio_from_end = 0.6 ){
   if (na.rm) profile <- profile[!is.na(profile)]
   peak_val <- mean(profile[ (as.integer(length(profile)/2)-as.integer(length(profile)*ratio_from_peak)) : (as.integer(length(profile))/2+as.integer(length(profile)*ratio_from_peak)) ])
   bg <- mean(profile[c( 1:(as.integer(length(profile)/2)-as.integer(length(profile)*(1-ratio_from_end))), (as.integer(length(profile))/2+as.integer(length(profile)*(1-ratio_from_end))):length(profile) )])

   fc_over_bg <- peak_val/bg - 1
   return(fc_over_bg)
}

mean_pval_of_peak <- function( profile, plot_dist = FALSE, na.rm = TRUE, ratio_from_peak = 0.01, ratio_from_end = 0.6 ){
   if (na.rm) profile <- profile[!is.na(profile)]
   peak_val <- mean(profile[ (as.integer(length(profile)/2)-as.integer(length(profile)*ratio_from_peak)) : (as.integer(length(profile))/2+as.integer(length(profile)*ratio_from_peak)) ])
   bg <- profile[c( 1:(as.integer(length(profile)/2)-as.integer(length(profile)*(1-ratio_from_end))), (as.integer(length(profile))/2+as.integer(length(profile)*(1-ratio_from_end))):length(profile) )]

   if (plot_dist)
   {
      hist(bg, 30, xlim = c(min(bg), max(peak_val, max(bg) )),  main = "Distribution of background and value of peak")
      abline(v = peak_val, col = "red")
   }

   pval_of_peak <- pnorm( peak_val, mean = mean(bg), sd = sd(bg), lower.tail = FALSE )
   return(pval_of_peak)
}

get_multiStructEnrichmentPlot <- function(all_tad_dt, resolution_kb,
                                  proteins,
                                  # proteins = default_prot, 
                                  # ProteinPeaks_Folder = "data/struct_prot/",
                                  ProteinPeaks_Files, # named list to datapath
                                  res_step = 5, 
                                     plot_tit=NULL, 
                                  regexsub=NULL
                                  ) {

   
   stopifnot(proteins %in% names(ProteinPeaks_Files))
   
   print(paste0("plot_tit multi = ", plot_tit))
   
   
   
   # plot(x=Tads_chr$V3, y =Tads_chr$V2, main=paste0(ChrSize))
   # return()
   
   
   
   if(is.null(names(all_tad_dt))){
      names(all_tad_dt) <- paste0("dataset_", seq_along(all_tad_dt))
   }
   
   nProfileCol <- as.integer(1000/res_step)+3
   
   
   # Final_Ratios_df <- data.frame(TadsFile = TadsFile, resolution_kb = NA, protein = NA, domains_ratio = NA)
   Final_Ratios_df <- data.frame(
      TadsFile = "foo",
      # TadsFile = dsname,
      resolution_kb = NA, protein = NA, domains_ratio = NA)
   Final_EnrichPeaks_df <- as.data.frame(matrix(nrow = 0, ncol = nProfileCol ))
   colnames(Final_EnrichPeaks_df) <- c("TadsFile", "resolution_kb", "protein", as.character(c(1:as.integer(1000/res_step))))
   
   
    
   Final_Ratios_df$fc_over_bg <- NA
   Final_Ratios_df$pval_of_peak <- NA
   
   all_chrs <- list()
   
   for(i in seq_along(all_tad_dt)) {
      
      tad_dt <- all_tad_dt[[i]]
      dsname <- names(all_tad_dt)[i]
      
      Tads_chr <- tad_dt
      
      head(Tads_chr)
      
      # chr = as.numeric(substr(unique(Tads_chr$V1),4,nchar(unique(Tads_chr$V1))))
      # change MZ: not use as numeric because won't work if chrX
      chr = as.character(substr(unique(Tads_chr$V1),4,nchar(unique(Tads_chr$V1))))
      stopifnot(!is.na(chr))
      
      all_chrs[[i]] <- paste0("chr", chr)
      
      ChrSize <- max(as.numeric(as.character(Tads_chr$V3)))
      stopifnot(ChrSize >= as.numeric(as.character(Tads_chr$V2)))
      
      head(Tads_chr)
      
      save(Tads_chr, file="Tads_chr.Rdata", version=2)
      
      Allprots_enrichs <- data.frame(row.names = c(1:as.integer(1000/res_step)))
      
      for (protein in proteins)
      {
         
         print(paste0("> START - ", dsname, " - ", protein))
         
         peakFile <- ProteinPeaks_Files[[paste0(protein)]]
         stopifnot(length(peakFile) > 0 & !is.null(peakFile))
         print(paste0("... reading: ",peakFile))
         # Peaks <- read.table(paste0(ProteinPeaks_Folder, "chr6_", protein,'_peaks.bed'), quote = '', stringsAsFactors=FALSE, sep = "\t", header = FALSE)
         Peaks <- read.table(peakFile, quote = '', stringsAsFactors=FALSE, sep = "\t", header = FALSE)
         
         Peaks_chr <- Peaks[Peaks[,1]==paste0("chr",as.character(chr)),c(2,3)]
         Enrich_df = data.frame(matrix(vector(), 0, as.integer(1000/res_step)))
         
         Tad_edges <- unique(c(Tads_chr[,2], (Tads_chr[,3]+1)))
         
         Tad_edges=Tad_edges[1:10] ################# !!!!!!!!!!!!!!!!!!!!!!!!!!!!! COMMENT HERE FOR FINAL VERSION <<<<<<<<<<<<<<<<<
         
         for (tad_edge in Tad_edges)
         {
            if (((tad_edge-500*1000)<0) | ((tad_edge+500*1000)>ChrSize)) next
            Tad_surroundings <- seq(from = tad_edge-500*1000, to = tad_edge+500*1000, by = res_step*1000)
            Peaks_tad <- Peaks_chr[Peaks_chr[,2]>=(tad_edge-500*1000) & Peaks_chr[,1]<=(tad_edge+500*1000),]
            Vec_tad <- vector(mode = "numeric", length = as.integer(1000/res_step))
            for (pt in 1:dim(Peaks_tad)[1])
            {
               query <- as.numeric(Peaks_tad[pt,])
               query <- query[(query>Tad_surroundings[1]) & (query<Tad_surroundings[length(Tad_surroundings)]) ]
               a <- hist(query, breaks = Tad_surroundings, plot = FALSE)
               Vec_tad <- Vec_tad + as.numeric(a$counts>0)
            }
            Enrich_df <- rbind(Enrich_df, Vec_tad)
         }
         Ratio_tagged <- sum(as.numeric(rowSums(Enrich_df[,c((length(Vec_tad)/2-resolution_kb/res_step+1):(length(Vec_tad)/2+resolution_kb/res_step))])>0))/dim(Enrich_df)[1]
         Allprots_enrichs[,protein] <- colMeans(Enrich_df)
         Ratios_temp_df <- data.frame(
                                    # TadsFile = TadsFile, 
                                    TadsFile = dsname,
                                      # resolution_kb = resolution_kb, protein = protein, domains_ratio = Ratio_tagged)
                                    resolution_kb = resolution_kb, protein = protein, domains_ratio = Ratio_tagged,
                                    fc_over_bg = NA,pval_of_peak=NA
                                    )
         Final_Ratios_df <- rbind(Final_Ratios_df, Ratios_temp_df)
         save(Final_Ratios_df, file=paste0( dsname, "_Final_Ratios_df.Rdata"), version=2)
         print(paste0("... saved: ", dsname, "_Final_Ratios_df.Rdata"))
         # TadsFile resolution_kb protein domains_ratio
         # 1 1_arrowhead_100kb_ICE_final_domains.txt            NA    <NA>            NA
         # 2 1_arrowhead_100kb_ICE_final_domains.txt            25    CTCF             1
         # 
         Final_EnrichPeaks_temp_annotation_df <- data.frame(
                                                            # TadsFile = TadsFile, 
                                                            TadsFile = dsname, 
                                                            resolution_kb = resolution_kb, 
                                                            # protein = proteins)
                                                                protein = protein)
         # Final_EnrichPeaks_temp_values_df <- as.data.frame(t(Allprots_enrichs))
         Final_EnrichPeaks_temp_values_df <- as.data.frame(t(Allprots_enrichs[,protein]))
         Final_EnrichPeaks_temp_df <- cbind(Final_EnrichPeaks_temp_annotation_df, Final_EnrichPeaks_temp_values_df)
         Final_EnrichPeaks_df <- rbind(Final_EnrichPeaks_df, Final_EnrichPeaks_temp_df)
         print(paste0("dim Final_EnrichPeaks_df = ", dim(Final_EnrichPeaks_df)))
         
         save(Final_EnrichPeaks_temp_df, file=paste0( dsname, "_Final_EnrichPeaks_temp_df.Rdata"), version=2)
         print(paste0("... saved: ", dsname, "_Final_EnrichPeaks_temp_df.Rdata"))
         # TadsFile resolution_kb protein 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
         # CTCF 1_arrowhead_100kb_ICE_final_domains.txt            25    CTCF 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
         # 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39  40 41 42 43 44 45 46 47 48
         # CTCF  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 0.1  0  0  0  0  0  0  0  0
         
         
         stopifnot(ncol(Final_EnrichPeaks_df) == nProfileCol)
         
            profile <- Final_EnrichPeaks_df[ (
                                    # (Final_EnrichPeaks_df$TadsFile==TadsFile) &
                                       (Final_EnrichPeaks_df$TadsFile==dsname) &
                                                 (Final_EnrichPeaks_df$resolution_kb==resolution_kb) & 
                                                 (Final_EnrichPeaks_df$protein==protein) ) , 
                                             # c(4:203)]
                                                   c(4:nProfileCol)]
            
            stopifnot(length(profile) ==  as.integer(1000/res_step))
            
            fc_over_bg_val <- mean_fc_over_bg(profile)
            pval_of_peak_val <- mean_pval_of_peak(profile, plot_dist = FALSE)

            Final_Ratios_df$fc_over_bg[ 
                                       # Final_Ratios_df$TadsFile==TadsFile & 
                                          Final_Ratios_df$TadsFile==dsname & 
                                           Final_Ratios_df$resolution_kb==resolution_kb &
                                           Final_Ratios_df$protein==protein ] <- fc_over_bg_val
            Final_Ratios_df$pval_of_peak[ 
                                             # Final_Ratios_df$TadsFile==TadsFile & 
                                                Final_Ratios_df$TadsFile==dsname & 
                                             Final_Ratios_df$resolution_kb==resolution_kb & 
                                             Final_Ratios_df$protein==protein ] <- pval_of_peak_val
         
         
         
      }   # end protein   
   } # end dataset

   save(Final_EnrichPeaks_df, file=paste0( "Final_EnrichPeaks_df.Rdata"), version=2)
   print("... saved: Final_EnrichPeaks_df.Rdata")
   
   save(Final_Ratios_df, file=paste0( "Final_Ratios_df.Rdata"), version=2)
   print("... saved: Final_Ratios_df.Rdata")
   
   print(paste0("length(all_tad_dt) = ", paste0(length(all_tad_dt), collapse=",")))
   print(paste0("length(proteins) = ", paste0(length(proteins), collapse=",")))
   
   print(paste0("dim(Final_EnrichPeaks_df) = ", paste0(dim(Final_EnrichPeaks_df), collapse=",")))
   stopifnot(nrow(Final_EnrichPeaks_df) == length(all_tad_dt)*length(proteins))
   stopifnot(ncol(Final_EnrichPeaks_df) == nProfileCol)
   
   print(paste0("done for:\n", paste0( names(all_tad_dt), collapse = " - "),"\n",paste0( proteins, collapse = " - ") ))
   
   plot_dt <- Final_Ratios_df[Final_Ratios_df$TadsFile %in% names(all_tad_dt),]
   save(plot_dt, file="plot_dt.Rdata", version=2)
   print(paste0("... written: ", "plot_dt.Rdata"))
   
   if(!is.null(regexsub) & regexsub != ""){
      plot_dt$TadsFile <- gsub(regexsub, "", plot_dt$TadsFile)
   }
   save(plot_dt, file="plot_dt2.Rdata", version=2)
   print(paste0("... written: ", "plot_dt2.Rdata"))
   
   nona_dt <- plot_dt[!is.na(plot_dt$fc_over_bg),]
   na_dt <- plot_dt[is.na(plot_dt$fc_over_bg),]
   na_levels <- sort(unique(as.character(na_dt$TadsFile)))

   tmp_dt <- aggregate(fc_over_bg~TadsFile, data=nona_dt, FUN=mean)
   tmp_dt <- tmp_dt[sort(tmp_dt$fc_over_bg, decreasing = TRUE),]

   plot_dt$TadsFile <- factor(as.character(plot_dt$TadsFile), levels=c(as.character(tmp_dt$TadsFile), na_levels))
   stopifnot(!is.na(plot_dt$TadsFile))
   
   myxlab <- ""
   
   plot_chr <- paste0(unique(unlist(all_chrs)), collapse = ", ")
   
   p1 <- ggbarplot(plot_dt,   
                   position = position_dodge(0.9),
                   x="TadsFile", 
                   y="fc_over_bg", 
                   xlab=myxlab, ylab="FC over background", color="protein", fill="protein") +
      # geom_bar(stat="count", position="dodge2")+
      scale_color_manual(values=mycols) +
      scale_fill_manual(values=mycols) +
      labs(color="", fill="")

   p2 <- ggbarplot(plot_dt, 
                   position = position_dodge(0.9),
                   x="TadsFile", y="pval_of_peak", 
                   xlab=myxlab, ylab="Peak p-val", color="protein", fill="protein") +
      scale_color_manual(values=mycols) +
      scale_fill_manual(values=mycols) + 
      labs(color="", fill="")

   patchwork <- (p1 + p2)
   
   patchwork + 
      plot_layout(guides = 'collect')+
      plot_annotation(
      title = plot_tit,
      subtitle = paste0(plot_chr, "; agg. resol. = ", res_step, " kb; Hi-C resol. = ", resolution_kb, " kb") 
   ) & theme(plot.title=element_text(hjust=0.5, size=14, face = "bold"), 
             plot.subtitle=element_text(hjust=0.5, size=12, face = "italic"),
             legend.position = 'top')
   
   
}



# 
# # pdf(paste0(OutFolder, "Enrich_StructProt_chr",chr,"_res",resolution_kb,"kb.pdf"))
# ymax <- max(Allprots_enrichs) 
# ymax = max(ymax, 0.14)
# ymin <- min(Allprots_enrichs) 
# plot_df = data.frame(row.names = c(1:as.integer(1000/res_step)), x = c(1:as.integer(1000/res_step)), CTCF = Allprots_enrichs[,"CTCF"], RAD21 = Allprots_enrichs[,"RAD21"], SMC3 = Allprots_enrichs[,"SMC3"] )
# # write.table(plot_df, file = paste0(OutFolder, "Enrich_StructProt_chr",chr,"_res",resolution_kb,"kb.txt"), row.names = F, col.names = T, quote = F )
# plot(x=c(1:as.integer(1000/res_step)), y=Allprots_enrichs[,"CTCF"], type = "l", col = color1, main = paste0("chr",chr,", resolution ",resolution_kb,"kb"), xlab = "Distance from TAD boundary", ylab = paste0("Average # of peaks / ", as.character(res_step)," kb"), xaxt='n', cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
# axis(1, at=seq(from = 0, to = as.integer(1000/res_step), by = as.integer(1000/res_step)/4), labels=c("-500 kb","-250 kb", "0", "+250 kb", "+500 kb"), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
# lines(x=c(1:as.integer(1000/res_step)), y=Allprots_enrichs[,"RAD21"], col = color2)
# lines(x=c(1:as.integer(1000/res_step)), y=Allprots_enrichs[,"SMC3"], col = color3)
# legend(x="topright", legend = proteins, col = c(color1, color2, color3), lty = 1, cex = 1.5)
# dev.off()
# 
# Final_EnrichPeaks_temp_annotation_df <- data.frame(TadsFile = TadsFile, resolution_kb = resolution_kb, protein = proteins)
# Final_EnrichPeaks_temp_values_df <- as.data.frame(t(Allprots_enrichs))
# Final_EnrichPeaks_temp_df <- cbind(Final_EnrichPeaks_temp_annotation_df, Final_EnrichPeaks_temp_values_df)
# Final_EnrichPeaks_df <- rbind(Final_EnrichPeaks_df, Final_EnrichPeaks_temp_df)
# 
# ##########################################################
# 
# ################# QUANTIFYING ENRICHMENT #################
# 

# 
# Final_Ratios_df$fc_over_bg <- NA
# Final_Ratios_df$pval_of_peak <- NA
# 
# for (protein in proteins)
# {
#    profile <- Final_EnrichPeaks_df[ ((Final_EnrichPeaks_df$TadsFile==TadsFile) & (Final_EnrichPeaks_df$resolution_kb==resolution_kb) & (Final_EnrichPeaks_df$protein==protein) ) , c(4:203)]
# 
#    fc_over_bg_val <- mean_fc_over_bg(profile)
#    pval_of_peak_val <- mean_pval_of_peak(profile, plot_dist = FALSE)
# 
#    Final_Ratios_df$fc_over_bg[ Final_Ratios_df$TadsFile==TadsFile & Final_Ratios_df$resolution_kb==resolution_kb & Final_Ratios_df$protein==protein ] <- fc_over_bg_val
#    Final_Ratios_df$pval_of_peak[ Final_Ratios_df$TadsFile==TadsFile & Final_Ratios_df$resolution_kb==resolution_kb & Final_Ratios_df$protein==protein ] <- pval_of_peak_val
# }
# Final_Ratios_df = Final_Ratios_df[!is.na(Final_Ratios_df$protein),]
# write.table(Final_Ratios_df, file = paste0( OutFolder, "StructProteins_chr",chr,"_res",resolution_kb,"kb.txt" ), quote = FALSE, sep = "\t", row.names = FALSE)
# 
# ##########################################################
# 
# ####### END #######
# end_time <- Sys.time()
# time_taken <- end_time - start_time
# print(time_taken)
###################