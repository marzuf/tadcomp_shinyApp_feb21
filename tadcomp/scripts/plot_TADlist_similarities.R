# plot_TADlist_similarities.R

plot_TADlist_comparison <- function(TAD_list, metric, lowColor="blue", highColor="red", nCpu=1, ...) {
  symmMetric <- TRUE
  require(doMC)
  require(foreach)
  require(ggplot2)
  registerDoMC(nCpu)
  metric_args <- list(...)
  
  all_cmbs <- combn(names(TAD_list), 2)
  
  if(metric != "getMoc") {
    if(metric_args["matchFor"] == "set1" | metric_args["matchFor"] == "set2")  {
      symmMetric <- FALSE
    }
  }
  if(metric == "get_variationInformation")
    symmMetric <- FALSE
  
  if(!symmMetric)
    all_cmbs <- cbind(all_cmbs, all_cmbs[c(2,1),])  
  
  if(metric == "get_variationInformation") {  # for VI, comparison with itselfs does not yield 1
    selfcmbs <- matrix( c(names(TAD_list), names(TAD_list)), nrow=2, byrow=TRUE)
    all_cmbs <- cbind(all_cmbs, selfcmbs)
  }
  dt <- foreach(i=1:ncol(all_cmbs), .combine = 'rbind') %dopar% {
    set1_name <- paste0(all_cmbs[1,i])
    set2_name <- paste0(all_cmbs[2,i])
    set1_dt <- TAD_list[[paste0(set1_name)]]
    set2_dt <- TAD_list[[paste0(set2_name)]]
    cmpValue <- do.call(metric, list(set1_dt,set2_dt, ...))
    data.frame(
      set1 = set1_name,
      set2 = set2_name,
      cmpValue = cmpValue,
      stringsAsFactors = FALSE
    )
  }
  dt_save <- dt
  if(metric != "get_variationInformation") {
    selfdt <- data.frame(set1=names(TAD_list), set2=names(TAD_list), cmpValue=1, stringsAsFactors = FALSE)
    dt <- rbind(dt, selfdt)
  }

  dt$cmpValue_rd <- round(dt$cmpValue, 2)
  
  tmp <- aggregate(cmpValue ~ set1, FUN=sum, data = dt)
  tmp <- tmp[order(tmp$cmpValue, decreasing = TRUE),]
  set_levels <- tmp$set1
  
  dt$set1 <- factor(dt$set1, levels=(set_levels))
  dt$set2 <- factor(dt$set2, levels=rev(set_levels))
  stopifnot(!is.na(dt$set1))
  stopifnot(!is.na(dt$set2))

  dt_save$set1 <- factor(dt_save$set1, levels=(set_levels))
  dt_save$set2 <- factor(dt_save$set2, levels=(set_levels))
  dt_save <- dt_save[order(as.numeric(dt_save$set1), as.numeric(dt_save$set2)),]

  
  subTit <- metric_args[["matchFor"]]
  if(is.null(subTit)) {
    subTit <- ""
  }
  plotTit <- gsub("get_","", metric)
  
  heatplot <- ggplot(data = dt, aes(x=set1, y=set2, fill=cmpValue)) + 
    ggtitle(paste0(plotTit), subtitle=paste0(subTit))+
    geom_tile()+
    geom_text(aes(label=cmpValue_rd),color = "black", size = 6, fontface="bold") +
    labs(fill=paste0(plotTit))+
    scale_fill_gradient(low=lowColor, high=highColor) +
    theme(
      plot.title = element_text(size=16, hjust=0.5, face = "bold"),
      plot.subtitle = element_text(size=14, hjust=0.5, face = "italic"),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size=12)) 
  if(symmMetric) {
    heatplot <- heatplot +
      theme(
        legend.justification = c(1,0),
        legend.position = c(0.9,0.9),
        legend.direction = "horizontal"
      ) + 
      guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                   title.position = "top", title.hjust = 0.5))
  } else {
    heatplot <- heatplot + 
      guides(fill = guide_colorbar(barwidth = 1, barheight = 7,
                                   title.position = "top", title.hjust = 0.5))
  }
  # heatplot
  return(list(heatplot, dt_save))
}
