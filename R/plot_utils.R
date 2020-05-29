
create.plots <- function(data, pvals, fdr){
  #list to store plots
  plots = list(sum(pvals <= 1))

  for(i in 1:nrow(pvals)){
    marker = rownames(pvals)[i]
    marker_pvals = pvals[i,]

    #select significant results
    sig = colnames(pvals)[marker_pvals<fdr] %>% stringr::str_remove("adjp_")

    for(variable in sig){
      #create plot for every significant variable
      name = paste(variable, marker, sep="_")
      plots[[name]] <- local({
        new_plot <- plot.data(data, marker, variable)
        print(new_plot)
      })
    }
  }
  return(plots)
}

plot.data <- function(data, marker, variable){
  counts=data$counts
  parent_counts=data$parent_counts
  md=data$metadata

  #create
  data_tmp <- data.frame(y = as.numeric(counts[md$ID, marker]), md)
  data_tmp <- data_tmp[!is.na(data_tmp$y),]
  weights <- parent_counts[data_tmp$ID, marker]

  #take out anything with parent population <50 cells
  selected = weights>50
  selected[is.na(selected)] <- FALSE
  data_tmp <- data_tmp[selected,]

  #if variable is age, scale back up by 100
  if(variable=="Age") variable="Age_group2"

  #create boxplot
  plot <- data_tmp %>% ggplot2::ggplot(aes_string(x=as.character(variable), y="y")) +
    geom_boxplot() +
    geom_boxplot(outlier.size = 0) + #turn off outliers on boxplot so they don't overlap with dotplot
    geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.75), dotsize=0.5) +
    ggtitle(marker) +
    ylab("Cell proportions") +
    theme_bw() +
    theme(axis.text.y = element_blank(), #remove ticks on y-axis
          axis.ticks.y = element_blank(), #remove tick labels on y-axis
          axis.title.y = element_blank(), #remove y-axis label
          axis.title.x =  element_blank()) #remove x-axis label

  return(plot)
}

save.plots <- function(plots, filepath){
  sig_names = names(plots)
  for(name in sig_names){
    name = sub("/", "_", name)
    file=paste(filepath, "/", name, ".pdf", sep="")
    ggsave(file, plots[[name]])
  }
}

show.plots <- function(plots, name){
  sig_names = names(plots)
  for(name in sig_names){
    print(plots[[name]])
  }
}
