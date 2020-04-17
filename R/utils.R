
load.data <- function(panel){
  datafile = paste("../data/", panel, ".rds", sep = "")
  parentfile = paste("../data/", panel, "_pcounts.rds", sep = "")

  #load metadata
  metadata <- readRDS("../data/metadata.rds")

  #load data
  counts <- readRDS(datafile)

  #load parent population cell counts
  parent_counts <- readRDS(parentfile)

  return(list(counts=counts, parent_counts=parent_counts, metadata=metadata))
}

process.data <- function(data, variables, markers){
  metadata=data$metadata
  counts=data$counts

  #scale Age for model
  metadata$Age <- metadata$Age/100

  #remove unwanted samples
  metadata = metadata[(is.na(metadata$Person)==FALSE),]
  metadata = metadata[(metadata$Person=='P' & metadata$Age<3.96),]
  ids = metadata[metadata$ID %in% rownames(counts),]$ID %>% as.data.frame()
  colnames(ids) = c("ID")

  metadata = left_join(ids, metadata, by="ID")
  counts = counts[match(ids$ID, rownames(counts)),]

  #remove markers that are unwanted
  counts = dplyr::select(counts, -markers)

  #remove metadata columns not in variables
  extra_variables = c('ID', 'PIN', 'Date', 'Age_group', 'Age_group2')
  all_variables = do.call(c, list(extra_variables, variables))
  metadata = dplyr::select(metadata, all_variables)

  return(list(counts=counts, parent_counts=data$parent_counts, metadata=metadata))
}

remove.markers <- function(panel, data){
  if(panel=='bcell'){
    markers = colnames(data[[1]])[grepl('CD1c', colnames(data[[1]]))]}
  else if(panel=='tcell'){
    markers = c('CD3+', 'CD3-circle', 'CD3+Tetramer', 'CD3+Tetramer+CD8-CD4-/CD3+Tetramer+CD25+', 'CD3+Tetramer+CD8-CD4-/CD3+Tetramer+PD1+',
                'CD3+Tetramer+CD8+CD4-/CD3+Tetramer+CD25+', 'CD3+Tetramer+CD8-CD4+/CD3+Tetramer+PD1+', 'CD3+Tetramer+CD8-CD4+/CD3+Tetramer+CD25+',
                'CD3+Tetramer+CD8-CD4-/CD3+Tetramer+PD1+', 'CD3+Tetramer+CD8+CD4-/CD3+Tetramer+PD1+')
  }
  else{
    markers = c('CD3+', 'Total NK', 'CD57+NKG2Chi', 'CD56hiCD16+/NKG2C+ (Total NK)/CD107a-',	'CD56hiCD16+/NKG2C+ (Total NK)/CD107a+', 'CD56hiCD16+/NKG2C+ (Total NK)/CD57-', 'CD56hiCD16+/NKG2C+ (Total NK)/CD57+', 'CD56hiCD16-/NKG2C+ (Total NK)/CD107a-', 'CD56hiCD16-/NKG2C+ (Total NK)/CD107a+')
  }

  return(markers)
}

run.analysis <- function(panel, variables){
  print("Running analysis")
  print(panel)
  #load data
  data = load.data(panel)

  #select which markers to remove:
  markers = remove.markers(panel, data)
  print("Markers to remove:")
  print(markers)

  #prepare data; keep only primaries
  data = process.data(data, variables, markers)

  #run analysis for every marker using variables
  results = perform.fit(data, variables)

  #create plots
  fdr=0.05
  plots = create.plots(data, results$adjp, fdr)
  results = c(results, plots=list(plots))

  results = c(results, panel=panel)
  return(results)
}

save.results <- function(results, directory){
 #create directory for results
  filepath <- file.path(paste("../results/", results$panel, sep=""), directory)
  dir.create(filepath)

  #save p-values to excel file
  pvals_out = cbind(results$pvals, results$adjp, results$coeff)
  write.file(pvals_out, results$formula, filepath)

  #save plots
  save.plots(results$plots, filepath)
}

write.file <- function(pvals, formula, filepath){
  out_file = paste(filepath, "1pvalues.csv", sep="/")
  line1 = paste("Formula: ", formula, "\n", "\n")
  cat(line1, file=out_file)
  write.table(pvals, out_file, sep=",", append=TRUE, row.names=TRUE, col.names=NA)
}


