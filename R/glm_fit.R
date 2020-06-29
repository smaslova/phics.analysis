
perform.fit <- function(data, variables){
  counts = data$counts %>% t() %>% data.frame()
  colnames(counts) <- rownames(data$counts)
  rownames(counts) <- colnames(data$counts)

  results = differential_abbundance(counts, data$parent_counts, data$metadata, variables)

  return(results)
}


differential_abbundance <- function(counts, parent_counts, md, variables){

  #create formula
  formula <- as.formula(paste('y', paste(paste(variables, collapse=' + '), '(1|PIN) + (1|ID)', sep=' + '), sep=' ~ '))

  ## Fit the GLMM for each cluster separately
  fit_binomial <- lapply(1:nrow(counts), function(i){
    marker = rownames(counts[i,])
    marker = noquote(toString(marker))
    marker = gsub('/', "", marker, fixed = TRUE)

    data_tmp <- data.frame(y = as.numeric(counts[i, md$ID]), md)
    data_tmp <- data_tmp[!is.na(data_tmp$y),]
    weights <- parent_counts[data_tmp$ID, rownames(counts[i,data_tmp$ID])]

    #take out anything with parent population <50 cells
    selected = weights>50
    selected[is.na(selected)] <- FALSE
    data_tmp <- data_tmp[selected,]
    weights <- weights[selected]


    fit_tmp <- glmer(formula, weights=weights, family=binomial(link = "logit"), data = data_tmp,
                     control=glmerControl(optimizer="bobyqa", nAGQ=1, optCtrl=list(maxfun=2e8)))

    coeff_val = fixef(fit_tmp)

    #coefficients to test
    coeff = names(fixef(fit_tmp))
    K <- c()
    for (j in 2:length(coeff)){
      K = append(K, paste(coeff[j], "0", sep=" = "))
    }

    #extract p-values
    contr_tmp <- glht(fit_tmp, linfct = K)
    summ_tmp <- summary(contr_tmp)
    pval <- summ_tmp$test$pvalues

    #r squared
    print(str(summ_tmp))
    r2 = summ_tmp$r.squared

    results = list(pval, coeff_val, r2)
    return(results)
  })

  results <- fit_binomial

  #get coefficients
  coeff1 <- lapply(results, `[[`, 2)
  coeff <- do.call(rbind, coeff1)[,-1]
  colnames(coeff) <- paste0("coeff_", colnames(coeff))
  rownames(coeff) <- rownames(counts)

  #get rsquared values
  rsquared1 <- lapply(results, `[[`, 3)
  r2 <- do.call(rbind, rsquared1)[,-1]
  colnames(r2) <- paste0("R-squared")
  rownames(r2) <- rownames(counts)

  #get p-values
  pvals <- do.call(rbind, lapply(results, `[[`, 1))
  colnames(pvals) <- paste0("pval_", variables)
  rownames(pvals) <- rownames(counts)

  ## Adjust the p-values
  adjp <- apply(pvals, 2, p.adjust, method = "BH")
  colnames(adjp) <- paste0("adjp_", variables)

  return(list(pvals=pvals, adjp=adjp, coeff=coeff, rsquared=r2, formula = format(formula)))
}





