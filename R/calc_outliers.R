calc_outliers <- function(protein.data,clinical.info,ri,imputed_value = 0,ri.freq.cutoff = 0.25,col.annos = c(),sort.by = 'none'){
  # ==== sort clinical.info ====
  if(sort.by == 'none'){
    cat('Sort by "default".')
  }else{
    clinical.info = clinical.info[order(clinical.info[,which(colnames(clinical.info) == sort.by)]),]
  }
  
  expnames.info = clinical.info$Firmiana_ID
  cat('Clinical info:',length(expnames.info),'experiments.\n')
  expnames.protein = colnames(protein.data)
  cat('Protein data:',length(expnames.protein),'experiments.\n')
  expnames.intersected = intersect(expnames.info,expnames.protein)
  
  protein.data = protein.data[,colnames(protein.data) %in% expnames.intersected]
  clinical.info = clinical.info[clinical.info$Firmiana_ID %in% expnames.intersected,]
  cat('Exp kept:',length(expnames.intersected),'.\n')
  
  # ==== sort data by clinical.info ====
  
  protein.data = protein.data[,sort(colnames(protein.data))[rank(clinical.info$Firmiana_ID)]]
  
  if(length(clinical.info$Gender[!duplicated(clinical.info$Gender)]) > 1){
    cat('Error: There are two genders.')
    return(0)
  }
  gender = as.character(clinical.info$Gender[1])
  cat('Gender:',gender,'\n')
  if(gender == '男'){
    ri_4_calc = data.frame(GeneSymbol = rownames(ri),RI = ri$ri.male)
  }else if(gender == '女'){
    ri_4_calc = data.frame(GeneSymbol = rownames(ri),RI = ri$ri.female)
  }else{
    cat('Error: Gender error')
    return(0)
  }
  p_data = data.frame(GeneSymbol = rownames(protein.data),protein.data)
  merged_data = merge(ri_4_calc,p_data,by = 'GeneSymbol',all.y = T)
  merged_data[is.na(merged_data)] = imputed_value

  # ==== calc diff matrix ====
  diff_matrix = merged_data[,-1:-2]
  rownames(diff_matrix) = merged_data$GeneSymbol
  diff_matrix = diff_matrix - merged_data$RI
  diff_matrix[diff_matrix < 0] = 0
  diff_matrix[diff_matrix > 0] = 1

  outliers_of_each_sample = data.frame(Sample = colnames(diff_matrix),RI.count = apply(diff_matrix,2,sum))

  outlier_freq = data.frame(GeneSymbol = rownames(diff_matrix),RI.freq = apply(diff_matrix,1,sum))
  outlier_freq_table = merge(outlier_freq,merged_data,by = 'GeneSymbol',all = T)
  outlier_freq_table$GeneSymbol = as.character(outlier_freq_table$GeneSymbol)
  result_info = outlier_freq_table

  x = outliers_of_each_sample$RI.count
  x = c('','','Outlier.count',x)
  result_info = rbind(x,result_info)
  
  proteins_detected = data.frame(Sample = colnames(diff_matrix),proteins.detected = apply(protein.data,2,function(x){return(length(x[x>imputed_value]))}))
  x = proteins_detected$proteins.detected
  x = c('','','Proteins.detected',x)
  result_info = rbind(x,result_info)
  
  if(length(col.annos) == 0){
    cat('No column annotation detected.')
  }else{
    for(col.anno in col.annos){
      x = as.vector(clinical.info[,which(colnames(clinical.info) == col.anno)])
      x = c('','',col.anno,x)
      result_info = rbind(x,result_info)
    }
  }

  result = list(diff_matrix = diff_matrix,
                outliers_of_each_sample = outliers_of_each_sample,
                outlier_freq_table = outlier_freq_table,
                result_info = result_info)
  return(result)
}

