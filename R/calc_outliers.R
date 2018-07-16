calc_outliers <- function(protein.data,clinical.info,ri,imputed_value = 0,ri.freq.cutoff = 0.25,col.annos = c(),sort.by = 'none'){
  # ==== convert to lower letter ====
  colnames(clinical.info) = tolower(colnames(clinical.info))
  clinical.info$firmiana_id = tolower(clinical.info$firmiana_id)
  
  colnames(protein.data) = tolower(colnames(protein.data))
  # ==== sort clinical.info ====
  if(sort.by == 'none'){
    cat('Sort by "default".')
  }else{
    clinical.info = clinical.info[order(clinical.info[,which(colnames(clinical.info) == sort.by)]),]
  }
  
  expnames.info = clinical.info$firmiana_id
  cat('Clinical info:',length(expnames.info),'experiments.\n')
  expnames.protein = colnames(protein.data)
  cat('Protein data:',length(expnames.protein),'experiments.\n')
  expnames.intersected = intersect(expnames.info,expnames.protein)
  
  protein.data = protein.data[,colnames(protein.data) %in% expnames.intersected]
  clinical.info = clinical.info[clinical.info$firmiana_id %in% expnames.intersected,]
  cat('Exp kept:',length(expnames.intersected),'.\n')
  
  # ==== sort data by clinical.info ====
  
  protein.data = protein.data[,sort(colnames(protein.data))[rank(clinical.info$firmiana_id)]]
  
  if(length(clinical.info$gender[!duplicated(clinical.info$gender)]) > 1){
    cat('Error: There are two genders.')
    return(0)
  }
  gender = as.character(clinical.info$gender[1])
  cat('Gender:',gender,'\n')
  if(gender == '男'|gender == 'm'){
    ri_4_calc = data.frame(genesymbol = rownames(ri),ri = ri$ri.male)
  }else if(gender == '女'|gender == 'f'){
    ri_4_calc = data.frame(genesymbol = rownames(ri),ri = ri$ri.female)
  }else{
    cat('Error: Gender error')
    return(0)
  }
  p_data = data.frame(genesymbol = rownames(protein.data),protein.data)
  merged_data = merge(ri_4_calc,p_data,by = 'genesymbol',all.y = T)
  merged_data[is.na(merged_data)] = imputed_value

  # ==== calc diff matrix ====
  diff_matrix = merged_data[,-1:-2]
  rownames(diff_matrix) = merged_data$genesymbol
  diff_matrix = diff_matrix - merged_data$ri
  diff_matrix[diff_matrix < 0] = 0
  diff_matrix[diff_matrix > 0] = 1

  outliers_of_each_sample = data.frame(Sample = colnames(diff_matrix),ri.count = apply(diff_matrix,2,sum))

  outlier_freq = data.frame(genesymbol = rownames(diff_matrix),ri.freq = apply(diff_matrix,1,sum))
  outlier_freq_table = merge(outlier_freq,merged_data,by = 'genesymbol',all = T)
  outlier_freq_table$genesymbol = as.character(outlier_freq_table$genesymbol)
  result_info = outlier_freq_table

  x = outliers_of_each_sample$ri.count
  x = c('outlier.count','','',x)
  result_info = rbind(x,result_info)
  
  proteins_detected = data.frame(Sample = colnames(diff_matrix),proteins.detected = apply(protein.data,2,function(x){return(length(x[x>imputed_value]))}))
  x = proteins_detected$proteins.detected
  x = c('proteins.detected','','',x)
  result_info = rbind(x,result_info)
  
  if(length(col.annos) == 0){
    cat('No column annotation detected.')
  }else{
    for(col.anno in col.annos){
      x = as.vector(clinical.info[,which(colnames(clinical.info) == col.anno)])
      x = c(col.anno,'','',x)
      result_info = rbind(x,result_info)
    }
  }

  result = list(diff_matrix = diff_matrix,
                outliers_of_each_sample = outliers_of_each_sample,
                outlier_freq_table = outlier_freq_table,
                result_info = result_info)
  return(result)
}

