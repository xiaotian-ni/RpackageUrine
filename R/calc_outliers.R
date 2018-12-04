calc_outliers <- function(protein.data,clinical.info,ri,sri,imputed_value = 0,ri.freq.cutoff = 0.25,col.annos = c(),sort.by = 'none'){
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
    sri_4_calc = data.frame(genesymbol = rownames(sri),sri = sri$ri.male)
  }else if(gender == '女'|gender == 'f'){
    ri_4_calc = data.frame(genesymbol = rownames(ri),ri = ri$ri.female)
    sri_4_calc = data.frame(genesymbol = rownames(sri),sri = sri$ri.female)
  }else{
    cat('Error: Gender error')
    return(0)
  }
  p_data = data.frame(genesymbol = rownames(protein.data),protein.data)
  
  #SRI
  merged_data = merge(sri_4_calc,p_data,by = 'genesymbol',all.y = T)
  merged_data[is.na(merged_data)] = imputed_value
  # ==== calc diff matrix ====
  diff_matrix = merged_data[,-1:-2]
  rownames(diff_matrix) = merged_data$genesymbol
  SRI_diff_matrix = diff_matrix - merged_data$sri
  SRI_diff_matrix[SRI_diff_matrix < 0] = 0
  SRI_diff_matrix[SRI_diff_matrix > 0] = 1
  SRI_outliers_of_each_sample = data.frame(Sample = colnames(SRI_diff_matrix),sri.count = apply(SRI_diff_matrix,2,sum))
  SRI_outlier_freq = data.frame(genesymbol = rownames(SRI_diff_matrix),sri.freq = apply(SRI_diff_matrix,1,sum))
  merged_data = merge(SRI_outlier_freq,merged_data,by = 'genesymbol',all = T)
  merged_data$genesymbol = as.character(merged_data$genesymbol)
  # result_info = SRI_outlier_freq
  # ==== calc diff matrix ====
  
  #RI
  merged_data = merge(ri_4_calc,merged_data,by = 'genesymbol',all.y = T)
  merged_data[is.na(merged_data)] = imputed_value
  # ==== calc diff matrix ====
  diff_matrix = merged_data[,-1:-4]
  rownames(diff_matrix) = merged_data$genesymbol
  RI_diff_matrix = diff_matrix - merged_data$ri
  RI_diff_matrix[RI_diff_matrix < 0] = 0
  RI_diff_matrix[RI_diff_matrix > 0] = 1
  RI_outliers_of_each_sample = data.frame(Sample = colnames(RI_diff_matrix),ri.count = apply(RI_diff_matrix,2,sum))
  RI_outlier_freq = data.frame(genesymbol = rownames(RI_diff_matrix),ri.freq = apply(RI_diff_matrix,1,sum))
  RI_outlier_freq_table = merge(RI_outlier_freq,merged_data,by = 'genesymbol',all = T)
  RI_outlier_freq_table$genesymbol = as.character(RI_outlier_freq_table$genesymbol)
  result_info = RI_outlier_freq_table
  # ==== calc diff matrix ====
  
  
  # ====  add ROW====
  x=apply(protein.data,2,function(x){return(length(x[x>imputed_value]))})
  y=RI_outliers_of_each_sample$ri.count
  z=SRI_outliers_of_each_sample$sri.count
  proteins.detected = c('proteins.detected','','','','',x)
  RI.outlier.count= c('RI.outlier.count','','','','',y)
  RI.outlier.freq=c('RI.outlier.freq','','','','',round(y/x,3))
  SRI.outlier.count= c('SRI.outlier.count','','','','',z)
  SRI.outlier.freq=c('SRI.outlier.freq','','','','',round(z/x,3))
  
  
  result_info = rbind(SRI.outlier.freq,result_info)
  result_info = rbind(SRI.outlier.count,result_info)
  result_info = rbind(RI.outlier.freq,result_info)
  result_info = rbind(RI.outlier.count,result_info)
  result_info = rbind(proteins.detected,result_info)
  
  if(length(col.annos) == 0){
    cat('No column annotation detected.')
  }else{
    for(col.anno in col.annos){
      x = as.vector(clinical.info[,which(colnames(clinical.info) == col.anno)])
      x = c(col.anno,'','',x)
      result_info = rbind(x,result_info)
    }
  }
  result_info=data.frame(result_info,row.names = 1)
  result = list(diff_matrix = diff_matrix,
                RI_outliers_of_each_sample = RI_outliers_of_each_sample,
                SRI_outliers_of_each_sample = SRI_outliers_of_each_sample,
                result_info = result_info)
  return(result)
}

