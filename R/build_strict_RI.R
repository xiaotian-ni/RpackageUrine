build_strict_ri <- function(male,female,upper.ri = 0.975){
  male_ri=c()
  for(i in 1:nrow(male)){
    ri=quantile(male[i,][which(as.numeric(male[i,])>0)],upper.ri,na.rm = T)
    ri=data.frame(genesymbol=row.names(male[i,]),ri)
    male_ri=rbind(male_ri,ri)
  }
  female_ri=c()
  for(i in 1:nrow(female)){
    ri=quantile(female[i,][which(as.numeric(female[i,])>0)],upper.ri,na.rm = T)
    ri=data.frame(genesymbol=row.names(female[i,]),ri)
    female_ri=rbind(female_ri,ri)
  }
  ri = merge(male_ri,female_ri,by = 'genesymbol',all = T)
  colnames(ri) = c('genesymbol','ri.male','ri.female')
  rownames(ri) = ri$genesymbol
  ri = ri[,-1]
  ri = round(ri,2)
  return(ri)
}
