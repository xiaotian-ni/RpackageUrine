# load('Normal_male_Top1000.RData')

build_strict_ri <- function(male,female,upper.ri = 0.975){
  male_ri=c()
  for(i in 1:nrow(male)){
    ifot=as.numeric(male[i,])
    count=apply(male[i,],1,function(x){return(length(x[x>0]))})
    if(count<10){
      ri=0
    }else if(count>=10 & count<120){
      if(upper.ri<0.975){
        ri=quantile(ifot[which(ifot>0)],upper.ri)
        }
      else{
        ri=max(male[i,])
        }
    }else{
      ri=quantile(ifot[which(ifot>0)],upper.ri)
    }
    ri=data.frame(genesymbol=row.names(male[i,]),ri)
    male_ri=rbind(male_ri,ri)
  }
  
  female_ri=c()
  for(i in 1:nrow(female)){
    ifot=as.numeric(female[i,])
    count=apply(female[i,],1,function(x){return(length(x[x>0]))})
    if(count<10){
      ri=0
    }else if(count>=10 & count<120){
      if(upper.ri<0.975){
        ri=quantile(ifot[which(ifot>0)],upper.ri)
      }
      else{
        ri=max(female[i,])
      }
    }else{
      ri=quantile(ifot[which(ifot>0)],upper.ri)
    }    
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













