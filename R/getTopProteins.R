getTopProteins <- function(data,gps.cutoff = 700){
  data[is.na(data)]=0
  protein_count=data.frame(colnames(data),protein_count=apply(data,2,function(x){return(length(x[x>0]))}))
  sample=row.names(protein_count[which(protein_count$protein_count>gps.cutoff),])
  sample=data[,which(colnames(data) %in% sample)]
<<<<<<< HEAD
  getTopProteins=data.frame(genesymbol=character(0))
  for(i in 1:ncol(sample)){
    #i=3
    tmp=data.frame(genesymbol=row.names(sample),ibaq=data[,i])
    topproteins=tmp[order(tmp$ibaq,decreasing=T),][c(1:gps.cutoff),]
    topproteins$ifot=topproteins$ibaq/sum(topproteins$ibaq)*100000
    topproteins=data.frame(genesymbol=topproteins$genesymbol,ifot=topproteins$ifot)
    getTopProteins=merge(getTopProteins,topproteins,by='genesymbol',all=T)
  }
  getTopProteins[is.na(getTopProteins)]=0
  colnames(getTopProteins)=c('genesymbol',colnames(sample))
=======
  sample=sample[which(apply(sample,1,max)>0),]
  getTopProteins=data.frame(GeneSymbol=character(0))
  for(i in 1:ncol(sample)){
    tmp=data.frame(GeneSymbol=row.names(sample),iBaq=sample[,i])
    TopProteins=tmp[order(tmp$iBaq,decreasing=T),][c(1:gps.cutoff),]
    TopProteins$iFOT=(TopProteins$iBaq/sum(TopProteins$iBaq))*100000
    TopProteins=data.frame(GeneSymbol=TopProteins$GeneSymbol,iFOT=TopProteins$iFOT)
    getTopProteins=merge(getTopProteins,TopProteins,by='GeneSymbol',all=T)
  }
  getTopProteins[is.na(getTopProteins)]=0
  colnames(getTopProteins)=c('GeneSymbol',colnames(sample))
  getTopProteins=data.frame(getTopProteins,row.names = 1)
>>>>>>> 3dd6b7968e1209e69684a23d4938d1981ca3a365
  return(getTopProteins)
}




