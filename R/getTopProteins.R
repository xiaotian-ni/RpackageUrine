getTopProteins <- function(data,gps.cutoff = 700){
  data[is.na(data)]=0
  protein_count=data.frame(colnames(data),protein_count=apply(data,2,function(x){return(length(x[x>0]))}))
  sample=row.names(protein_count[which(protein_count$protein_count>gps.cutoff),])
  sample=data[,which(colnames(data) %in% sample)]
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
  return(getTopProteins)
}




