#get top 1000
getTopProteins <- function(data,gps.cutoff = 1000){
  data[is.na(data)]=0
  getTopProteins=data.frame(GeneSymbol=character(0))
  for(i in 1:ncol(data)){
    tmp=data.frame(GeneSymbol=row.names(data),iBaq=data[,i])
    TopProteins=tmp[order(tmp$iBaq,decreasing=T),][c(1:gps.cutoff),]
    TopProteins$iFOT=(TopProteins$iBaq/sum(TopProteins$iBaq))*100000
    TopProteins=data.frame(GeneSymbol=TopProteins$GeneSymbol,iFOT=TopProteins$iFOT)
    getTopProteins=merge(getTopProteins,TopProteins,by='GeneSymbol',all=T)
  }
  getTopProteins[is.na(getTopProteins)]=0
  colnames(getTopProteins)=c('GeneSymbol',colnames(data))
  getTopProteins=data.frame(getTopProteins,row.names = 1)
  return(getTopProteins)
}