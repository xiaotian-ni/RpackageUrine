# library(Urine)
#data = filterUS(targetDIR = 'C:\\Users\\songlan\\Desktop\\凤凰中心\\实践\\2018.7.9-7.15\\getTopProteins\\VIP_Firmiana',type = 'ibaq')
#setwd('C:\\Users\\songlan\\Desktop\\凤凰中心\\实践\\2018.7.9-7.15\\getTopProteins')
getTopProteins <- function(data,gps.cutoff = 700){
  data=data.frame(data,row.names = 1)
  data[is.na(data)]=0
  protein_count=data.frame(colnames(data),protein_count=apply(data,2,function(x){return(length(x[x>0]))}))
  sample=row.names(protein_count[which(protein_count$protein_count>gps.cutoff),])
  sample=data[,which(colnames(data) %in% sample)]
  getTopProteins=data.frame(GeneSymbol=character(0))
  for(i in 1:ncol(sample)){
    #i=3
    tmp=data.frame(GeneSymbol=row.names(sample),iBaq=data[,i])
    TopProteins=tmp[order(tmp$iBaq,decreasing=T),][c(1:gps.cutoff),]
    TopProteins$iFOT=TopProteins$iBaq/sum(TopProteins$iBaq)*100000
    TopProteins=data.frame(GeneSymbol=TopProteins$GeneSymbol,iFOT=TopProteins$iFOT)
    getTopProteins=merge(getTopProteins,TopProteins,by='GeneSymbol',all=T)
  }
  getTopProteins[is.na(getTopProteins)]=0
  colnames(getTopProteins)=c('GeneSymbol',colnames(sample))
  return(getTopProteins)
}




