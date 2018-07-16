# library(Urine)
#data = filterUS(targetDIR = 'C:\\Users\\songlan\\Desktop\\凤凰中心\\实践\\2018.7.9-7.15\\getTopProteins\\VIP_Firmiana',type = 'ibaq')
#setwd('C:\\Users\\songlan\\Desktop\\凤凰中心\\实践\\2018.7.9-7.15\\getTopProteins')
getTopProteins <- function(data,gps.cutoff = 700){
  data=data.frame(data,row.names = 1)
  data[is.na(data)]=0
  protein_count=data.frame(colnames(data),protein_count=apply(data,2,function(x){return(length(x[x>0]))}))
  sample=row.names(protein_count[which(protein_count$protein_count>gps.cutoff),])
  sample=data[,which(colnames(data) %in% sample)]
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
  return(getTopProteins)
}




