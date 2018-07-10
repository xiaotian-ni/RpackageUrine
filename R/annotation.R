# protein.data = read.table('~/Desktop/learnFunction/PHU00001_protein_data.txt',header = T,row.names = 1)
# rownames(protein.data) = protein.data$GeneSymbol
# protein.data = protein.data[,-1]

annotate <- function(protein.data,gene.description = T,tissue.specificity = T){
  protein.data = data.frame(GeneSymbol = rownames(protein.data),protein.data)
  data("annotationData")
  annotated.data = data.frame(GeneSymbol = anno$GeneSymbol)
  if(gene.description == T){
    annotated.data = data.frame(annotated.data,Gene.description = anno$Gene.description,Protein.class = anno$Protein.class)
  }
  if(tissue.specificity == T){
    annotated.data = data.frame(annotated.data,
                      RNA.tissue.category = anno$RNA.tissue.category,
                      RNA.TS = anno$RNA.TS,
                      RNA.TS.TPM = anno$RNA.TS.TPM
                      )
  }
  annotated.data = merge(annotated.data,protein.data,by = 'GeneSymbol',all.y = T)
  return(annotated.data)
}