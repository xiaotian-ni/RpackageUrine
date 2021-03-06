annotate <- function(protein.data,gene.description = T,tissue.specificity = T){
  protein.data = data.frame(genesymbol = rownames(protein.data),protein.data)
  data("annotationData")
  data("type_annotation")
  colnames(type_annotation)[1]='genesymbol'
  type_annotation=data.frame(genesymbol=type_annotation$genesymbol,type=type_annotation$Type.1)
  
  annotated.data = data.frame(genesymbol = anno$GeneSymbol)
  
  if(tissue.specificity == T){
    annotated.data = data.frame(annotated.data,
                                #rna.tissue.category = anno$category,
                                number.of.detected.tissues = anno$no..detected.tissues,
                                rna.TS = anno$elevated.tissue.s.,
                                rna.tissue.specific.score = anno$tissue.specific.score,
                                rna.group.specific.score = anno$group.specific.score,
                                rna.enhanced.score = anno$enhanced.score
    )
  }
  if(gene.description == T){
    annotated.data = data.frame(annotated.data,
                                gene.description = anno$description,
                                gene.summary = anno$ncbi_gene_summary)
  }
  annotated.data = merge(annotated.data,type_annotation,by='genesymbol',all.x=T)
  annotated.data = merge(annotated.data,protein.data,by = 'genesymbol',all.y = T,sort = F)
  x = annotated.data$genesymbol
  annotated.data = annotated.data[!duplicated(x),]
  rownames(annotated.data) = annotated.data$genesymbol
  annotated.data = annotated.data[sort(rownames(annotated.data))[rank(rownames(protein.data))],]
  return(annotated.data)
}


