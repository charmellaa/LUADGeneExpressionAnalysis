getEnrichmentPlot <- function(annotation,type,top_term,thr_pval,dirEnrich){
  
  annotation <- annotation[annotation$Adjusted.P.value < thr_pval, ]
  annotation <- annotation[order(annotation$Adjusted.P.value, decreasing = F),]
  
  annotation$Gene_count <- sapply(annotation$Genes, function(x){
    
    tmp <- unlist(strsplit(x, split = ";"))
    count <- length(tmp)
    
  })
  
  annotation$Gene_ratio <- unlist(lapply(annotation$Overlap, function(x){
    
    total <- as.numeric(strsplit(x,"/")[[1]][2])
    count <- as.numeric(strsplit(x,"/")[[1]][1])
    
    Gene_ratio <- count/total
    
  } ))
  
  if(length(top_term) != 0 & top_term <= nrow(annotation)){
    annotation_top <- annotation[1:top_term,]
  }else{
    annotation_top <- annotation
  }
  
  # bar plot
  g1 <- ggplot(annotation_top, 
              aes_string(x = "Gene_count", y = fct_reorder(annotation_top$Term, annotation_top$Gene_count), fill = "Adjusted.P.value" )) +
    geom_bar(stat="identity") +
    scale_fill_continuous(low="red", high="blue", name = "Adjusted.P.value",
                          guide=guide_colorbar(reverse=TRUE)) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
    theme_bw(base_size = 10) +
    ylab(NULL)
  
  pdf(paste0(dirEnrich,type,"_barplot.pdf"))
  print(g1)
  dev.off()
  
  # dot plot
  g2 <- ggplot(annotation_top, 
               aes_string(x = "Gene_count", y = fct_reorder(annotation_top$Term, annotation_top$Gene_count))) +
    geom_point(aes(size = Gene_ratio, color = Adjusted.P.value )) +
    scale_colour_gradient(limits=c(min(annotation_top$Adjusted.P.value), max(annotation_top$Adjusted.P.value)), low = "red", high = "blue") +
    theme_bw(base_size = 10) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
    ylab(NULL)
  
  pdf(paste0(dirEnrich,type,"_dotplot.pdf"))
  print(g2)
  dev.off()
  
  write.table(annotation[,c("Term","Overlap","P.value","Adjusted.P.value","Gene_count","Gene_ratio","Genes")],
              paste0(dirEnrich,type,"_adj_pval_",0.01,".txt"), sep = "\t", quote = F, col.names = T, row.names = F)
  
  
}