####################################################
######              PREPARE DATA              ###### 
####################################################

DF <- list()
for (i in 1:length(SNPs_OR_type)){
  SNP_QC_param_OR <- SNPs_OR_type[[i]]
  max <- max(c(max(Annotation[,3]), max(SNP_QC_param_OR$position))) #Define limits
  max<- max(max, virus_size)
  
  #Add gene annotation if it hasn't be done before.
  if (is.null(SNP_QC_param_OR$gene)== TRUE){
    gene_annotation <- read.annotation(Annotation, virus_size)
    SNP_QC_param_OR <- gene.annotation(SNP_QC_param_OR, gene_annotation) #snps ha de ser rownames
    genes <- SNP_QC_param_OR$gene
  }else{
    genes <- SNP_QC_param_OR$gene
  }
  
  if (only.genes == TRUE) {
    #if TRUE: delete SNPs that do not belong to any gnee
    SNP_QC_param_OR2 <- SNP_QC_param_OR[!SNP_QC_param_OR$gene == "NA",]
    SNP_QC_param_OR <- SNP_QC_param_OR[!SNP_QC_param_OR$gene == "",]
  }else if(only.genes == FALSE){
    #add NA postions to Genebar:
    a <- c(1:max)
    B <- NULL
    for (d in 1:nrow(Annotation)){
      b <- c(Annotation[d,2]:Annotation[d,3])
      B <- c(B,b)
    }
    B <- unique(B)
    c <- a[-B]
    sta <- c(B+1,c)
    sta <- c(0,sta[duplicated(sta)])
    fin <- c(B-1,c)
    fin <- c(fin[duplicated(fin)],max)
    add <- as.data.frame(cbind("NA", sta, fin))
    colnames(add) <- colnames(Annotation)
    Annotation <- rbind(Annotation, add)
  }
  
  #Plot points in their position (TRUE) or inside the gene group (FALSE)
  if (position == TRUE){
    #if TRUE: plot each SNP at each position
    x <- as.integer(SNP_QC_param_OR$position)
  }else{
    #if FALSE: plot each SNP at their gene group
    x <- droplevels(SNP_QC_param_OR$gene)
    Genebar <- FALSE
  }
  
  # PLOT points with different size according the OR (TRUE) or with the same size (FALSE)
  if (OR.size == TRUE){
    #if TRUE: change point size according to the risk (major vs minor allele)
    OddsRatio <- as.numeric(as.character(SNP_QC_param_OR$OR))
  }else{
    OddsRatio <- rep(1,nrow(SNP_QC_param_OR))
  }
  y <- -log10(as.numeric(as.character(SNP_QC_param_OR$pval)))
  # print(as.numeric(as.character(SNP_QC_param_OR$pval)) > 1)
  
  # Genes <- as.character(SNP_QC_param_OR$gene)
  
  z <- rep(names(SNPs_OR_type)[i], length(x))
  df1 <- data.frame(x, y, z, OddsRatio, genes)
  DF[[i]] <- df1
}

df_manhattan <- do.call(rbind.data.frame, DF)
df_manhattan <- df_manhattan[!df_manhattan$genes == "NA", ]
SNP_QC_param_OR <- df_manhattan
colnames(SNP_QC_param_OR) <- c("position", "pval", "Lineages", "OR", "gene")


####################################################
######             PLOT MANHATTAN             ###### 
####################################################
suppressMessages(library("viridis"))
suppressMessages(library(cowplot))
suppressMessages(library(ggplot2))

max <- max(c(max(Annotation[,3]), max(SNP_QC_param_OR$position)))
max<- max(max, virus_size)

#Add gene annotation if it hasn't be done before.
if (is.null(SNP_QC_param_OR$gene)== TRUE){
  gene_annotation <- read.annotation(Annotation, virus_size)
  SNP_QC_param_OR <- gene.annotation(SNP_QC_param_OR, gene_annotation) #snps ha de ser rownames
}

if (only.genes == TRUE) {
  #if TRUE: delete SNPs with no gene annotation
  SNP_QC_param_OR2 <- SNP_QC_param_OR[!SNP_QC_param_OR$gene == "NA",]
  SNP_QC_param_OR <- SNP_QC_param_OR[!SNP_QC_param_OR$gene == "",]
}else if(only.genes == FALSE){
  #add NA postions to Genebar:
  a <- c(1:max)
  B <- NULL
  for (d in 1:nrow(Annotation)){
    b <- c(Annotation[d,2]:Annotation[d,3])
    B <- c(B,b)
  }
  B <- unique(B)
  c <- a[-B]
  sta <- c(B+1,c)
  sta <- c(0,sta[duplicated(sta)])
  fin <- c(B-1,c)
  fin <- c(fin[duplicated(fin)],max)
  add <- as.data.frame(cbind("NA", sta, fin))
  colnames(add) <- colnames(Annotation)
  Annotation <- rbind(Annotation, add)
}

if (position == TRUE){
  #if TRUE: plot each SNP at each position
  x <- as.integer(SNP_QC_param_OR$position)
}else{
  #if FALSE: plot each SNP at their gene group
  x <- droplevels(SNP_QC_param_OR$gene)
  Genebar <- FALSE
  
}

if (OR.size == TRUE){
  #if TRUE: change point size according to the risk (major vs minor allele)
  OddsRatio <- as.numeric(as.character(SNP_QC_param_OR$OR))
}else{
  OddsRatio <- rep(1,nrow(SNP_QC_param_OR))
}


#pval
# y <- -log10(as.numeric(as.character(SNP_QC_param_OR$pval)))
y <- SNP_QC_param_OR$pval
Genes <- as.character(SNP_QC_param_OR$gene)
y_OR <- round(max(y))

if (Genebar == TRUE) {
  if (is.data.frame(Annotation)==FALSE){
    stop("Gene annotation data.frame is required.")
  }
    if (Lineages == TRUE) {
      z <- SNP_QC_param_OR$Lineages
      df1 <- data.frame(x,y,z)
      len <- nrow(df1)
    }else{
      #define well genes to have the same color than gene bar
      Genes <- as.character(Genes)
      len <- length(Genes)
      lev <- levels(Annotation[,1])
      Genes <- factor(c(Genes, lev))
      add <- rep(0,length(lev))
      y <- c(y, add)
      OddsRatio <- c(OddsRatio, add)
      x <- c(x, add)
      z <- SNP_QC_param_OR$gene
      df1 <- data.frame(x,y,z)
    }
  
  #Define Genebar parameters: 
  max <- max(Annotation[,3])
  Genes <- Genes[1:len] 
  gene <- NULL; start <- NULL; end <- NULL; df_genebar <- NULL; ccc <- 0
  for (i in 1:nrow(Annotation)){
    g <- as.character(Annotation[i,1]) #gene name
    # w <- Annotation[i,1] #y value
    s <- Annotation[i,2] #start
    e <- Annotation[i,3] #end
    ccc <- ccc+1
    gene <- c(gene, g)
    start <- c(start, s)
    end <- c(end, e)
    vec <- rbind(c(g,s,ccc), c(g,e,ccc))
    df_genebar <- rbind(df_genebar, vec)
  }
  df_genebar <- as.data.frame(df_genebar)
  colnames(df_genebar) <- c("y","x","group")
  
  #plot genomic ranges
  if (is.na(gene_bar_color) == TRUE){
    genebar <- ggplot(data=df_genebar, aes(x=as.integer(as.character(x)), y=y, group=group)) +
      geom_line(aes(color = y, size = 3))+
      geom_point(aes(color = df_genebar$y, size = 0)) +
      scale_color_viridis(discrete = TRUE, option = col) +
      xlab(xlab) +
      xlim(0, max) +
      theme(legend.position="none",
            legend.text = element_text(color = "white", size = 1),
            axis.text.y=element_text(size=6),
            axis.title.y=element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"))
  }else{
    genebar <- ggplot(data=df_genebar, aes(x=as.integer(as.character(x)), y=y, group=group)) +
      geom_line(aes(color = y, size = 3))+
      geom_point(aes(color = y, size = 0)) +
      scale_color_manual(values=c(rep(gene_bar_color, length(levels(df_genebar$y))))) +
      xlab(xlab) +
      xlim(0, max) +
      theme(legend.position="right",
            legend.text = element_text(color = "white", size = 1),
            axis.text.y=element_text(size=8),
            axis.text.x=element_text(size=8),
            axis.title.y=element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"))
  }
  genebar <- genebar + geom_text(x=6700, y=9, label="URR", cex = 3) +
    geom_text(x=3900, y=8, label="L2", cex = 3) +
    geom_text(x=5300, y=7, label="L1", cex = 3) +
    geom_text(x=1200, y=6, label="E7", cex = 3) +
    geom_text(x=900, y=5, label="E6", cex = 3) +
    geom_text(x=4500, y=4, label="E5", cex = 3) +
    geom_text(x=3900, y=3, label="E4", cex = 3) +
    geom_text(x=4100, y=2, label="E2", cex = 3) +
    geom_text(x=3100, y=1, label="E1", cex = 3)
    
    

  man <- ggplot(df1, aes(x=x, y=y, group = z)) +
    geom_point(aes(color = z, size = 4)) +
    scale_color_viridis(discrete = TRUE, option = "B") +
    geom_hline(yintercept=-log10(0.05), col = "black", size = 1) +
    ylab("-log10(Pval)") +
    xlim(0, max) +
    geom_point(aes(x=0, y=0), col="purple", size = 0.01) +
    geom_point(aes(x=(max), y=0), col="yellow", size = 0.01) +
    theme(legend.position= "right",
          axis.text.y=element_text(size=10),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle(main)
  
  man <- plot_grid(man, genebar,ncol = 1, nrow = 2,  rel_heights = c(6/9, 3/9), align = "v")
  
  
  
}else{
  if (Lineages == TRUE) {
    z <- SNP_QC_param_OR$Lineages
    df1 <- data.frame(x,y,z)
  }else{
    #define well genes to have the same color than gene bar
    z <- SNP_QC_param_OR$gene
    df1 <- data.frame(x,y,z)
  }
  if (OR.size == T){
    man <- ggplot(df1, aes(x=x, y=y)) +
      geom_point(aes(size=OddsRatio, color = z)) +
      geom_point(aes(x=max, y=(y_OR +0.5), size=round((3/4)*(max(OddsRatio)))), colour="black") +
      annotate("text", label = paste("OR size:",round((3/4)*(max(OddsRatio))), sep = " "), x = (x_max-200), y = (y_OR +0.5), size = 3) +
      geom_point(aes(x=max, y=y_OR, size=round(mean(OddsRatio))), colour="black") +
      annotate("text", label = as.character(round(mean(OddsRatio))), x = (x_max), y = (y_OR), size = 3) +
      geom_point(aes(x=max, y=(y_OR - 0.5), size=1), colour="black") +
      annotate("text", label = "1", x = (x_max), y = (y_OR -0.5), size = 3) +
      scale_color_viridis(discrete = TRUE, option = col) +
      geom_hline(yintercept=-log10(0.05), color = "black", size = 1) +
      ylab("-log10(Pval)") +
      xlim(0, max) +
      geom_point(aes(x=0, y=0), colour="white") +
      geom_point(aes(x=(max), y=0), colour="white") +
      theme(legend.position="none",
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      ggtitle(main)
  }else{
    man <- ggplot(df1, aes(x=x, y=y)) +
      geom_point(aes(size=OddsRatio, color =  z)) +
      scale_color_viridis(discrete = TRUE, option = col) +
      geom_hline(yintercept=-log10(0.05), color = "black", size = 1) +
      ylab("-log10(Pval)") +
      xlim(0, max) +
      geom_point(aes(x=0, y=0), colour="white") +
      geom_point(aes(x=(max), y=0), colour="white") +
      theme(legend.position="none",
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      ggtitle(main)
  }
  
}
  
man
