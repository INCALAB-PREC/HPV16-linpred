##################################################################
###                           R-code                           ###
###                                                            ###
###                         FUNCTIONS                          ###
###                                                            ###
##################################################################
###                                                            ###
###                   LAURA ASENSIO PUIG                       ###
###                  lasensio@IDIBELL.cat                      ###
###                                                            ###
##################################################################

####   FUNCTION 1: setup.genome.df   ####
#CREATE A DATAFRAME WITH ALL THE FASTA SEQUENCES

setup.genome.df <- function(sequences){
  ID <- names(sequences)
  genome <- NULL
  for(i in 1:length(sequences)){
    genome <- cbind(genome, sequences[[i]])
  }
  colnames(genome) <- ID
  rownames(genome) <- c(1:nrow(genome))
  genome <- t(genome)
  return(genome)
}

## USAGE:
# library(sequinr)
# seq <- seqinr::read.fasta("sequences.fasta", seqtype = "DNA")
# sequences <- getSequence(seq, as.string = FALSE)
# GENOME <- setup.genome.df(sequences)


# ************************************************************* #
####   FUNCTION 2: findSNPs   ####
#FIND SNPs candidates (mutations)

findSNPs <- function(GENOME, info_seq = NULL) {
  #info_seq vector must follow the same order: ID, lineage, casecon
  SNP.candidates <- function(tbl_nt){
    SNP_candidates <- NULL
    for (i in 1:length(tbl_nt)){
      if(length(tbl_nt[[i]])==1){ 
      }else if(length(tbl_nt[[i]])==2){
        if(length(grep("n", names(tbl_nt[[i]])))==0){
          SNP_candidates <- c(SNP_candidates,tbl_nt[i])
        }
      }else if(length(tbl_nt[[i]])>=3){ 
        if(length(grep("n", names(tbl_nt[[i]])))==1){
          pos <- grep("n",names(tbl_nt[[i]]))
          N <- tbl_nt[[i]][-pos]
          N <- list(N)
          names(N) <- names(tbl_nt[i])
          SNP_candidates <- c(SNP_candidates, N)
        }else{
          SNP_candidates <- c(SNP_candidates,tbl_nt[i])
        }
      }
    }
    SNP_candidates[sapply(SNP_candidates, is.null)] <- NULL
    return(SNP_candidates)
  }
  
  if (is.null(info_seq) ==TRUE){
    tbl_nt <- apply(GENOME, 2, table) 
    names(tbl_nt) <- c(paste("snp",colnames(GENOME),sep=""))
    SNPs_candidates <- SNP.candidates(tbl_nt)
  }else{
    LINEAGE <- levels(droplevels(info_seq[,1]))
    SNPs_candidates <- NULL
    for (i in 1:length(LINEAGE)){
      info_lineage <- grep(LINEAGE[i], info_seq[,1])
      genome_lin <- GENOME[info_lineage,]
      tbl_nt <- apply(genome_lin, 2, table)
      names(tbl_nt) <- c(paste("snp",1:length(tbl_nt),sep=""))
      SNPs_cand <- SNP.candidates(tbl_nt)
      SNPs_candidates[[i]] <- SNPs_cand
    }
    names(SNPs_candidates) <- LINEAGE
  }
  return(SNPs_candidates)
}

## USAGE:
# SNP_candidates <- findSNPs(GENOME)




# ************************************************************* #
####   FUNCTION 3: MVF   ####
#Calculate MVF (minim variant frequency) and call rate. 

MVF <- function(SNP_candidates){
  mvf.calcul <- function(table){
    if(length(table)==2){
      A <- table[1]
      B <- table[2]
      if(A == B){
        mvf <- 0.5
        max_a <- names(table[1])
        minor_a <- names(table[2])
        count_max <-as.numeric(A)
        count_min <-as.numeric(B)
        row <- c(max_a, count_max, minor_a, count_min, mvf)
      }else{
        mvf <- min(c(A,B))/(max(c(A,B)+min(c(A,B))))
        minor_a <- rownames(table)[grep(paste("^",min(c(A,B)),"$",sep = ""),table)]
        count_min <- as.numeric(table[grep(paste("^",min(c(A,B)),"$",sep = ""),table)])
        max_a <- rownames(table)[grep(paste("^",max(c(A,B)),"$",sep = ""),table)]
        count_max <- as.numeric(table[grep(paste("^",max(c(A,B)),"$",sep = ""),table)])
        row <- c(max_a, count_max, minor_a, count_min, mvf)
      }
    }else if(length(table)==4){
      a <- table[1]
      b <- table[2]
      c <- table[3]
      d <- table[4]
      order <- table[order(c(a,b,c,d))]
      min <- as.numeric(order[1])
      max <- as.numeric(order[4])
      mvf <- min/(max+min)
      minor_a <- rownames(table)[grep(paste("^",min,"$",sep = ""),table)]
      minor_a <- paste(minor_a,collapse="/")
      count_min <- as.numeric(table[grep(paste("^",min,"$",sep = ""),table)])
      count_min <- max(count_min)
      max_a <- rownames(table)[grep(paste("^",max,"$",sep = ""),table)]
      max_a <- paste(max_a,collapse="/")
      count_max <- as.numeric(table[grep(paste("^",max,"$",sep = ""),table)])
      count_max<- max(count_max)
      row <- c(max_a, count_max, minor_a, count_min, mvf)
    }else{
      if((length(grep("-",names(table))))==0){
        a <- table[1]
        b <- table[2]
        c <- table[3]
        order <- table[order(c(a,b,c))]
        min <- as.numeric(order[2])
        max <- as.numeric(order[3])
        mvf <- min/(max+min)
        minor_a <- rownames(table)[grep(paste("^",min,"$",sep = ""),table)]
        count_min <- as.numeric(table[grep(paste("^",min,"$",sep = ""),table)])
        if (length(minor_a)==2){
          minor_a <- paste(minor_a[1],minor_a[2], sep= "/")
          count_min <- count_min[1]+count_min[2]
        }
        max_a <- rownames(table)[grep(paste("^",max,"$",sep = ""),table)]
        count_max <- as.numeric(table[grep(paste("^",max,"$",sep = ""),table)])
        if (length(max_a)==2){
          max_a <- paste(max_a[1],max_a[2], sep= "/")
          count_max <- count_max[1]+count_max[2]
        }
        delete <- as.numeric(order[1])
        delete_a <- rownames(table)[grep(delete,table)]
        row <- c(max_a, count_max, minor_a, count_min, mvf)
      }else{
        table <- table[!names(table)=="-"]
        a <- table[1]
        b <- table[2]
        mvf <- min(c(a,b))/(max(c(a,b)+min(c(a,b))))
        minor_a <- rownames(table)[grep(paste("^",min(c(a,b)),"$",sep = ""),table)]
        count_min <- as.numeric(table[grep(paste("^",min(c(a,b)),"$",sep = ""),table)])
        max_a <- rownames(table)[grep(paste("^",max(c(a,b)),"$",sep = ""),table)]
        count_max <- as.numeric(table[grep(paste("^",max(c(a,b)),"$",sep = ""),table)])
        row <- c(max_a, count_max, minor_a, count_min, mvf)
      }
    }
    return(row)
  }
  mvf <- lapply(SNP_candidates, mvf.calcul)
  SNP_name <- names(mvf)
  SNP_pos <- as.integer(substr(SNP_name,4,nchar(SNP_name)))
  mvf_df <- NULL
  for(i in 1:length(mvf)){
    row <- t(mvf[[i]])
    if (length(row) == 5){
      mvf_df <- rbind(mvf_df, row)
    }else{
      warning(paste("It may be a error with this position", names(mvf[i]), sep = " "))
      show(i)
      show(mvf[i])
      problem <- c("X","None","All","0","0")
      mvf_df <- rbind(mvf_df, problem)
    }
  }
  mvf_df <- as.data.frame(mvf_df)
  rownames(mvf_df) <- names(mvf)
  SNP_ID <- paste(SNP_name,toupper(mvf_df[,1]), sep = "")
  mvf_df <- data.frame(SNP_pos, mvf_df)
  rownames(mvf_df) <- SNP_ID
  colnames(mvf_df) <- c("position","major", "count_max", "minor", "count_min", "mvf")
  #CALL
  mvf_df <- transform(mvf_df, count_max = as.numeric(as.character(count_max)),
                      count_min = as.numeric(as.character(count_min)),
                      mvf = as.numeric(as.character(mvf)))
  real_call <- lapply(SNP_candidates, sum)#real
  mvf_df$call <- unlist(real_call)
  N <- max(mvf_df$call)
  mvf_df$call.rate <- floor(mvf_df$call/N*100)
  return(mvf_df)
}

## USAGE
# SNP_parameters <- MVF(SNP_candidates)
## MORE: Add gene annotation for each position: (function )
# SNP_parameters <- annotation(SNP_parameters, virus= "HPV16")




# ************************************************************* #
####   FUNCTION 4: filter.snps   ####
#Filter SNPs according to MVF and callrate parameters

filter.snps <- function(SNP_parameters, mvf = 0.05, callrate = 95, show = TRUE){
  SNP_i <- nrow(SNP_parameters)
  QC_call.rate <- sum(!SNP_parameters$call.rate>callrate)
  SNP_parameters <- SNP_parameters[SNP_parameters$call.rate>callrate,]
  QC_mvf <- sum(!SNP_parameters$mvf>mvf)
  SNP_parameters <- SNP_parameters[SNP_parameters$mvf>mvf,]
  QC_NA <- sum(SNP_parameters$minor=="-")
  SNP_parameters <- SNP_parameters[!SNP_parameters$minor=="-",]
  SNP_f <- SNP_i - nrow(SNP_parameters)
  percentage <-  as.numeric(format(round(SNP_f/SNP_i*100, 2), nsmall = 2))
  SNPs_deleted_QC <- QC_call.rate + QC_mvf + QC_NA
  deleted_SNPs <- rbind(QC_call.rate,QC_mvf,QC_NA, SNPs_deleted_QC)
  deleted_SNPs <- as.data.frame(deleted_SNPs)
  rownames(deleted_SNPs) <- c("Call rate", "MVF", "NA", "TOTAL:")
  colnames(deleted_SNPs) <- "SNPs deleted"
  
  message0 <- paste("Parameters used are MVF = ", mvf, " and callrate = ", callrate, "%.", sep = "")
  message <- paste(SNP_f, " SNPs deleted out of ", SNP_i, " Candidates -  ", percentage, "% of deleted SNPs.",sep= "")
  if (show == TRUE){
    show(message0)
    show(message)
    show(deleted_SNPs)
  }
  return(SNP_parameters)
}
## USAGE
# filtered_SNP_parameters<- filter.snps(SNP_parameters)


###########################
####   THE FUNCTION    ####
###########################
# Summary of all the previous functions (2,3,4) to find the SNPs in one single step. 

SNPs <- function(GENOME, wd_output, mvf, callrate){
  #count nt
  x0 <- ncol(GENOME)
  
  #find SNPs
  SNP_candidates <- findSNPs(GENOME)
  x1 <- length(SNP_candidates)
  
  #Calculate MVF and callrate
  SNP_parameters <- MVF(SNP_candidates)
  
  #FILTER by MVF and callrate
  SNP_QC_param <- filter.snps(SNP_parameters, mvf, callrate)
  x2 <- nrow(SNP_QC_param)
  
  #Take each SNP position nucleotide:
  sbst <- as.integer(substr(rownames(SNP_QC_param),4,nchar(rownames(SNP_QC_param))-1))
  SNP_geno <- GENOME[,(colnames(GENOME) %in% sbst)]
  
  colnames(SNP_geno) <- rownames(SNP_QC_param)
  SNP_geno <- as.data.frame(SNP_geno)
  dim(SNP_geno)
  # write.table(SNP_geno, paste(wd_output, "SNP_geno_HC.txt", sep = ""))
  return(SNP_QC_param)
}

#########################################
####         GENE ANNOTATION         ####
#########################################
#Use the following functions to know which gene belongs each position in HPV16 (aligned to reference genome)

#### FUNCTION 5: read.annotation
#From a gene Annotation file (3 columns: gene, start position, end position), create a new dataframe with the information for each position

read.annotation <- function(Annotation, virus_size){
  # gene_position <- read.table(filepath)
  df <- data.frame(matrix(NA, nrow = virus_size, ncol =nrow(Annotation)))
  rownames(df)<- c(1:virus_size)
  colnames(df)<- Annotation[,1]
  genes <- as.character(Annotation[,1])
  for(i in 1:nrow(Annotation)){
    vec <- c(Annotation[i,2]:Annotation[i,3])
    df[vec,i]<- genes[i]
  }
  j <- 2
  new <- df[,1]
  while(j < length(genes)+1){
    new <- paste(new,df[,j],sep="/")
    j <- j +1
  }
  new <- gsub("NA/","", new)
  new <- gsub(" ","NA",new)
  new <- gsub("/NA","", new)
  
  gene_annotation <- data.frame(c(1:virus_size),new)
  colnames(gene_annotation) <- c("position", "gene")
  return(gene_annotation)
}

# # USAGE:
# annotation <- read.table("/path/Gene_annotation.txt")
# gene_annotation <- read.annotation(annotation, virus_size = 5009)


#### FUNCTION 6: annotation
# From a position vector or column, creates a new vector or column with the gene name for each position. 
gene.annotation <- function(data, gene_annotation, Ncol = 0){
  if ((is.vector(data)==TRUE | is.character(data)==TRUE)==TRUE){
    if(substr(data,1,3)[1]== "snp"){
      logical <- gene_annotation[,1] %in% substr(data,4,nchar(data)-1)
    }else{
      logical <- gene_annotation[,1] %in% data
    }
    position <- subset(gene_annotation, logical)
    rownames(position) <- data
    data <- position
  }else if ((is.data.frame(data)==TRUE | is.matrix(data)==TRUE)==TRUE){
    if (Ncol == 0){
      logical <- as.character(gene_annotation[,1]) %in% substr(rownames(data),4,nchar(rownames(data))-1)
      position <- subset(gene_annotation, logical)
      data$SNPs <- rownames(data)
      data <- merge(data, position, by = "position")
      rownames(data) <-data$SNPs 
      data <- subset(data, select = -SNPs)
    }else{
      logical <- as.character(gene_annotation[,1]) %in% substr(as.character(data[,Ncol]),4,nchar(as.character(data[,Ncol]))-1)
      position <- subset(gene_annotation, logical)
      data$SNPs <- rownames(data)
      data <- merge(data, position, by = "position")
      rownames(data) <-data$SNPs 
      data <- subset(data, select = -SNPs)
    }
  }else {
    stop("Error: vector, matrix or data.frame is required.")
  }
  return(data)
}


pvalue.return <- function(pval){
  if( pval < 0.001){pval <- "<0.001"}else{pval <- round(pval,3)}
  return(pval)}



##################################################
### FUNCTION 7.1 (internal function)
#Internal function for function 7: Single.SNP.association

trans <- function(GLM){
  answer <- NULL
  for (i in 1:nrow(GLM)){
    #remove any contribution from n
    snp <- rownames(GLM)[i]
    sel <- GLM[i,grep("pVal", colnames(GLM))]
    min_pval <- as.numeric(min(sel[!is.na(sel)]))
    pval_pos <- which(GLM[i,] == min_pval[1]) #modificat el [1]
    or_pos <- abs(pval_pos[1]-1)
    or <-  as.numeric(as.character(GLM[i,or_pos]))
    allele <- as.character(GLM[i,"alleles"])
    allele<- paste(c(substr(allele, 1, 1),substr(allele, pval_pos+1, pval_pos+1)),collapse = "/")
    ans <- c(snp, or, min_pval, allele)
    answer <- rbind(answer, ans)
  }
  answer <- as.data.frame(answer)
  colnames(answer) <- c("SNPs", "OR","pval","alleles")
  rownames(answer) <- answer$SNPs
  answer <- answer[,-1]
  
  return(answer)
}
##################################################



#### FUNCTION 7: Single SNP association ####
# Calculates the OR between Individual SNPs and Lineage
Single.SNP.association <- function(SNP_geno, info_seq, method = "GLM", rare.variants = FALSE) {
  #are samples names identical?
  if(identical(as.character(rownames(SNP_geno)), as.character(rownames(info_seq)))==FALSE){
    stop("ID samples doesn't match. Use numerical ID")}
  #is the variable of study a two categorigal factor?
  if (is.factor(info_seq[,1]) == TRUE | is.integer(info_seq[,1]) == TRUE ){
    if(length(unique(info_seq[,1]))==2){
      casecon <- info_seq[,1]
      df <- cbind(casecon, SNP_geno) 
    }else{
      stop("ERROR: The variable of study must contain two categorical elements")
    }  
  }else{
    stop("ERROR: The variable of study must be a factor or integer class")
  }
  
  if (rare.variants == FALSE){ #IF FALSE, mvf = 0.05, no rare-variants included
    mvf_parameter <- 0.05
  }else{ #If TRUE, mvf = 0.01, rare variants included
    mvf_parameter <- 0.01
  }
  #Prepare data for the analysis
  SETS <- list() 
  info_SNP <- NULL
  for (i in 2:ncol(df)) { 
    set <- as.data.frame(cbind(df[,1],as.character(df[,i])))
    set <- set[!set[,2]=="n",]
    set$V2 <- droplevels(set$V2)
    #borrar nt si  MAF < 0.05
    mvf <- prop.table(table(set[,2]))
    if (length(grep(TRUE,mvf < mvf_parameter))>=1) { #if there is a TRUE, then we need to remove it
      pos <- grep(TRUE,mvf < mvf_parameter)
      rmv <- names(mvf)[pos]
      set <- set[!set[,2]==rmv,]
      set$V2 <- droplevels(set$V2)
    }
    colnames(set) <- c(colnames(df[1]),colnames(df[i]))
    x <- table(set[,2])
    
    if (nrow(x)==1) {
      major <- names(x)
      minor <- "rare_variant"
      set2 <- NULL
    } else if (nrow(x)==2){
      major <- names(x)[grep(max(x), x)]
      minor <- names(x)[grep(min(x), x)]
      if (length(minor)>=2){
        minor <- minor[2]
      }
      set2 <- set
    } else if (nrow(x)==3){
      ord <- order(x)
      major <- names(x)[grep(max(x), x)]
      rmv <- names(x)[grep(min(x), x)]
      minor <- names(x)[grep("2", ord)]
      if (length(rmv)>=2){
        minor <- rmv[2]
        rmv <- rmv[1]
      }
      set2 <- set[!set[,2]==rmv,]
    }
    info <- c(colnames(df[i]), major, minor)
    info_SNP <- rbind(info_SNP, info)
    SETS[[i]] <- set2
  }
  #Rewrite nicely info_SNPs and SETS
  rownames(info_SNP) <- info_SNP[,1]
  info_SNP <- as.data.frame(info_SNP)[,-1]
  colnames(info_SNP) <- c("major", "minor")
  
  names(SETS) <- colnames(df)[1:ncol(df)]
  SETS[[1]] <- NULL
  SETS[sapply(SETS, is.null)] <- NULL
  snp_names_rm <- rownames(info_SNP[grep("rare_variant", info_SNP$minor),])
  if (length(snp_names_rm) > 0){
    info_SNP <- info_SNP[-grep("rare_variant", info_SNP$minor),]
  }
  df <- df[,(colnames(df) %in% snp_names_rm)==FALSE]
  ####### remove info_SNP with rare-variants. 
  #Chose method: 
  if(method == "OR"){
    library(epitools)
    show(paste ("Reference levels -", paste(levels(as.factor(df[,1])), collapse = " vs "), sep = " "))
    OR <- NULL
    for(j in 1:length(SETS)){
      set <- SETS[[j]]
      maj <- as.character(info_SNP[j,1])
      table <- table(set)
      table <- t(table)
      if(nrow(table) == 2){ 
        if(rownames(table)[1]==maj){
          reverse = "neither"
        }else{
          reverse = "rows"
        }
        or <- oddsratio.fisher(table, rev = reverse)
        alleles <- paste(rownames(or$p.value), collapse = "/")
        val <- c(or$measure[2,1], or$p.value[2,1], 0,0, alleles) 
        OR <- rbind(OR, val)
        
      }else if(nrow(table) == 3){ 
        table <- table[order(-table[,2]),]
        or <- oddsratio.fisher(table[,1:2])
        alleles <- paste(rownames(or$p.value), collapse = "/")
        val <- c(or$measure[2,1], or$p.value[2,1], or$measure[3,1], or$p.value[3,1], alleles) 
        OR <- rbind(OR, val)
      }
    }
    #make nicer the OR data.frame
    OR <- replace(OR, OR == 0, NA)
    min_pval <- as.numeric(apply(OR, 1, function(x) min(x, na.rm = TRUE)))
    OR_short <- NULL
    for (l in 1:nrow(OR)) {
      x <- grep(min_pval[l], as.matrix(OR[l,]))+1
      minor <- substr((OR[l,5]),x,x)
      major <- substr((OR[l,5]),1,1)
      alleles <- paste(major, minor, sep ="/")
      pval <- OR[l,(x-1)]
      est <- OR[l,(x-2)]
      or <- c(est, pval, alleles)
      OR_short <- rbind(OR_short, or)
    }             
    rownames(OR_short) <- names(SETS)
    colnames(OR_short) <- c("OR", "pval", "alleles")
    OR <- as.data.frame(OR_short)
  }else if (method == "GLM"){
    show(paste ("Reference levels -", paste(levels(as.factor(df[,1])), collapse = " vs "), sep = " "))
    OR_list <- list()
    for(k in 1:length(SETS)){ 
      set <- SETS[[k]]
      allele <-  colnames(table(set))
      maj <- as.character(info_SNP[k,1])
      other <- allele[-grep(maj, allele)]
      set[,2] <- factor(set[,2], levels = c(maj, other))
      model <- glm(set[,1] ~ set[,2], data = set,family=binomial(link="logit"))
      pval <- t(as.matrix(coef(summary(model))[,4]))
      while(length(pval)!=4){
        pval <- c(pval, NA)
      }
      estimate <- t(as.matrix(exp(coef(model))))
      while(length(estimate)!=4){
        estimate <- c(estimate, NA)
      }
      other_alleles <- substr(rownames(coef(summary(model))),nchar(rownames(coef(summary(model)))),
                              nchar(rownames(coef(summary(model)))))[-1]
      alleles <- paste(maj, paste(other_alleles, collapse = "/"), sep= "/")
      or <- c(estimate, pval, alleles)
      OR_list[[k]] <- or
    }
    names(OR_list) <- colnames(df)[2:ncol(df)]
    GLM <- data.frame(do.call(rbind, OR_list))
    colnames(GLM) <- c("I","OR", "OR.2", "OR.3","I2", "pVal", "pVal.2", "pVal.3","alleles")
    GLM <- as.data.frame(GLM[,c(2,6,3,7,4,8,9)])
    OR <- trans(GLM)
  }else{
    stop("Unknown method: Implemented methods are OR and GLM")
  }
  return(OR)
}

#### FUNCTION 8: Plot Confussion Matrix ####
# Shows the confusion matrix as a table with green (matches) and red (mismatches) colors

plot.cm <- function(predicted, real, rescales = T, xlab = "real", ylab = "predicted"){
  if(length(table(is.na(predicted))) == 2){
    predicted <- as.character(predicted)
    predicted[is.na(predicted)] <- "n"
    predicted <- as.factor(predicted) 
    # predicted <- droplevels(predicted)
  }
  if(length(table(is.na(real))) == 2){
    real <- as.character(real)
    real[is.na(real)] <- "n"
    real <- as.factor(real)
    # real <- droplevels(real)
  }
  
  suppressMessages(library(mlearning))
  Conf <- confusion(predicted, real)
  if (rescales == TRUE){
    prior(Conf) <- 100
  }
  # The above rescales the confusion matrix such that columns sum to 100.
  opar <- par(mar=c(5.1, 6.1, 2, 2))
  x <- x.orig <- unclass(Conf)
  x <- log(x + 0.5) * 2.33
  x[x < 0] <- NA
  x[x > 10] <- 10
  diag(x) <- -diag(x)
  image(1:ncol(x), 1:ncol(x),
        -(x[, nrow(x):1]), xlab= xlab, ylab='',
        col=colorRampPalette(c(hsv(h = 0, s = 0.9, v = 0.9, alpha = 1), 
                               hsv(h = 0, s = 0, v = 0.9, alpha = 1), 
                               hsv(h = 2/6, s = 0.9, v = 0.9, alpha = 1)))(41), 
        xaxt='n', yaxt='n', zlim=c(-10, 10))
  axis(1, at=1:ncol(x), labels=colnames(x), cex.axis=0.8)
  axis(2, at=ncol(x):1, labels=colnames(x), las=1, cex.axis=0.8)
  title(ylab= ylab, line=4.5)
  abline(h = 0:ncol(x) + 0.5, col = 'gray')
  abline(v = 0:ncol(x) + 0.5, col = 'gray')
  nn <- length(levels(predicted))
  text(1:nn, rep(nn:1, each=nn), 
       labels = sub('^0$', '', round(c(x.orig), 0)))
  box(lwd=2)
  par(opar)
}

#### USAGE ####
# plot.cm(predicted, real, rescales = F, xlab = "MLT", ylab = "RF")
####################

