#prevent t-test errors from crashing
my.t.test <- function(...) {
    obj<-try(t.test(...), silent=TRUE)
    if (is(obj, "try-error")) {
		obj$p.value=NA
		return(obj) 
	}else{ 
		return(obj)
	}
}

cor.test.somevsall <- function(mat2, rows, method) {
  sourceCpp('deleterow.cpp')      
  mat <- deleterow(as.matrix(mat2));
        N = length(mat[,1]);
        M = length(rows);
	library('doParallel'); # Works perfectly well but obviously want it to be more comaptible with BioConductor
	n_cores<-detectCores();
	cl<-makeCluster(n_cores);
	registerDoParallel(cl)
  
 x<-    foreach (i=1:M, .combine='rbind', .packages='foreach') %dopar% {
                row = rows[i];
        foreach (j=1:N, .combine='c') %dopar% {
                test <- cor.test(unlist(mat[row,]),unlist(mat[j,]), method=method)
                if (i == j){data.frame(pv=1, st=test$estimate)} # Gene does not target itself
                else{data.frame(pv=test$p.val, st=test$estimate)}
        }
  }

        if(M == 1){x <- data.frame(x)} #output is vector otherwise
        stopCluster(cl)
        pvals <- x[,grep("pv", colnames(x))];
		    strength <- x[,grep("st", colnames(x))];
        output = list();
        output$pvals = pvals;
        output$strength = strength;
        colnames(output$strength) = rownames(mat)
        rownames(output$strength) = rownames(mat)[rows]
        colnames(output$pvals) = rownames(mat)
        rownames(output$pvals) = rownames(mat)[rows]
        output
}


# Even with or instead of and always get nothing, is something wrong with the pdcors? -> needed to unlist each row of data.
cor.test_pair_interaction <- function (x,method,Adj,Mat) {
  source("ppcor.R")
  gene1 = x[1];
  gene2 = x[2];
  
  gene1_rows = which(as.character(Adj[,1]) == gene1)
  gene2_rows = which(as.character(Adj[,1]) == gene2)
  
  gene1_targets = Adj[gene1_rows,2];
  gene2_targets = Adj[gene2_rows,2];
  
  sharedtargets = unique(gene1_targets[gene1_targets %in% gene2_targets])
  interact_targets = vector();
  for (t in sharedtargets) {
    dcor_1t <- Adj[(Adj[,1] == gene1 & Adj[,2] == t),]$strength
    pdcor1t_2 <- pcor.test(unlist(Mat[rownames(Mat)==gene1,]), unlist(Mat[rownames(Mat)==t,]),unlist(Mat[rownames(Mat)==gene2,]), method=method)$estimate
    dcor_2t <- Adj[(Adj[,1] == gene2 & Adj[,2] == t),]$strength
    pdcor2t_1 <- pcor.test(unlist(Mat[rownames(Mat)==gene2,]), unlist(Mat[rownames(Mat)==t,]),unlist(Mat[rownames(Mat)==gene1,]), method=method)$estimate
    
    if (abs(dcor_1t) < abs(pdcor1t_2) & abs(dcor_2t) < abs(pdcor2t_1)) {
      # INTERACTION
      if (pdcor1t_2 > 0 & pdcor2t_1 > 0) {
        interact_targets = c(interact_targets,paste("+",t,sep=""));
      } else if (pdcor1t_2 < 0 & pdcor2t_1 < 0) {
        interact_targets = c(interact_targets,paste("-",t,sep=""));
      } else {
        interact_targets = c(interact_targets,t);
      }
    }
  }
  if (length(interact_targets) >= 1) {
    return(list(pair = c(gene1,gene2), targets = interact_targets));
  } else {
    return(NA)
  }
}


dcor.test.somevsall <- function(mat, rows) {
  library(Rcpp)
  sourceCpp('deleterow.cpp')
  mat2=deleterow(mat);
  N = length(mat2[,1]);
  M = length(rows);
  library('doParallel');
  n_cores<-detectCores()-1;
  cl<-makeCluster(n_cores);
  registerDoParallel(cl)
  x <- foreach(i=1:M, .combine='rbind', .packages='foreach') %dopar% {
    row = rows[i];
    
    bins = quantile(mat[row,],c(0.25,0.75))
    low = unlist(quantile(mat[row,],c(0.25)))
    high = unlist(quantile(mat[row,],c(0.75)))
    
    foreach(j=1:N, .combine='c') %dopar% {
      require("energy")
      require("pdcor2")
      test = dcor.ttest(unlist(mat2[row,]),unlist(mat2[j,]))
      
      my.t.test <- function(...) {
        obj<-try(t.test(...), silent=TRUE)
        if (is(obj, "try-error")) {
          obj$p.value=NA
          return(obj) 
        }else{ 
          return(obj)
        }
      }
      
      test_dir = my.t.test(mat2[j,(mat2[row,]<=low)],mat2[j,(mat2[row,]>=high)])
      if (j == row) {
        if (is.na(test_dir$p.value)) {
          data.frame(pv=1,st=test$estimate,d=0);
        } else {
          if (low != high & test_dir$p.value < 0.05) {
            if (test_dir$estimate[1] < test_dir$estimate[2]) {
              data.frame(pv=1,st=test$estimate,d=1);
            } else {
              data.frame(pv=1,st=test$estimate,d=-1);
            }
          } else {
            data.frame(pv=1,st=test$estimate,d=0);
          }
        }
      }# Gene no target itself.
      else{
        if (is.na(test_dir$p.value)) {
          data.frame(pv=test$p.val,st=test$estimate,d=0);
        } else {
          if (low != high & test_dir$p.value < 0.05) {
            if (test_dir$estimate[1] < test_dir$estimate[2]) {
              data.frame(pv=test$p.val,st=test$estimate,d=1);
            } else {
              data.frame(pv=test$p.val,st=test$estimate,d=-1);
            }
          } else {
            data.frame(pv=test$p.val,st=test$estimate,d=0);
          }
        }
      }
    }
  }
  
  stopCluster(cl)
  pvals<-x[,grep("pv",colnames(x))]
  strength<-x[,grep("st",colnames(x))]
  direction<-x[,grep("d",colnames(x))]
  output = list();
  output$pvals = pvals;
  output$strength = strength;
  output$dir = direction
  colnames(output$strength) = rownames(mat2)
  rownames(output$strength) = rownames(mat2)[rows]
  colnames(output$pvals) = rownames(mat2)
  rownames(output$pvals) = rownames(mat2)[rows]
  colnames(output$dir) = rownames(mat2)
  rownames(output$dir) = rownames(mat2)[rows]
  output
}

# Even with or instead of and always get nothing, is something wrong with the pdcors? -> needed to unlist each row of data.
dcor.test_pair_interaction <- function (x, Adj, Mat, margin = 0) {
    gene1 = x[1]
    gene2 = x[2]
    gene1_rows = which(as.character(Adj[,1]) == gene1)
    gene2_rows = which(as.character(Adj[,1]) == gene2)
    
    gene1_targets = Adj[gene1_rows,2];
    gene2_targets = Adj[gene2_rows,2];
    
    sharedtargets = unique(gene1_targets[gene1_targets %in% gene2_targets])
    sharedtargets <- as.matrix(table(factor(sharedtargets, exclude =0)));
    interact_targets = vector();
    pathway_targets = vector(); seen = 0;
    
    if(ncol(sharedtargets) > 0){
        #library('doParallel');
        #n_cores<-detectCores()-1;
        #cl<-makeCluster(n_cores);
        #registerDoParallel(cl)
    
    for(i in 1:ncol(sharedtargets)){
        
        t <- rownames(sharedtargets[i]);
        if (sum(as.character(Adj[,1]) == gene1 & as.character(Adj[,2]) == t) > 0 & sum(as.character(Adj[,1]) == gene2 & as.character(Adj[,2]) == t) > 0)
        {
            if ( gene1 != gene2 & gene1 != t & gene2 != t) {
                dcor_1t = Adj[(Adj[,1] == gene1 & Adj[,2] == t),]$strength
                dcor_2t = Adj[(Adj[,1] == gene2 & Adj[,2] == t),]$strength
                
                library('Rcpp')
                sourceCpp('partialdcor.cpp')
                
                Mat_1 <- as.matrix(unlist(Mat[rownames(Mat)==gene1,]));
                Mat_2 <- as.matrix(unlist(Mat[rownames(Mat)==gene2,]));
                Mat_t <- as.matrix(unlist(Mat[rownames(Mat)==t,]));
                
                tmpdataframe <- pair_interaction(t, Adj, Mat_1, Mat_2, Mat_t,
                dcor_1t, dcor_2t, margin = 0, gene1, gene2, seen);
                
                interact_target <- tmpdataframe$interact_target;
                pathway_target <- tmpdata$pathway_target;
                seen <- tmpdataframe$seen;
                first <- tmpdataframe$first;
                second <- tmpdataframe$second;
                
                interact_targets = c(interact_target, interact_targets);
                pathway_targets = c(pathway_target, pathway_targets);
            }
        }
    }
      interact_targets <- na.omit(interact_targets);
      pathway_targets <- na.omit(pathway_targets);
    
    
    if ((length(pathway_targets) >= 1 & seen == 1 ) & length(interact_targets) >= 1) {
        return(list(complicated = c(first,second),targets=c(pathway_targets,interact_targets)));
    } else if (length(pathway_targets) >= 1 & seen == 1) {
        return(list(pathway = c(first, second), targets = pathway_targets));
    } else if (length(interact_targets) >= 1) {
        return(list(pair = c(gene1,gene2), targets = interact_targets));
    } else {
        return(NA)
    }
    #stopCluster(cl)
    }else{
        return(NA)
    }
}


dcor.test_pair_pathway <- function (x, Adj, Mat, margin=0.001){
    gene1 = x[1]
    gene2 = x[2]
    gene1_rows = which(as.character(Adj[,1]) == gene1)
    gene2_rows = which(as.character(Adj[,1]) == gene2)
    
    gene1_targets = Adj[gene1_rows,2];
    gene2_targets = Adj[gene2_rows,2];
    
    sharedtargets = unique(gene1_targets[gene1_targets %in% gene2_targets])
    sharedtargets = factor(sharedtargets, exclude = 0)
    interact_targets = vector();
    seen = 0;
    
    for (t in sharedtargets) {
        
        dcor_1t = Adj[(Adj[,1] == gene1 & Adj[,2] == t),]$strength
        dcor_2t = Adj[(Adj[,1] == gene2 & Adj[,2] == t),]$strength
        
        library('Rcpp')
        sourceCpp('partialdcor.cpp')
        
        Mat_1 <- as.matrix(unlist(Mat[rownames(Mat)==gene1,]));
        Mat_2 <- as.matrix(unlist(Mat[rownames(Mat)==gene2,]));
        Mat_t <- as.matrix(unlist(Mat[rownames(Mat)==t,]));
        
        tmpdataframe <- pair_pathway(t, Adj, Mat_1, Mat_2, Mat_t,
        dcor_1t, dcor_2t, margin = 0.001, gene1, gene2, seen)
        
        interact_target <- tmpdataframe$interact_target;
        seen <- tmpdataframe$seen;
        first <- tmpdataframe$first;
        second <- tmpdataframe$second;
        
        interact_targets = c(interact_target, interact_targets);
    }
    interact_targets = na.omit(interact_targets);
    
    if (length(interact_targets) >= 1 & seen == 1) { #only report those with consistent direction
        return(list(pair = c(first,second), targets = interact_targets));
    } else {
        return(NA)
    }
}

do_interactions <- function(Mat,candidates, cor_type, multi_test, margin=0) {
  candidates = as.character(candidates);
  candidates.rows = which(rownames(Mat) %in% candidates);
  candidates.rows = candidates.rows[!is.na(candidates.rows)];
  if (length(candidates.rows) < length(candidates)) {
    warning(paste("Only ",length(candidates.rows)," candidates have expression data.\n", sep=""))
  }
  if (length(candidates.rows) == 0) {stop("No expression data for any candidate genes. Are you sure gene IDs are consistent?");}
  
  # Get Pair-wise dependencies
  if (cor_type == "dcor") {
    out <- dcor.test.somevsall(Mat,candidates.rows)
  } else if (grep("cor=", cor_type))  {
    method = unlist(strsplit(cor_type,"="));
    out <- cor.test.somevsall(Mat,candidates.rows,method[2]);
  } else {
    warning("Did not recognize specified correlation method, using default (dcor)")
    out <- dcor.test.somevsall(Mat,candidates.rows)
  }
  # Tidy up result (find significant interactions)
  pvals <- as.matrix(out$pvals)
  strength <- as.matrix(out$strength) 
  direction <- out$dir
  
  # Apply Multiple Testing Correction
  if (multi_test == "bon") {
    Sig <- which(pvals < 0.05/length(pvals[1,]),arr.ind=T)
  } else if (grep("fdr=", multi_test))  {
    method = unlist(strsplit(multi_test,"="));
    tmp <- pvals[p.adjust(pvals,method="fdr") < method[2]];
    Sig <- which(pvals <= max(tmp), arr.ind=T);
  } else if (grep("str=", multi_test)) {
    method = unlist(strsplit(multi_test,"="));
    Sig <- which(strength > method[2],arr.ind=T);
  } else {
    warning("Did not recognize specified multiple-testing correction, using default (bon)")
    Sig<- which(pvals < 0.05/length(pvals[1,]),arr.ind=T)
  }
  
  # recall cor.test.somevsall does not output a direction matrix
  Sig=data.matrix(Sig)
  library('Rcpp');
  sourceCpp('adj.cpp');
  if(cor_type == "dcor"){
    Adj <- adj_with_direction(Sig,pvals,strength,direction);
  }else if (grep("cor=", cor_type)) {
    Adj <- adj_without_direction(Sig,pvals,strength)
  }else{
    Adj <- adj_with_direction(Sig,pvals,strength,direction);
  }
  
  Adj <- na.omit(Adj); #interaction calculations cannot proceed with NAs

  # Get Interactions by conditioning pairwise interactions
  #	on each other candidate gene
  pairs <- t(combn(candidates,2))
  pairs.list <- split(pairs, seq(nrow(pairs)))
  library(BiocParallel)
  param <- MulticoreParam(workers = 3);

  if (cor_type == "dcor") {
    interactions <- bplapply(pairs.list,dcor.test_pair_interaction,Adj=Adj,Mat=Mat,margin=margin, BPPARAM = param)
  } else if (grep("cor=", cor_type))  {
    #method = unlist(strsplit(cor_type,"="));
    interactions <- bplapply(pairs.list,cor.test_pair_interaction,method=method[2],Adj=Adj,Mat=Mat, BPPARAM = param)
  } else {
    warning("Did not recognize specified partial correlation method, using default (pdcor)")
    interactions <- bplapply(pairs.list,dcor.test_pair_interaction,Adj=Adj, BPPARAM = param)
  }

  #Clean & Return results
  interactions <- matrix(unlist(interactions), ncol = 2, byrow = TRUE, use.names = FALSE)
  interactions <- interactions[!is.na(interactions)]
  tmp <- list(Adj = Adj,interactions = interactions)
  saveRDS(tmp, 'do_interactions_output.rds')
  return(tmp);
}

Adj_intersection <- function(list_matrices){
  #Input list of adjacency matrices, returns their intersection
  #To be used when comparing significant results over different tests 
  n = length(list_matrices);
  if(n > 2){
  M <- rbind(list_matrices[[1]], list_matrices[[2]]);
  for (i in 1:n-2){
    M <- rbind(M, list_matrices[[i+2]]);
  }
  }
  else{
    M <- rbind(list_matrices[[1]], list_matrices[[2]]);
  }
  return(unique(M))
}

Adj_not_intersection <- function(list_matrices){
    #Input list of adjacency matrices, returns list of same length giving those not in the intersection that appeared in each list
    M <- Adj_intersection(list_matrices);
    N <- list();
    n =length(list_matrices);
    for (i in 1:n){
        tmp <- listmatrices[[i]][-rownames(M),];
    }
}
