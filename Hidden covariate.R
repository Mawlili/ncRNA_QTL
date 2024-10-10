########################################################################
# parse arguments
########################################################################

args <- commandArgs(trailingOnly = TRUE)
sim_num <- as.numeric(args[1])

########################################################################
# load dependencies
########################################################################

library(data.table)
library(MASS)
library(pracma)

### HCP 80 with seqPCs
hidden_convariate_linear <- function(F,Y,k,lambda,lambda2,lambda3,iter) {
  ## Use Example
  # hcp = hidden_convariate_linear(standardize(datSeq), standardize(t(datExpr)),k=10,iter = 100)
  #
  # function [Z,B,U,o,error,error1,error2,dz,db,du] = hidden_covariate_linear(F,Y,k,lambda,lambda2,lambda3,iter);
  # input:
  #      F: a matrix nxd of known covariates, where n is the number of
  #      subjects and d is the number of known covariates. *must be standardize (columns have 0 mean and constant SS).
  #      Y: a matrix of nxg of expression data (must be standardized (columns
  #      scaled to have constant SS and mean 0). ** use standardize function to standardize F and Y.
  #      k: number of inferred hidden components (k is an integer)
  #      lambda, lambda2, lambda3 are model parameters
  #      (optional) iter: number of iterations (default = 100);
  #
  #      note: k>0, lambda>0, lambda2>0, lambda3>0 must be set by the user based on the data at hand. one can set these values
  #      using cross-validation, by evaluating the "performance" of the  resulting residual data on a desired task.
  #      typically, if lambda>5, then hidden factors match the known covariates closely.
  #
  # objective:
  #
  # this function solves the following problem:
  # argmin_{Z,B,U}   ||Y-Z*B||_2 + \lambda*||Z-F*U||_2 + \lambda2*||B||_2 + \lambda_3||U||_2
  #
  # output:
  #      Z: matrix of hidden components, dimensionality: nxk
  #      B: matrix of effects of hidden components, dimensionality: kxg
  #      o: value of objective function on consecutive iterations.
  #
  # to use the residual data: Residual = Y - Z*B

  
  tol = 1e-6;
  
  U = matrix(0, nrow=dim(F)[2],k)
  Z = matrix(0, nrow=dim(F)[1],k)
  B = matrix(runif(dim(Z)[2]*dim(Y)[2]), nrow=dim(Z)[2], ncol=dim(Y)[2])
  F = as.matrix(F)
  
  n1 = dim(F)[1]
  d1 = dim(F)[2]
  
  n2 = dim(Y)[1]
  d2 = dim(Y)[2]
  
  if(n1!=n2)    stop("number of rows in F and Y must agree")
  
  if (k<1 | lambda<1e-6 | lambda2<1e-6 | lambda3<1e-6 ) {
    stop("lambda, lambda2, lambda3 must be positive and/or k must be an integer");
  }
  
  o = vector(length=iter)
  
  for (ii in 1:iter) {
    print(ii)
    o[ii] = sum((Y - Z%*%B)^2) + sum((Z -  F%*%U)^2)*lambda + 
      (sum(B^2))*lambda2 + lambda3*(sum(U^2));
    Z = (Y %*% t(B) + lambda * F %*%U) %*% ginv(B %*% t(B) + lambda * diag(dim(B)[1]))
    B = mldivide(t(Z) %*% Z + lambda2 * diag(dim(Z)[2]), (t(Z) %*% Y))
    U = mldivide(t(F) %*% F * lambda + lambda3 * diag(dim(U)[1]), lambda * t(F) %*% Z)
    
    if(ii > 1 &&  (abs(o[ii]-o[ii-1])/o[ii]) < tol)  break
  }
  
  error =  sum((Y - Z%*%B)^2) / sum(Y^2)  + sum((Z - F%*%U)^2)/sum((F%*%U)^2)
  error1 = sum((Y - Z%*%B)^2) / sum(Y^2);
  error2 = sum((Z - F%*%U)^2) / sum((F%*%U)^2);
  
  dz = Z%*%(B%*%t(B) + lambda*diag(dim(B)[1]))-(Y%*%t(B) + lambda*F%*%U);
  db = (t(Z)%*%Z + lambda2*diag(dim(Z)[2]))%*%B - t(Z)%*%Y;
  du = (t(F)%*%F*lambda + lambda3*diag(dim(U)[1]))%*%U-lambda*t(F)%*%Z;
  
  
  dataout = list(Z = Z, B = B, U = U)
  return(dataout)
}


standardize<- function(X){
  X = as.matrix(X)
  # n = dim(X)[1]
  # p = dim(X)[2]
  
  X = scale(X, center = TRUE, scale = FALSE)
  # X = scale(X,center=FALSE, scale=sqrt(apply(X^2,2,sum)))
  
  # m = apply(X,2,mean)
  # st = sqrt(apply(X^2,2,sum));
  # st_mat = matrix(st, nrow = length(st), ncol = dim(X)[2], byrow=FALSE)
  # X2 = X / st_mat
  return (X)
}


gene_bed <- fread("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/TCGA_BRCA_BED/named_lifted_log2_input_TCGA_BRCA.bed")
gene_cbt = as.matrix(bed_data[,7:ncol(gene_bed)])

#for (k in seq(100,300,by=25))
  setwd('/rsrch5/home/epi/bhattacharya_lab/users/whwu1/temp/hidden_cov')
  #file.remove('/u/scratch/a/abtbhatt/bigbrain_rnaseqqc/all_seq_stats.tsv')
  system('cat flipped* >> /rsrch5/home/epi/bhattacharya_lab/data/TCGA/BRCA/merged_cov_with_barcodes.tsv')
  seqStats = fread('/rsrch5/home/epi/bhattacharya_lab/data/TCGA/BRCA/merged_cov_with_barcodes.tsv')
  #seqStats = subset(seqStats,
  #                  individualID != 'individualID')
  seqStatMat = as.matrix(seqStats[, -c(1:5)])
  class(seqStatMat) = 'numeric'
  rownames(seqStatMat) = seqStats$barcode
  seqStatMat = seqStatMat[,which(apply(seqStatMat,2,var) != 0)]
  seqStatMat = seqStatMat[colnames(gene_bed)[7:ncol(gene_bed)],]

  F <- standardize(seqStatMat)
  Y <- standardize(t(gene_cbt))


  
# HCP parameters to iterate through 
k_values <- c(10,25,50) # number of hidden components to be estimated
lambda_values <- c(0.01,0.05,1,5)
lambda2_values <- c(0.01,0.05,1,5)
lambda3_values <- c(0.01,0.05,1,5)
parameter_space <- expand.grid(k_values,
                               lambda_values,
                               lambda2_values,
                               lambda3_values)
  
# assign parameters based on simulation number ( should be a number between 1 and nrow(parameter_space) )
# simulation number passed via command line ( see line 16 )

k_values = parameter_space[sim_num,1] 
lambda_values = parameter_space[sim_num,2]
lambda2_values = parameter_space[sim_num,3]
lambda3_values = parameter_space[sim_num,4]

# function to perform cross-validation
cross_validate <- function(F, Y, k_values, lambda_values, lambda2_values, lambda3_values, iter = 100, n_folds = 5, seed=555) {
  set.seed(seed)  # For reproducibility
  folds <- sample(rep(1:n_folds, length.out = nrow(F)))
  
  results <- data.frame()
  
  # loop through all parameter combinations (should really just be 1 in this script)
  for (k in k_values) {
    for (lambda in lambda_values) {
      for (lambda2 in lambda2_values) {
        for (lambda3 in lambda3_values) {
          
          errors <- numeric(n_folds)
          
          for (fold in 1:n_folds) {
            # split the data into training and validation sets
            train_idx <- which(folds != fold)
            val_idx <- which(folds == fold)
            
            F_train <- F[train_idx, ]
            Y_train <- Y[train_idx, ]
            F_val <- F[val_idx, ]
            Y_val <- Y[val_idx, ]
            
            # train the model on the training set
            model <- hidden_convariate_linear(F_train, Y_train, k = k, lambda = lambda, lambda2 = lambda2, lambda3 = lambda3, iter = iter)
            
            # predict on the validation set
            Z_val <- (Y_val %*% t(model$B) + lambda * F_val %*% model$U) %*% ginv(model$B %*% t(model$B) + lambda * diag(dim(model$B)[1]))
            error_val <- sum((Y_val - Z_val %*% model$B)^2) / sum(Y_val^2)
            
            errors[fold] <- error_val
          }
          
          # save results
          avg_error <- mean(errors)
          results <- rbind(results, data.frame(k = k, lambda = lambda, lambda2 = lambda2, lambda3 = lambda3, avg_error = avg_error))
        }
      }
    }
  }
  
  return(results[order(results$avg_error), ])
}

# perform cross-validation
cv_results <- cross_validate(F, Y, k_values, lambda_values, lambda2_values, lambda3_values, iter = 100, n_folds = 5, seed = sim_num)

# save the results
save(cv_results,file=paste0("5fold_CV_res_",sim_num,".RData"))


  
 # outFolder = '/rsrch5/home/epi/bhattacharya_lab/users/whwu1/temp/hidden_cov'
  #hcp <- hidden_convariate_linear(standardize(seqStatMat), 
  #                                standardize(t(gene_cbt)), 
   #                               lambda=5,lambda2=1, lambda3=1, 
   #                               k=k, iter=100)
  

  
