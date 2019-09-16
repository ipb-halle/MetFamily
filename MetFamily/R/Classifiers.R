
colSumsPos_classifier <- list(
  train = function(matrix_train, classes_pm_train){
    posRows <- classes_pm_train=="+"
    colSumPos <- Matrix::colSums(x = matrix_train[ posRows, ])
    colSumNeg <- Matrix::colSums(x = matrix_train[!posRows, ])
    
    colSums_PosNeg <- (colSumPos / sum( posRows)) - (colSumNeg / sum(!posRows))
    colSums_PosNeg[colSums_PosNeg < 0] <- 0
    #names(colSums_PosNeg) <- colnames(matrix_train)
    
    ## storage
    colSums_PosNeg <- unname(colSums_PosNeg)
    
    return(colSums_PosNeg)
  },
  classify = function(classifier, matrix_test){
    ## convert matrix format
    dgTMatrix <- as(matrix_test, "dgTMatrix")
    matrixRows <- dgTMatrix@i
    matrixCols <- dgTMatrix@j
    matrixVals <- dgTMatrix@x
    stMatrix_test <- simple_triplet_matrix(
      i    = matrixRows + 1, 
      j    = matrixCols + 1, 
      v    = matrixVals, 
      nrow = nrow(matrix_test), 
      ncol = ncol(matrix_test)
    )
    
    ## score test matrix row by row
    scores <- rowapply_simple_triplet_matrix(x = stMatrix_test, FUN = function(x){
      sum(x * classifier)
    })
    
    return(scores)
  }
)

predict_ColSums <- function(matrix_train, classes_pm_train, matrix_test, classes_pm_test){
  posRows <- classes_pm_train=="+"
  #colSumPos <- apply(X = matrix_train[ posRows, ], MARGIN = 2, FUN = sum)
  #colSumNeg <- apply(X = matrix_train[!posRows, ], MARGIN = 2, FUN = sum)
  colSumPos <- Matrix::colSums(x = matrix_train[ posRows, ])
  colSumNeg <- Matrix::colSums(x = matrix_train[!posRows, ])
  
  colSums <- (colSumPos / sum( posRows)) - (colSumNeg / sum(!posRows))
  
  if(FALSE){
    ## top ten
    names(colSums) <- colnames(matrix_train)
    tail(x = sort(colSums), n = 10)
  }
  
  #scores <- apply(X = matrix_test, MARGIN = 1, FUN = function(x){
  #  sum(x * colSums)
  #})
  
  dgTMatrix <- as(matrix_test, "dgTMatrix")
  matrixRows <- dgTMatrix@i
  matrixCols <- dgTMatrix@j
  matrixVals <- dgTMatrix@x
  stMatrix <- simple_triplet_matrix(i = matrixRows + 1, j = matrixCols + 1, v = matrixVals, nrow=nrow(matrix_test), ncol=ncol(matrix_test))
  
  scores <- rowapply_simple_triplet_matrix(x = stMatrix, FUN = function(x){
    sum(x * colSums)
  })
  
  return(scores)
}
predict_CosinusDistance <- function(matrix_train, classes_pm_train, matrix_test, classes_pm_test, ratio){
  ## sum(a*b) / (sqrt(sum(a*a)) * sqrt(sum(b*b)))
  posRows <- classes_pm_train=="+"
  colSumPos <- apply(X = matrix_train[ posRows, ], MARGIN = 2, FUN = sum)
  colSumNeg <- apply(X = matrix_train[!posRows, ], MARGIN = 2, FUN = sum)
  
  scores_pos <- apply(X = matrix_test, MARGIN = 1, FUN = function(x){
    sum(x * colSumPos) / (sqrt(sum(x * x)) * sqrt(sum(colSumPos * colSumPos)))
  })
  scores_neg <- apply(X = matrix_test, MARGIN = 1, FUN = function(x){
    sum(x * colSumNeg) / (sqrt(sum(x * x)) * sqrt(sum(colSumNeg * colSumNeg)))
  })
  
  if(!ratio){
    scores <- scores_pos
  } else {
    scores <- scores_pos / scores_neg
  }
  
  return(scores)
}
predict_Prod <- function(matrix_train, classes_pm_train, matrix_test, classes_pm_test, ratio){
  ## sum(a*norm(b))
  
  posRows <- classes_pm_train=="+"
  colSumPos <- apply(X = matrix_train[ posRows, ], MARGIN = 2, FUN = sum)
  colSumNeg <- apply(X = matrix_train[!posRows, ], MARGIN = 2, FUN = sum)
  
  colSumPos <- colSumPos / sum(colSumPos)
  colSumNeg <- colSumNeg / sum(colSumNeg)
  
  scores_pos <- apply(X = matrix_test, MARGIN = 1, FUN = function(x){
    sum(x * colSumPos)
  })
  scores_neg <- apply(X = matrix_test, MARGIN = 1, FUN = function(x){
    sum(x * colSumNeg)
  })
  
  if(!ratio){
    scores <- scores_pos
  } else {
    scores <- scores_pos / scores_neg
  }
  
  return(scores)
}
predict_Jaccard <- function(matrix_train, classes_pm_train, matrix_test, classes_pm_test, ratio){
  ## sum(intersection(a,b))/sum(union(a,b))
  
  posRows <- classes_pm_train=="+"
  colSumPos <- apply(X = matrix_train[ posRows, ], MARGIN = 2, FUN = sum)
  colSumNeg <- apply(X = matrix_train[!posRows, ], MARGIN = 2, FUN = sum)
  
  scores_pos <- apply(X = matrix_test, MARGIN = 1, FUN = function(x){
    intersection <- colSumPos > 0 & x > 0
    union        <- colSumPos > 0 | x > 0
    sum(intersection) / sum(union)
  })
  scores_neg <- apply(X = matrix_test, MARGIN = 1, FUN = function(x){
    intersection <- colSumNeg > 0 & x > 0
    union        <- colSumNeg > 0 | x > 0
    sum(intersection) / sum(union)
  })
  
  if(!ratio){
    scores <- scores_pos
  } else {
    scores <- scores_pos / scores_neg
  }
  
  return(scores)
}
predict_JaccardWeighted <- function(matrix_train, classes_pm_train, matrix_test, classes_pm_test, ratio){
  ## sum(intersection(aW,bW))/sum(union(aW,bW))
  
  posRows <- classes_pm_train=="+"
  colSumPos <- apply(X = matrix_train[ posRows, ], MARGIN = 2, FUN = sum)
  colSumNeg <- apply(X = matrix_train[!posRows, ], MARGIN = 2, FUN = sum)
  
  scores_pos <- apply(X = matrix_test, MARGIN = 1, FUN = function(x){
    intersection <- colSumPos > 0 & x > 0
    union        <- colSumPos > 0 | x > 0
    sum(colSumPos[intersection]) / sum(colSumPos[union])
  })
  scores_neg <- apply(X = matrix_test, MARGIN = 1, FUN = function(x){
    intersection <- colSumNeg > 0 & x > 0
    union        <- colSumNeg > 0 | x > 0
    sum(colSumNeg[intersection]) / sum(colSumNeg[union])
  })
  
  if(!ratio){
    scores <- scores_pos
  } else {
    scores <- scores_pos / scores_neg
  }
  
  return(scores)
}

predict_Correlation <- function(matrix_train, classes_pm_train, matrix_test, classes_pm_test, ratio, corMethod = c("pearson", "kendall", "spearman"), linkage = c("single", "average", "centroid")){
  if(corMethod == "kendall")
    stop("Operation not supported")
  
  correlations <- cor(x = t(as.matrix(matrix_test)), y = t(as.matrix(matrix_train)), method = corMethod)
  
  posItems <- classes_pm_train=="+"
  negItems <- classes_pm_train=="-"
  switch(linkage,
         "single"={
           ## single linkage
           scoresPosNeg <- apply(X = correlations, MARGIN = 1, FUN = function(x){
             if(all(is.na(x)))
               return(c("+" = 0, "-" = 1))
             else{
               maxPos <- max(x[posItems], na.rm = TRUE)
               maxNeg <- max(x[negItems], na.rm = TRUE)
               return(c("+" = maxPos, "-" = maxNeg))
             }
           })
         },
         "average"={
           ## average linkage
           scoresPosNeg <- apply(X = correlations, MARGIN = 1, FUN = function(x){
             if(all(is.na(x)))
               return(c("+" = 0, "-" = 1))
             else{
               meanPos <- mean(posItems, na.rm = TRUE)
               meanNeg <- mean(negItems, na.rm = TRUE)
               return(c("+" = meanPos, "-" = meanNeg))
             }
           })
         },
         "centroid"={
           ## centroid linkage
           distPlus  <- as.matrix(dist(matrix_train[posItems, ]))
           #distMinus <- as.matrix(dist(matrix_train[negItems, ]))
           centroidPlus  <- which.min(apply(X = distPlus , MARGIN = 1, FUN = sum))
           #centroidMinus <- which.min(apply(X = distMinus, MARGIN = 1, FUN = sum))
           centroidPlus  <- which(posItems)[[centroidPlus]]
           #centroidMinus <- which(negItems)[[centroidMinus]]
           
           scoresPosNeg <- apply(X = correlations, MARGIN = 1, FUN = function(x){
             if(all(is.na(x)))
               return(c("+" = 0, "-" = 1))
             else{
               centroidPos <- x[[centroidPlus]]
               meanNeg <- mean(x[negItems], na.rm = TRUE)
               return(c("+" = centroidPos, "-" = meanNeg))
             }
           })
         },
         stop(paste("Unknown linkage (", linkage, ")!", sep = ""))
  )
  
  scores_pos <- scoresPosNeg["+", ]
  scores_neg <- scoresPosNeg["-", ]
  
  if(!ratio){
    scores <- scores_pos
  } else {
    scores <- scores_pos / scores_neg
  }
  
  return(scores)
  #return(predicted_classes_pm)
}

colSums_classifier <- list(
  train = function(matrix_train, classes_pm_train){
    posRows <- classes_pm_train=="+"
    colSumPos <- Matrix::colSums(x = matrix_train[ posRows, ])
    colSumNeg <- Matrix::colSums(x = matrix_train[!posRows, ])
    
    colSums_PosNeg <- (colSumPos / sum( posRows)) - (colSumNeg / sum(!posRows))
    #names(colSums_PosNeg) <- colnames(matrix_train)
    
    return(colSums_PosNeg)
  },
  classify = function(classifier, matrix_test){
    ## convert matrix format
    dgTMatrix <- as(matrix_test, "dgTMatrix")
    matrixRows <- dgTMatrix@i
    matrixCols <- dgTMatrix@j
    matrixVals <- dgTMatrix@x
    stMatrix_test <- simple_triplet_matrix(
      i    = matrixRows + 1, 
      j    = matrixCols + 1, 
      v    = matrixVals, 
      nrow = nrow(matrix_test), 
      ncol = ncol(matrix_test)
    )
    
    ## score test matrix
    scores <- rowapply_simple_triplet_matrix(x = stMatrix_test, FUN = function(x){
      sum(x * classifier)
    })
    
    return(scores)
  }
)