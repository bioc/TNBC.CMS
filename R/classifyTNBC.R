#'TNBC classification
#'
#'This function assigns the TNBC consensus molecular subtype to TNBC samples.
#'
#'@param mat A matrix of gene expression with genes in rows and samples in columns (rownames corresopnding to gene symbol).
#'@return A vector of assigned subtypes.
#'@export
#'@importFrom sva ComBat
#'@importFrom e1071 impute
#'@importFrom stats predict median
classifyTNBC <- function(mat){

  #Avoid test-set bias
  commongenes = intersect(rownames(exprPublic), rownames(mat))
  input = cbind(exprPublic[commongenes, ], mat[commongenes, ])
  output = sva::ComBat(dat = as.matrix(input), batch = c(rep(0, ncol(exprPublic)), rep(1, ncol(mat))),
                       ref.batch = 0)
  mat = output[, ((ncol(exprPublic)+1):ncol(output))]

  #Impute missing values
  missings = rownames(mat)[!rownames(mat) %in% genelist]
  if(length(missings) > 0){
    mat = rbind(mat, matrix(NA, nrow = length(missings), ncol = ncol(mat),
                            dimnames = list(missings, colnames(mat))))
    mat = impute(mat, what = "median")
    mat = mat[genelist,]
  }

  #Predict CMS of new samples
  mat = mat - apply(mat, 1, median)
  predict(svmModel, data.frame(t(mat), check.names = F))

}
