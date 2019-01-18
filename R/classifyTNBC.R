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

  #Adjust scales of input data
  commongenes = intersect(rownames(expr_core), rownames(mat))
  combat_input = cbind(expr_core[commongenes, ], mat[commongenes, ])
  batch = c(rep(0, ncol(expr_core)), rep(1, ncol(mat)))
  combat_output = ComBat(dat = as.matrix(combat_input), batch = batch, ref.batch = 0)
  mat_scaled = combat_output[, ((ncol(expr_core)+1):ncol(combat_output))]

  #Impute missing values
  missings = rownames(mat_scaled)[!rownames(mat_scaled) %in% genelist]
  if(length(missings) > 0){
    mat_scaled = rbind(mat_scaled, matrix(NA, nrow = length(missings), ncol = ncol(mat_scaled),
                            dimnames = list(missings, colnames(mat_scaled))))
    mat_scaled = impute(mat_scaled, what = "mean")
    mat_scaled = mat_scaled[genelist,]
  }

  #Predict CMS of new samples
  predict(svmModel, data.frame(t(mat_scaled), check.names = F))

}
