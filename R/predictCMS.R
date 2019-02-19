#'TNBC consensus molecular subtype prediction
#'
#'This function predicts the TNBC consensus molecular subtype of TNBC samples.
#'
#'@param exp.mat A matrix of gene expression with genes in rows and samples in columns (rownames corresopnding to gene symbol).
#'@return A vector of assigned subtypes.
#'@export
#'@importFrom e1071 impute
#'@importFrom stats predict median
predictCMS <- function(exp.mat){

  #Impute missing values
  missings = genelist[!genelist %in% rownames(exp.mat)]
  if(length(missings) > 0){
    input = rbind(exp.mat, matrix(NA, nrow = length(missings), ncol = ncol(exp.mat),
                                          dimnames = list(missings, colnames(exp.mat))))
    input = impute(input, what = "median")
  } else {
    input = exp.mat
  }

  #Adjust scales of input data
  input = log2(input + 1)
  input = input - apply(input, 1, median)

  #Predict CMS of new samples
  pred = predict(SVM.model, data.frame(t(input[genelist,]), check.names = F), probability = T)
  attr(pred, "probabilities") = attr(pred, "probabilities")[,c("MSL", "IM", "LAR", "SL")]

  return(pred)

}
