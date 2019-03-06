#'TNBC consensus molecular subtype prediction
#'
#'Predicts the TNBC consensus molecular subtype of TNBC samples.
#'
#'@param expr A \code{SummarizedExperiment} object or a matrix containig gene
#'expression profiles. If input is a \code{SummarizedExperiment}, the first
#'element in the assays list should be a matrix of gene expression.
#'Rows and columns of the gene expression matrix correspond to genes and
#'samples, respectively (rownames must be to gene symbols).
#'@return A vector of assigned subtypes.
#'@export
#'@importFrom e1071 impute
#'@importFrom stats predict median
#'@importFrom methods is
#'@import SummarizedExperiment
#'@examples
#'# Load gene expression profiles of TNBC samples
#'data(GSE25055)
#'
#'# Predict consensus molecular subtypes of TNBC samples
#'prediction <- predictCMS(expr = GSE25055)
#'table(prediction)
predictCMS <- function(expr){

  if (is(expr, "SummarizedExperiment")){
    exp.mat <- assays(expr)[[1]]
  } else{
    exp.mat <- expr
  }

  #Impute missing values
  missings <- genelist[!genelist %in% rownames(exp.mat)]
  nmissing <- length(missings)
  if(nmissing > 0){
    missing.mat <- matrix(NA, nrow = nmissing, ncol = ncol(exp.mat))
    rownames(missing.mat) <- missings
    colnames(missing.mat) <- colnames(exp.mat)
    input <- rbind(exp.mat, missing.mat)
    input <- impute(input, what = "median")
  } else {
    input <- exp.mat
  }

  #Adjust scales of input data
  input <- log2(input + 1)
  input <- input - apply(input, 1, median)
  input <- data.frame(t(input[genelist,]), check.names = FALSE)

  #Predict CMS of new samples
  pred <- predict(SVM.model, input, probability = TRUE)
  prob <- attr(pred, "probabilities")
  attr(pred, "probabilities") <- prob[,c("MSL", "IM", "LAR", "SL")]

  return(pred)

}
