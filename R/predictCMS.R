#'TNBC consensus molecular subtype prediction
#'
#'Predicts the TNBC consensus molecular subtype of TNBC samples.
#'
#'@param exp.mat A matrix of gene expression with genes in rows and samples
#'in columns (rownames corresopnding to gene symbol).
#'@return A vector of assigned subtypes.
#'@export
#'@importFrom e1071 impute
#'@importFrom stats predict median
#'@examples
#'# Load gene expression profiles of TNBC samples
#'data(GSE25055.exprs)
#'
#'# Predict consensus molecular subtypes of TNBC samples
#'predictions <- predictCMS(exp.mat = GSE25055.exprs)
#'table(predictions)
predictCMS <- function(exp.mat){

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
