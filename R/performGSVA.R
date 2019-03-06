#'Gene set variation analysis
#'
#'Performs GSVA on gene sets. Also draws a heatmap representing GSVA scores.
#'
#'@param expr A \code{SummarizedExperiment} object or a matrix containig gene
#'expression profiles. If input is a \code{SummarizedExperiment}, the first
#'element in the assays list should be a matrix of gene expression.
#'Rows and columns of the gene expression matrix correspond to genes and
#'samples, respectively (rownames must be to gene symbols).
#'@param pred A vector of predicted consensus molecular subtypes.
#'@param gene.set Gene sets provided as a list. If NULL,
#'the hallmark pathway gene sets are used.
#'@param gsva.kcdf Kernel to be used in the estimation of
#'the cumulative distribution function. By default,
#'this is set to \code{"Gaussian"} which is suitable
#'for continuous expression values. If expression values
#'are counts, \code{"Poisson"} is recommended.
#'@return A matrix of GSVA enrichment scores.
#'@details This is a wrapper function of the \code{gsva}
#'function in the \code{GSVA} package to compute GSVA
#'enrichment scores per sample and produce a heatmap
#'comparing them across consensus molecular subtypes.
#'@importFrom GSVA gsva
#'@importFrom pheatmap pheatmap
#'@importFrom grDevices colorRampPalette
#'@importFrom RColorBrewer brewer.pal
#'@importFrom methods is
#'@import SummarizedExperiment
#'@references Liberzon, A. et al. (2015). The molecular signatures
#'database hallmark gene set collection. \emph{Cell systems}, 1, 417-425.
#'@export
#'@examples
#'# Load gene expression profiles of TNBC samples
#'data(GSE25055)
#'
#'# Predict consensus molecular subtypes of TNBC samples
#'prediction <- predictCMS(expr = GSE25055)
#'
#'# Perform GSVA on the hallmark pathway gene sets
#'resultGSVA <- performGSVA(expr = GSE25055, pred = prediction)
performGSVA <- function(expr, pred, gene.set = NULL, gsva.kcdf = "Gaussian"){

  if (is(expr, "SummarizedExperiment")){
    exp.mat <- assays(expr)[[1]]
  } else{
    exp.mat <- expr
  }

  pred <- factor(pred, levels = c("MSL", "IM", "LAR", "SL"))
  CMS_palette <- c("MSL" = "brown2", "IM" = "gold2",
                   "LAR" = "yellowgreen", "SL" = "midnightblue")
  rbpal <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)
  if(!is.null(gene.set)){
    GSVA.score <- gsva(as.matrix(exp.mat), gene.set,
                       kcdf = "Gaussian", verbose = FALSE)
    annotation_col <- data.frame(row.names =
                                   colnames(exp.mat)[order(pred)],
                                 CMS = pred[order(pred)])
    pheatmap(GSVA.score[, colnames(exp.mat)[order(pred)]],
             scale = "row", color = rbpal,
             annotation_col = annotation_col, gaps_row = c(12, 17, 24),
             annotation_colors = list(CMS = CMS_palette),
             cluster_cols = FALSE, cluster_rows = FALSE, legend = TRUE,
             annotation_names_col = FALSE, show_colnames = FALSE)
    return(GSVA.score)
  }

  GSVA.score <- gsva(as.matrix(exp.mat), Hallmark.geneset, parallel.sz = 1,
                     kcdf = "Gaussian", verbose = FALSE)
  rownames(GSVA.score) <- sapply(rownames(GSVA.score),
                                 function(x) substr(x, 10, nchar(x)))
  MSL.pathway <- c("KRAS_SIGNALING", "COAGULATION",
                   "EPITHELIAL_MESENCHYMAL_TRANSITION", "MYOGENESIS",
                   "ANGIOGENESIS", "APICAL_JUNCTION",
                   "TNFA_SIGNALING_VIA_NFKB", "APOPTOSIS",
                   "TGF_BETA_SIGNALING", "P53_PATHWAY",
                   "HYPOXIA", "HEDGEHOG_SIGNALING")
  IM.pathway <- c("INTERFERON_GAMMA_RESPONSE", "IL6_JAK_STAT3_SIGNALING",
                  "IL2_STAT5_SIGNALING", "INFLAMMATORY_RESPONSE",
                  "INTERFERON_ALPHA_RESPONSE")
  LAR.pathway <- c("ESTROGEN_RESPONSE_EARLY", "ESTROGEN_RESPONSE_LATE",
                   "ANDROGEN_RESPONSE", "BILE_ACID_METABOLISM",
                   "ADIPOGENESIS", "FATTY_ACID_METABOLISM",
                   "PROTEIN_SECRETION")
  SL.pathway <- c("G2M_CHECKPOINT", "E2F_TARGETS", "MYC_TARGETS_V2",
                  "MYC_TARGETS_V1", "MITOTIC_SPINDLE",
                  "DNA_REPAIR", "WNT_BETA_CATENIN_SIGNALING",
                  "NOTCH_SIGNALING", "GLYCOLYSIS")
  annotation_col <- data.frame(row.names = colnames(exp.mat)[order(pred)],
                               CMS = pred[order(pred)])
  pheatmap(GSVA.score[c(MSL.pathway, IM.pathway, LAR.pathway, SL.pathway),
                      colnames(exp.mat)[order(pred)]],
           scale = "row", color = rbpal,
           annotation_col = annotation_col, gaps_row = c(12, 17, 24),
           annotation_colors = list(CMS = CMS_palette),
           cluster_cols = FALSE, cluster_rows = FALSE, legend = TRUE,
           annotation_names_col = FALSE, show_colnames = FALSE)

  return(GSVA.score[c(MSL.pathway, IM.pathway, LAR.pathway, SL.pathway),])

}
