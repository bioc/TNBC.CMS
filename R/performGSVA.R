#'Gene set variation analysis
#'
#'Performs GSVA on the hallmark gene sets. Also draws a heatmap representing GSVA scores.
#'
#'@param exp.mat A matrix of gene expression with genes in rows and samples in columns (rownames corresopnding to gene symbol).
#'@param pred A vector of predicted consensus molecular subtypes.
#'@return A matrix of GSVA enrichment scores.
#'@importFrom GSVA gsva
#'@importFrom pheatmap pheatmap
#'@importFrom grDevices colorRampPalette
#'@importFrom RColorBrewer brewer.pal
#'@references Liberzon,A. et al. (2015) The molecular signatures database hallmark gene set collection. \emph{Cell systems},1,417-425.
#'@export
performGSVA = function(exp.mat, pred){

  GSVA.score = gsva(as.matrix(exp.mat), Hallmark.geneset, kcdf = "Gaussian", verbose = F)
  rownames(GSVA.score) = sapply(rownames(GSVA.score), function(x) substr(x, 10, nchar(x)))
  MSL.pathway = c("KRAS_SIGNALING", "COAGULATION", "EPITHELIAL_MESENCHYMAL_TRANSITION",
                  "MYOGENESIS", "ANGIOGENESIS", "APICAL_JUNCTION", "TNFA_SIGNALING_VIA_NFKB",
                  "APOPTOSIS", "TGF_BETA_SIGNALING", "P53_PATHWAY", "HYPOXIA", "HEDGEHOG_SIGNALING")
  IM.pathway = c("INTERFERON_GAMMA_RESPONSE", "IL6_JAK_STAT3_SIGNALING", "IL2_STAT5_SIGNALING",
                 "INFLAMMATORY_RESPONSE", "INTERFERON_ALPHA_RESPONSE")
  LAR.pathway = c("ESTROGEN_RESPONSE_EARLY", "ESTROGEN_RESPONSE_LATE", "ANDROGEN_RESPONSE",
                  "BILE_ACID_METABOLISM", "ADIPOGENESIS", "FATTY_ACID_METABOLISM",
                  "PROTEIN_SECRETION")
  SL.pathway = c("G2M_CHECKPOINT", "E2F_TARGETS", "MYC_TARGETS_V2", "MYC_TARGETS_V1", "MITOTIC_SPINDLE",
                 "DNA_REPAIR", "WNT_BETA_CATENIN_SIGNALING", "NOTCH_SIGNALING", "GLYCOLYSIS")
  annotation_col = data.frame(row.names = colnames(exp.mat)[order(pred)],
                              CMS = pred[order(pred)])
  pheatmap(GSVA.score[c(MSL.pathway, IM.pathway, LAR.pathway, SL.pathway), colnames(exp.mat)[order(pred)]],
           scale = "row", color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
           annotation_col = annotation_col, gaps_row = c(12, 17, 24),
           annotation_colors = list(CMS = c("MSL" = "brown2", "IM" = "gold2",
                                            "LAR" = "yellowgreen", "SL" = "midnightblue")),
           cluster_cols = F, cluster_rows = F, legend = T,
           annotation_names_col = F, show_colnames = F)

  return(GSVA.score[c(MSL.pathway, IM.pathway, LAR.pathway, SL.pathway),])

}
