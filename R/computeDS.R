#'Computation of drug signature scores.
#'
#'Computes drug signature scores. Also draws heatmap representing the average signature scores for each subtype.
#'
#'@param exp.mat A matrix of gene expression with genes in rows and samples in columns (rownames corresopnding to gene symbol).
#'@param pred A vector of predicted consensus molecular subtypes.
#'@return A matrix of drug signature scores.
#'@details Drug signature scores are the average of expression values of genes included in gene sets from MSigDB.
#'@references Liberzon,A. et al. (2011) Molecular signatures database (MSigDB) 3.0. \emph{Bioinformatics},27,1739-40.
#'@importFrom pheatmap pheatmap
#'@importFrom grDevices colorRampPalette
#'@importFrom RColorBrewer brewer.pal
#'@export
computeDS <- function(exp.mat, pred){

  pred = factor(pred, levels = c("MSL", "IM", "LAR", "SL"))

  computeSignature = function(ds){

    gs.resistant = ds[grep("RESISTANCE", ds)];gs.responsive = ds[grep("RESPONSE", ds)]
    gs.up = c(gs.resistant[grep("_DN", gs.resistant)], gs.responsive[grep("_UP", gs.responsive)])
    gs.down = c(gs.resistant[grep("_UP", gs.resistant)], gs.responsive[grep("_DN", gs.responsive)])

    gene.up = unique(unlist(CGP.geneset[gs.up]));gene.down = unique(unlist(CGP.geneset[gs.down]))
    score.up = colMeans(exp.mat[gene.up,], na.rm = T)
    score.down = colMeans(exp.mat[gene.down,], na.rm = T)
    na.up = sum(is.na(score.up)) > 0
    na.down = sum(is.na(score.down)) > 0
    if(na.up & na.down){
      out = rep(NA, ncol(exp.mat))
      names(out) = colnames(exp.mat)
      return(out)
    }
    else if(na.up){
      return(-score.down)
    }
    else if(na.down){
      return(score.up)
    }
    else{
      return(score.up - score.down)
    }
  }

  ds.dat = data.frame(row.names = colnames(exp.mat), stringsAsFactors = F)
  ds.dat$APLIDIN = computeSignature(c("MITSIADES_RESPONSE_TO_APLIDIN_UP",
                                      "MITSIADES_RESPONSE_TO_APLIDIN_DN"))
  ds.dat$CISPLATIN = computeSignature(c("KANG_CISPLATIN_RESISTANCE_DN", "LI_CISPLATIN_RESISTANCE_UP",
                                        "TSUNODA_CISPLATIN_RESISTANCE_UP",
                                        "KANG_CISPLATIN_RESISTANCE_UP"))
  ds.dat$DASATINIB = computeSignature(c("HUANG_DASATINIB_RESISTANCE_DN",
                                        "HUANG_DASATINIB_RESISTANCE_UP"))
  ds.dat$FORSKOLIN = computeSignature("SASSON_RESPONSE_TO_FORSKOLIN_UP")
  ds.dat$IMATINIB = computeSignature("MAHADEVAN_IMATINIB_RESISTANCE_UP")
  ds.dat$TROGLITAZONE = computeSignature("RUAN_RESPONSE_TO_TROGLITAZONE_UP")
  ds.dat$TZD = computeSignature("GERHOLD_RESPONSE_TO_TZD_UP")
  ds.dat$SB216763 = computeSignature("WANG_RESPONSE_TO_GSK3_INHIBITOR_SB216763_UP")
  ds.dat$ANDROGEN = computeSignature(c("DOANE_RESPONSE_TO_ANDROGEN_UP", "WANG_RESPONSE_TO_ANDROGEN_UP",
                                       "MOTAMED_RESPONSE_TO_ANDROGEN_UP",
                                       "NELSON_RESPONSE_TO_ANDROGEN_UP",
                                       "MOTAMED_RESPONSE_TO_ANDROGEN_DN"))
  ds.dat$TAMOXIFEN = computeSignature(c("BECKER_TAMOXIFEN_RESISTANCE_UP",
                                        "RIGGINS_TAMOXIFEN_RESISTANCE_UP",
                                        "MASSARWEH_TAMOXIFEN_RESISTANCE_UP"))
  ds.dat$BEXAROTENE = computeSignature(c("WANG_RESPONSE_TO_BEXAROTENE_UP",
                                         "WANG_RESPONSE_TO_BEXAROTENE_DN"))
  ds.dat$DOXORUBICIN = computeSignature("KANG_DOXORUBICIN_RESISTANCE_UP")

  ds.dat = t(ds.dat)
  annotation_col = data.frame(row.names = colnames(exp.mat)[order(pred)],
                              CMS = pred[order(pred)])
  pheatmap(ds.dat[,colnames(exp.mat)[order(pred)]],
           color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
           annotation_col = annotation_col, scale = "row",
           annotation_colors = list(CMS = c("MSL" = "brown2", "IM" = "gold2",
                                            "LAR" = "yellowgreen", "SL" = "midnightblue")),
           cluster_cols = F, cluster_rows = F, legend = T,
           annotation_names_col = F, show_colnames = F)

  return(ds.dat)

}
