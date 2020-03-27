#'Computation of drug signature scores.
#'
#'Computes drug signature scores. Also draws heatmap representing
#'the average signature scores for each subtype.
#'
#'@param expr A \code{SummarizedExperiment} object or a matrix containig gene
#'expression profiles. If input is a \code{SummarizedExperiment}, the first
#'element in the assays list should be a matrix of gene expression.
#'Rows and columns of the gene expression matrix correspond to genes and
#'samples, respectively (rownames must be to gene symbols).
#'@param pred A vector of predicted consensus molecular subtypes.
#'@param gene.set A user-provided list of gene sets associated with
#'drug response. Names of gene sets must follow the format of
#'[DRUG NAME]_[RESISTANCE/RESPONSE]_[UP/DN] (e.g. CISPLATIN_RESISTANCE_DN).
#'@return A matrix of drug signature scores.
#'@details Drug signature scores are the average of expression values of
#'genes included in gene sets from MSigDB.
#'@references Liberzon, A. et al. (2011). Molecular signatures database
#'(MSigDB) 3.0. \emph{Bioinformatics}, 27, 1739-40.
#'@importFrom pheatmap pheatmap
#'@importFrom grDevices colorRampPalette
#'@importFrom RColorBrewer brewer.pal
#'@importFrom methods is
#'@import SummarizedExperiment
#'@export
#'@examples
#'# Load gene expression profiles of TNBC samples
#'data(GSE25055)
#'
#'# Predict consensus molecular subtypes of TNBC samples
#'prediction <- predictCMS(expr = GSE25055)
#'
#'# Compute drug signature scores
#'resultDS <- computeDS(expr = GSE25055, pred = prediction)
computeDS <- function(expr, pred, gene.set = NULL){

  if (is(expr, "SummarizedExperiment")){
    exp.mat <- assays(expr)[[1]]
  } else{
    exp.mat <- expr
  }

  pred <- factor(pred, levels = c("MSL", "IM", "LAR", "SL"))
  rbpal <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)
  CMS_palette <- c("MSL" = "brown2", "IM" = "gold2",
                   "LAR" = "yellowgreen", "SL" = "midnightblue")

  exp.mat <- data.frame(exp.mat, stringsAsFactors = FALSE, check.names = FALSE)
  computeSignature <- function(ds, gs){

    gs.resistant <- ds[grep("RESISTANCE", ds)]
    gs.responsive = ds[grep("RESPONSE", ds)]
    gs.up <- c(gs.resistant[grep("_DN", gs.resistant)],
               gs.responsive[grep("_UP", gs.responsive)])
    gs.down <- c(gs.resistant[grep("_UP", gs.resistant)],
                 gs.responsive[grep("_DN", gs.responsive)])

    gene.up <- unique(unlist(gs[gs.up]))
    gene.down = unique(unlist(gs[gs.down]))
    score.up <- colMeans(exp.mat[gene.up,], na.rm = TRUE)
    score.down <- colMeans(exp.mat[gene.down,], na.rm = TRUE)
    na.up <- sum(is.na(score.up)) > 0
    na.down <- sum(is.na(score.down)) > 0
    if(na.up & na.down){
      out <- rep(NA, ncol(exp.mat))
      names(out) <- colnames(exp.mat)
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

  if(!is.null(gene.set)){
    gn <- names(gene.set)
    drug.names <- unlist(strsplit(gn, split = "\\_"))
    drug.names <- unique(drug.names[seq(0, length(gn) -1 ) * 3 + 1])
    ds.dat <- t(vapply(drug.names,
                       function(x) computeSignature(gn[grep(x, gn)],
                                                    gene.set),
                       rep(0, ncol(exp.mat))))


    annotation_col <- data.frame(row.names = colnames(exp.mat)[order(pred)],
                                 CMS = pred[order(pred)])

    pheatmap(ds.dat[,colnames(exp.mat)[order(pred)]], color = rbpal,
             annotation_col = annotation_col, scale = "row",
             annotation_colors = list(CMS = CMS_palette),
             cluster_cols = FALSE, cluster_rows = FALSE, legend = TRUE,
             annotation_names_col = FALSE, show_colnames = FALSE)

    return(ds.dat)

  }

  ds.dat <- data.frame(row.names = colnames(exp.mat),
                       stringsAsFactors = FALSE)
  aplidin.gs <- c("MITSIADES_RESPONSE_TO_APLIDIN_UP",
                  "MITSIADES_RESPONSE_TO_APLIDIN_DN")
  cisplatin.gs <- c("KANG_CISPLATIN_RESISTANCE_DN",
                    "LI_CISPLATIN_RESISTANCE_UP",
                    "TSUNODA_CISPLATIN_RESISTANCE_UP",
                    "KANG_CISPLATIN_RESISTANCE_UP")
  dasatinib.gs <- c("HUANG_DASATINIB_RESISTANCE_DN",
                    "HUANG_DASATINIB_RESISTANCE_UP")
  forskolin.gs <- "SASSON_RESPONSE_TO_FORSKOLIN_UP"
  imatinib.gs <- "MAHADEVAN_IMATINIB_RESISTANCE_UP"
  troglitazone.gs <- "RUAN_RESPONSE_TO_TROGLITAZONE_UP"
  tzd.gs <- "GERHOLD_RESPONSE_TO_TZD_UP"
  gsk.gs <- "WANG_RESPONSE_TO_GSK3_INHIBITOR_SB216763_UP"
  androgen.gs <- c("DOANE_RESPONSE_TO_ANDROGEN_UP",
                   "WANG_RESPONSE_TO_ANDROGEN_UP",
                   "MOTAMED_RESPONSE_TO_ANDROGEN_UP",
                   "NELSON_RESPONSE_TO_ANDROGEN_UP",
                   "MOTAMED_RESPONSE_TO_ANDROGEN_DN")
  tamoxifen.gs <- c("BECKER_TAMOXIFEN_RESISTANCE_UP",
                    "RIGGINS_TAMOXIFEN_RESISTANCE_UP",
                    "MASSARWEH_TAMOXIFEN_RESISTANCE_UP")
  bexarotene.gs <- c("WANG_RESPONSE_TO_BEXAROTENE_UP",
                     "WANG_RESPONSE_TO_BEXAROTENE_DN")
  doxorubicin.gs <- "KANG_DOXORUBICIN_RESISTANCE_UP"
  ds.dat$APLIDIN <- computeSignature(aplidin.gs, CGP.geneset)
  ds.dat$CISPLATIN <- computeSignature(cisplatin.gs, CGP.geneset)
  ds.dat$DASATINIB <- computeSignature(dasatinib.gs, CGP.geneset)
  ds.dat$FORSKOLIN <- computeSignature(forskolin.gs, CGP.geneset)
  ds.dat$IMATINIB <- computeSignature(imatinib.gs, CGP.geneset)
  ds.dat$TROGLITAZONE <- computeSignature(troglitazone.gs, CGP.geneset)
  ds.dat$TZD <- computeSignature(tzd.gs, CGP.geneset)
  ds.dat$SB216763 <- computeSignature(gsk.gs, CGP.geneset)
  ds.dat$ANDROGEN <- computeSignature(androgen.gs, CGP.geneset)
  ds.dat$TAMOXIFEN <- computeSignature(tamoxifen.gs, CGP.geneset)
  ds.dat$BEXAROTENE <- computeSignature(bexarotene.gs, CGP.geneset)
  ds.dat$DOXORUBICIN <- computeSignature(doxorubicin.gs, CGP.geneset)

  ds.dat <- t(ds.dat)

  annotation_col <- data.frame(row.names = colnames(exp.mat)[order(pred)],
                               CMS = pred[order(pred)])

  pheatmap(ds.dat[,colnames(exp.mat)[order(pred)]], color = rbpal,
           annotation_col = annotation_col, scale = "row",
           annotation_colors = list(CMS = CMS_palette),
           cluster_cols = FALSE, cluster_rows = FALSE, legend = TRUE,
           annotation_names_col = FALSE, show_colnames = FALSE)

  return(ds.dat)

}
