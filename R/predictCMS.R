#'TNBC consensus molecular subtype prediction
#'
#'This function predicts the TNBC consensus molecular subtype of TNBC samples.
#'
#'@param exp.mat A matrix of gene expression with genes in rows and samples in columns (rownames corresopnding to gene symbol).
#'@param gene.sig logical to determine if gene expression signature scores are computed
#'@param gsa logical to determine if gene set analysis is performed
#'@param drug.sig logical to determine if drug signatures are computed
#'@param rnaseq logical to determine if input data is RNA-Seq gene expression profile
#'@return A vector of assigned subtypes.
#'@export
#'@importFrom e1071 impute
#'@importFrom stats predict median cor
#'@importFrom pheatmap pheatmap
#'@importFrom grDevices colorRampPalette
#'@importFrom RColorBrewer brewer.pal
predictCMS <- function(exp.mat, gene.sig = T, gsa = T, drug.sig = T, rnaseq = F){

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
  pred = predict(svmModel, data.frame(t(input[genelist,]), check.names = F), probability = T)

  #Calculate signature scores
  if(gene.sig){
    EMTScore = colMeans(exp.mat[rownames(exp.mat) %in% EMTSig,])
    CINScore = colMeans(exp.mat[rownames(exp.mat) %in% CINSig,])
    HormoneScore = colMeans(exp.mat[rownames(exp.mat) %in% c("AR", "ESR1", "ERBB2", "PGR"),])
    X = exp.mat[rownames(exp.mat) %in% names(StemSig),]
    StemSig = StemSig[rownames(X)]
    s = apply(X, 2, function(z) {cor(z, StemSig, method="sp", use="complete.obs" )} )
    s = s - min(s)
    StemnessScore = s / max(s)
    ESTIMATEScore = computeESTIMATEscore(exp.mat)
    StromalScore = unlist(ESTIMATEScore["StromalSignature",])
    ImmuneScore = unlist(ESTIMATEScore["ImmuneSignature",])
    MicroenvironmentScore = computexCellScore(exp.mat, rnaseq)
    sig.mat = rbind(EMTScore, StromalScore, ImmuneScore, MicroenvironmentScore, StemnessScore,
                    HormoneScore, CINScore)
    rownames(sig.mat) = c("EMT", "Stromal", "Immune", "Microenvironment", "Stemness", "Hormone", "CIN")
    sig.mat = t(apply(sig.mat,1, function(x) tapply(x, pred, mean)))
    sig.mat <- t(scale(t(sig.mat)))
    pheatmap(sig.mat, color = colorRampPalette(rev(brewer.pal(n = 11, name =
                                                               "RdBu")))(100),
             annotation_col = data.frame(row.names = levels(pred),
                                         CMS = levels(pred)),
             annotation_colors = list(CMS = c("MSL" = "brown2", "IM" = "gold2",
                                                  "LAR" = "yellowgreen", "SL" = "midnightblue")),
             cluster_cols = F, cluster_rows = F, legend = T,
             annotation_names_col = F, show_colnames = F)
  }

  #Gene set analysis using CAMERA
  if(gsa){

    camera.mat = performCAMERA(exp.mat, pred)
    camera.mat = camera.mat[apply(camera.mat, 1, max) > -log10(0.05),]
    camera.rank = apply(-camera.mat, 2, rank)
    camera.keep = unique(unlist(apply(camera.rank, 2, function(x) names(x[x <= 5]))))
    pheatmap(camera.mat[camera.keep,], color = colorRampPalette(rev(brewer.pal(n = 11, name =
                                                                          "RdBu")))(100),
             annotation_col = data.frame(row.names = levels(pred),
                                         CMS = levels(pred)),
             annotation_colors = list(CMS = c("MSL" = "brown2", "IM" = "gold2",
                                                  "LAR" = "yellowgreen", "SL" = "midnightblue")),
             annotation_names_col = F, show_colnames = F,
             cluster_cols = F, cluster_rows = F, legend = T)

  }

  #Calculate drug signatures
  if(drug.sig){

    drugset = lapply(drugset, function(x) x[x %in% rownames(exp.mat)])
    drugmat = t(sapply(drugset, function(x) colMeans(exp.mat[x,])))
    drugdat = data.frame(row.names = colnames(exp.mat))
    drugdat$ANDROGEN_RESPONSE = drugmat["DOANE_RESPONSE_TO_ANDROGEN_UP", ]
    drugdat$IMATINIB_RESISTANCE = drugmat["MAHADEVAN_IMATINIB_RESISTANCE_UP", ]
    drugdat$CISPLATIN_RESISTANCE = colMeans(drugmat[c("LI_CISPLATIN_RESISTANCE_UP",
                                                      "TSUNODA_CISPLATIN_RESISTANCE_UP",
                                                      "KANG_CISPLATIN_RESISTANCE_UP"), ])
    drugdat$BEXAROTENE_RESPONSE = drugmat["WANG_RESPONSE_TO_BEXAROTENE_UP", ]
    drugdat$TROGLITAZONE_RESPONSE = drugmat["RUAN_RESPONSE_TO_TROGLITAZONE_UP", ]
    drugdat$DOXORUBICIN_RESISTANCE = drugmat["KANG_DOXORUBICIN_RESISTANCE_UP", ]
    drugdat$APLIDIN_RESPONSE = drugmat["MITSIADES_RESPONSE_TO_APLIDIN_UP", ]
    drugdat$DASATINIB_RESISTANCE = drugmat["HUANG_DASATINIB_RESISTANCE_UP", ]
    drugdat$TZD_RESPONSE = drugmat["GERHOLD_RESPONSE_TO_TZD_UP", ]
    drugdat$SB216763_RESPONSE = drugmat["WANG_RESPONSE_TO_GSK3_INHIBITOR_SB216763_UP", ]
    drugdat = t(drugdat)
    drugdat_mean = t(apply(drugdat,1, function(x) tapply(x, pred, mean)))
    drugdat_mean = t(scale(t(drugdat_mean)))
    drugsig = c("APLIDIN_RESPONSE", "CISPLATIN_RESISTANCE", "DASATINIB_RESISTANCE", "IMATINIB_RESISTANCE",
                "TROGLITAZONE_RESPONSE", "TZD_RESPONSE", "SB216763_RESPONSE", "ANDROGEN_RESPONSE",
                "BEXAROTENE_RESPONSE", "DOXORUBICIN_RESISTANCE")
    annotation_col = data.frame(row.names = levels(pred),
                                CMS = levels(pred))
    pheatmap(drugdat_mean[drugsig,], color = colorRampPalette(rev(brewer.pal(n = 11, name =
                                                                               "RdBu")))(100),
             annotation_col = annotation_col,
             annotation_colors = list(CMS = c("MSL" = "brown2", "IM" = "gold2",
                                              "LAR" = "yellowgreen", "SL" = "midnightblue")),
             cluster_cols = F, cluster_rows = F, legend = T,
             annotation_names_col = F, show_colnames = F)
  }

  return(pred)

}

#'Computation of stromal and immune scores
#'
#'Computes stromal and immune scores. This function is borrowed from the \code{estimate} package and modified to accept R object as input.
#'
#'@param mat A matrix of gene expression with genes in rows and samples in columns (rownames corresopnding to gene symbol).
#'@return A data frame containing stromal and immune scores
computeESTIMATEscore <- function(mat){

  m <- mat
  gene.names <- rownames(mat)
  sample.names <- colnames(mat)
  Ns <- length(m[1, ])
  Ng <- length(m[, 1])
  for (j in 1:Ns) {
    m[, j] <- rank(m[, j], ties.method = "average")
  }
  m <- 10000 * m/Ng
  gs <- as.matrix(SI_geneset[, -1], dimnames = NULL)
  N.gs <- 2
  gs.names <- row.names(SI_geneset)
  score.matrix <- matrix(0, nrow = N.gs, ncol = Ns)
  for (gs.i in 1:N.gs) {
    gene.set <- gs[gs.i, ]
    gene.overlap <- intersect(gene.set, gene.names)
    if (length(gene.overlap) == 0) {
      score.matrix[gs.i, ] <- rep(NA, Ns)
      next
    }
    else {
      ES.vector <- vector(length = Ns)
      for (S.index in 1:Ns) {
        gene.list <- order(m[, S.index], decreasing = TRUE)
        gene.set2 <- match(gene.overlap, gene.names)
        correl.vector <- m[gene.list, S.index]
        TAG <- sign(match(gene.list, gene.set2, nomatch = 0))
        no.TAG <- 1 - TAG
        N <- length(gene.list)
        Nh <- length(gene.set2)
        Nm <- N - Nh
        correl.vector <- abs(correl.vector)^0.25
        sum.correl <- sum(correl.vector[TAG == 1])
        P0 <- no.TAG/Nm
        F0 <- cumsum(P0)
        Pn <- TAG * correl.vector/sum.correl
        Fn <- cumsum(Pn)
        RES <- Fn - F0
        max.ES <- max(RES)
        min.ES <- min(RES)
        if (max.ES > -min.ES) {
          arg.ES <- which.max(RES)
        }
        else {
          arg.ES <- which.min(RES)
        }
        ES <- sum(RES)
        EnrichmentScore <- list(ES = ES, arg.ES = arg.ES,
                                RES = RES, indicator = TAG)
        ES.vector[S.index] <- EnrichmentScore$ES
      }
      score.matrix[gs.i, ] <- ES.vector
    }
  }
  score.data <- data.frame(score.matrix)
  names(score.data) <- sample.names
  row.names(score.data) <- gs.names
  score.data

}

#'Computation of microenvironment score
#'
#'Computes microenvironment score. This function wraps around the \code{xCellAnalysis} function of the \code{xCell} package to compute microenvironment score.
#'
#'@param mat A matrix of gene expression with genes in rows and samples in columns (rownames corresopnding to gene symbol).
#'@param rnaseq logical to determine if input data is RNA-Seq gene expression profile
#'@return A data frame containing stromal and immune scores
#'@importFrom GSVA gsva
#'@importFrom pracma lsqlincon
#'@importFrom stats aggregate
computexCellScore <- function(mat, rnaseq){

  signatures = xCell.data$signatures
  genes = xCell.data$genes
  if(rnaseq){
    spill = xCell.data$spill
  } else{
    spill = xCell.data$spill.array
  }

  #Compute enrichment scores
  shared.genes = intersect(rownames(mat), genes)
  expr = mat[shared.genes,]
  expr = apply(expr, 2, rank)
  scores = gsva(expr, signatures, method = "ssgsea", ssgsea.norm = F,
                      parallel.sz = 4, verbose = F)
  scores = scores - apply(scores, 1, min)
  cell_types <- unlist(strsplit(rownames(scores), "%"))
  cell_types <- cell_types[seq(1, length(cell_types), 3)]
  agg <- aggregate(scores ~ cell_types, FUN = mean)
  rownames(agg) <- agg[, 1]
  scores <- agg[, -1]

  #Transform scores
  rows <- rownames(scores)[rownames(scores) %in% rownames(spill$fv)]
  tscores <- scores[rows, ]
  minX <- apply(tscores, 1, min)
  A <- rownames(tscores)
  tscores <- (as.matrix(tscores) - minX)/5000
  tscores[tscores < 0] <- 0
  tscores <- (tscores^spill$fv[A,2])/(spill$fv[A,3]*2)

  #Adjust scores
  K = spill$K * 0.5
  diag(K) = 1
  rows <- rownames(tscores)[rownames(tscores) %in%
                              rownames(K)]
  ascores <- apply(tscores[rows, ], 2, function(x) lsqlincon(K[rows,rows],
                                                                     x, lb = 0))
  ascores[ascores<0] = 0
  rownames(ascores) <- rows

  #Compute microenvironment scores
  ImmuneScore = apply(ascores[c('B-cells','CD4+ T-cells','CD8+ T-cells','DC','Eosinophils','Macrophages','Monocytes','Mast cells','Neutrophils','NK cells'),],2,sum)/1.5
  StromaScore = apply(ascores[c('Adipocytes','Endothelial cells','Fibroblasts'),],2,sum)/2
  MicroenvironmentScore = ImmuneScore+StromaScore

  MicroenvironmentScore
}

#'Gene set analysis using CAMERA
#'
#'Performs a competitive gene set test (CAMERA) for each subtype. This function is borrowed from the \code{CMScaller} package and modified to test for TNBC consensus molecular subtype.
#'
#'@param mat A log-transformed matrix of gene expression
#'@param predicted A vector of predicted subtypes
#'@return A matrix with adjusted p-value for each subtype
#'@importFrom stats model.matrix
#'@importFrom limma makeContrasts camera
performCAMERA <- function(mat, predicted){

  camera.result = vector("list", 4)
  names(camera.result) = levels(predicted)
  design = model.matrix(~ 0 + predicted)
  colnames(design) = levels(predicted)
  index = lapply(HALLMARK_geneset, function(x) rownames(mat) %in% x)
  for(i in 1:4){
    contrasts = paste0(levels(predicted)[i], "-(", paste(levels(predicted)[-i], collapse = "+"), ")/3")
    contrast = makeContrasts(contrasts = contrasts, levels = design)
    camera.result[[i]] = camera(y = log2(mat + 1), index = index, design = design,
                                contrast = contrast, inter.gene.cor = NA, sort = F)
  }
  camera.mat = matrix(NA, nrow = length(index), ncol = 4,
                      dimnames = list(names(index), levels(predicted)))
  for(i in 1:4){
    camera.mat[, i] = -log10(camera.result[[i]]$FDR)
    isDown = camera.result[[i]]$Direction == "Down"
    camera.mat[isDown, i] = log10(camera.result[[i]]$FDR[isDown])
  }
  camera.mat[is.na(camera.mat)] = 0

  camera.mat
}
