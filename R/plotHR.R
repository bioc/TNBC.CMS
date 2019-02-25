#'Forest plot of hazard ratios
#'
#'Produces a forest plot of hazard ratios for each gene. Also draws a forest plot of subtype-specific hazard ratios.
#'
#'@param exp.mat A matrix of gene expression with genes in rows and samples in columns (rownames corresopnding to gene symbol).
#'@param gene.symbol A vector of gene symbols for which hazard ratios are computed.
#'@param pred A vector of predicted consensus molecular subtypes.
#'@param time A vector of the follow-up time.
#'@param event A vector representing survival status (0 = alive, 1 = dead).
#'@param by.subtype A logical to determine if subtype-specific hazard ratios are computed (default is TRUE).
#'@return NULL
#'@importFrom survival Surv coxph
#'@importFrom R.utils insert
#'@import grid
#'@import ggplot2
#'@import forestplot
#'@export
#'@examples
#'# Load gene expression profiles and clinical information of TNBC samples
#'data(GSE25055.exprs)
#'data(GSE25055.clinical)
#'
#'# Predict consensus molecular subtypes of TNBC samples
#'predictions <- predictCMS(exp.mat = GSE25055.exprs)
#'
#'# Forest plot of hazard ratios for input genes
#'plotHR(exp.mat = GSE25055.exprs, gene.symbol = sample(rownames(GSE25055.exprssion), 4),
#'        pred = pred, time = GSE25055.clinical$DFS.month,
#'        event = GSE25055.clinical$DFS.status, by.subtype = FALSE)
#'
#'# Subtype-specific forest plot of hazard ratios for input genes
#'plotHR(exp.mat = GSE25055.exprs, gene.symbol = sample(rownames(GSE25055.exprssion), 4),
#'        pred = pred, time = GSE25055.clinical$DFS.month,
#'        event = GSE25055.clinical$DFS.status, by.subtype = TRUE)
plotHR <- function(exp.mat, gene.symbol, pred, time, event, by.subtype = TRUE){

  options(stringsAsFactors = FALSE)
  pred <- factor(pred, levels = c("MSL", "IM", "LAR", "SL"))
  CMS_palette <- c("MSL" = "brown2", "IM" = "gold2",
                   "LAR" = "yellowgreen", "SL" = "midnightblue")
  ng <- length(gene.symbol)
  survdat <- data.frame(row.names = colnames(exp.mat), pred, time, event)
  est.all <- c()
  se.all <- c()
  pcox.all <- c()
  est.ss <- c()
  se.ss <- c()
  pcox.ss <- c()
  msl <- colnames(exp.mat)[as.character(pred) == "MSL"]
  im <- colnames(exp.mat)[as.character(pred) == "IM"]
  lar <- colnames(exp.mat)[as.character(pred) == "LAR"]
  sl <- colnames(exp.mat)[as.character(pred) == "SL"]
  samples <- list(msl, im, lar, sl)
  names(samples) <- c("MSL", "IM", "LAR", "SL")
  samples <- samples[lapply(samples, length) > 0]
  ns <- length(samples)
  CMS_palette <- CMS_palette[names(samples)]

  for(x in gene.symbol){

    rg.by.gene <- ifelse(unlist(exp.mat[x,]) > median(unlist(exp.mat[x,])), "H", "L")
    fit.all <- coxph(Surv(time, event == 1) ~ unlist(exp.mat[x,]), data = survdat)
    est.all <- c(est.all, summary(fit.all)$coefficients[1])
    se.all <- c(se.all, summary(fit.all)$coefficients[3])
    pcox.all <- c(pcox.all, summary(fit.all)$coefficients[5])

    for(i in 1:ns){

      rg.by.gene <- ifelse(unlist(exp.mat[x,samples[[i]]]) > median(unlist(exp.mat[x,samples[[i]]])), "H", "L")
      fit.ss <- coxph(Surv(time, event == 1) ~ rg.by.gene, data = survdat[samples[[i]],])
      est.ss <- c(est.ss, summary(fit.ss)$coefficients[1])
      se.ss <- c(se.ss, summary(fit.ss)$coefficients[3])
      pcox.ss <- c(pcox.ss, summary(fit.ss)$coefficients[5])

    }

  }

  formatNumber <- function(x){
    res = c()
    for(num in x){
      if(abs(num) >= 1000){
        res <- c(res, format(num, scientific = TRUE, digits = 2, nsmall = 2))
      }else{
        res <- c(res, format(round(num, 2), nsmall = 2))
      }
    }
    return(res)
  }

  if(by.subtype){
    dat.ss <- data.frame(gene = rep(gene.symbol, each = ns),
                         subtype = rep(names(samples), times = ng),
                         estimate = est.ss, se = se.ss, pvalue = pcox.ss)
    dat.ss$subtype <- factor(dat.ss$subtype, levels = c("MSL", "IM", "LAR", "SL"))
    dat.ss$lower <- dat.ss$estimate - 1.96 * dat.ss$se
    dat.ss$upper <- dat.ss$estimate + 1.96 * dat.ss$se
    dat.ss$interval <- paste0(formatNumber(dat.ss$estimate), "[",
                              formatNumber(dat.ss$lower), ", ",
                              formatNumber(dat.ss$upper), "]")
    dat.ss$pvalue <- format(dat.ss$pvalue, scientific = T, digits = 2, nsmall = 2)

    tmp <- rep(NA, ng * (ns + 1))
    tmp[1] <- "Gene"
    tmp[(0:(ng - 1)) * (ns + 1) + 2] <- gene.symbol
    table.ss <- cbind(tmp, c("log HR", insert(dat.ss$interval, (1:(ng - 1)) * ns+1)),
                      c("P-value", insert(dat.ss$pvalue, (1:(ng - 1)) * ns+1)))

    input.ss <- data.frame(mean = insert(dat.ss$estimate, (0:(ng - 1)) * ns + 1),
                           lower = insert(dat.ss$lower, (0:(ng - 1)) * ns + 1),
                           upper = insert(dat.ss$upper, (0:(ng - 1)) * ns + 1))

    fn <- local({
      i <- 0
      b_clrs <- rep(CMS_palette, ng)
      l_clrs <- rep(CMS_palette, ng)

      function(..., clr.line, clr.marker){
        i <<- i + 1
        fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
      }
    })

    forestplot(table.ss, input.ss, is.summary = c(TRUE, rep(FALSE, ng * (ns + 1) - 1)),
               hrzl_lines = gpar(col = "#444444"), graph.pos = 2, clip = c(-2, 2),
               xticks = c(-2, -1, 0, 1, 2), lwd.ci = 2,
               txt_gp = fpTxtGp(ticks = gpar(cex = 0.9),
                                xlab = gpar(cex = 1, fontface = "bold")),
               fn.ci_norm = fn, xlab = "log hazard ratio")
  }else{

    dat.all <- data.frame(gene = gene.symbol, estimate = est.all, se = se.all, pvalue = pcox.all)
    dat.all$lower <- dat.all$estimate - 1.96 * dat.all$se
    dat.all$upper <- dat.all$estimate + 1.96 * dat.all$se
    dat.all$interval <- paste0(formatNumber(dat.all$estimate), "[",
                               formatNumber(dat.all$lower), ", ",
                               formatNumber(dat.all$upper), "]")
    dat.all$pvalue <- format(dat.all$pvalue, scientific = T, digits = 2, nsmall = 2)

    table.all <- cbind(c("Gene", dat.all$gene), c("log HR", dat.all$interval),
                       c("P-value", dat.all$pvalue))
    input.all <- data.frame(mean =  c(NA, dat.all$estimate), lower = c(NA, dat.all$lower),
                            upper = c(NA, dat.all$upper))

    forestplot(table.all, input.all, is.summary = c(TRUE, rep(FALSE, ng)),
               hrzl_lines = gpar(col = "#444444"), graph.pos = 2, lwd.ci = 2,
               clip = c(-2, 2), xticks = c(-2, -1, 0, 1, 2),
               txt_gp = fpTxtGp(ticks = gpar(cex = 0.9),
                                xlab = gpar(cex = 1, fontface = "bold")),
               col = fpColors(lines = "black"),
               xlab = "log hazard ratio")
  }

}
