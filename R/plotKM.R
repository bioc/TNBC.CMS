#'Subtype-specific survival curves
#'
#'Produces Kaplan-Meier survival curves for each subtype.
#'
#'@param pred A vector of predicted consensus molecular subtypes.
#'@param time A vector of the follow-up time.
#'@param event A vector representing survival status (0 = alive, 1 = dead).
#'@return A \code{ggplot} object.
#'@importFrom GGally ggsurv
#'@importFrom survival Surv survfit
#'@export
#'@examples
#'# Load clinical information of TNBC samples
#'data(GSE25055)
#'DFS.status <- colData(GSE25055)$DFS.status
#'DFS.month <- colData(GSE25055)$DFS.month
#'
#'# Predict consensus molecular subtypes of TNBC samples
#'prediction <- predictCMS(expr = GSE25055)
#'
#'# Plot Kaplan-Meier curves for each subtype
#'plotKM(pred = prediction, time = DFS.month, event = DFS.status)
plotKM <- function(pred, time, event){

  pred <- factor(pred, levels = c("MSL", "IM", "LAR", "SL"))
  CMS_palette <- c("MSL" = "brown2", "IM" = "gold2",
                   "LAR" = "yellowgreen", "SL" = "midnightblue")
  CMS_palette <- CMS_palette[unique(pred)]
  surv.dat <- data.frame(CMS = pred, time = time, event = event)
  km.CMS <- survfit(Surv(time, event == 1) ~ CMS, data = surv.dat)
  ggsurv(km.CMS, surv.col = CMS_palette,
         size.est = 1, order.legend = FALSE, cens.size = 0) +
    ylab("Eventless probability") + xlab("Follow-up time") +
    theme_classic()

}
