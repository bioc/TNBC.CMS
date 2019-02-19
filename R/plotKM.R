#'Subtype-specific survival curves
#'
#'Produces Kaplan-Meier survival curves for each subtype.
#'
#'@param pred A vector of predicted consensus molecular subtypes.
#'@param time A vector of the follow-up time.
#'@param event A vector representing survival status (0 = alive, 1 = dead).
#'@importFrom GGally ggsurv
#'@importFrom survival Surv survfit
#'@export
plotKM = function(pred, time, event){

  surv.dat = data.frame(CMS = pred, time = time, event = event)
  km.CMS = survfit(Surv(time, event == 1) ~ CMS, data = surv.dat)
  ggsurv(km.CMS, surv.col = c("brown2", "gold2", "yellowgreen", "midnightblue"),
         size.est = 1, order.legend = F, cens.size = 0) +
    ylab("Eventless probability") + xlab("Follow-up time") +
    theme_classic()

}
