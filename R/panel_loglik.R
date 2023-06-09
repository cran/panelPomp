#' @include get_col_row.R
NULL

## logLik,matrix-method
#' @title Handling of loglikelihood replicates
#' @description Handling of loglikelihood replicates.
#' @details When \code{se = TRUE}, the jackknife se's from
#' \code{pomp::logmeanexp} are squared, summed and the squared root is taken.
#' @param object Matrix with the same number of replicated estimates for each
#' panel unit loglikelihood.
#' @param repMargin The margin of the matrix having the replicates (1 for rows,
#'  2 for columns).
#' @param first Wether to \code{"aver"}(age replicates) or \code{"aggr"}(egate
#' units) before performing the other action.
#' @param aver How to average: \code{'logmeanexp'} to average on the likelihood
#' scale before taking logs or \code{'mean'} to average after taking logs (in
#' which case, which action is performed first does not change the result).
#' @param se logical; whether to give standard errors.
#' @return
#' \code{numeric} vector with the average panel log likelihood and, when
#' \code{se = TRUE}, the corresponding standard error.
#' @name panel_loglik
#' @aliases logLik,matrix-method
#' @family panelPomp workhorse functions
#' @author Carles \Breto
#' @examples
#' ulls <- matrix(c(1,1.1,10.1,10),nr=2)
#' # when combining log likelihood estimates, the order in which aggregation and
#' # averaging are done can make a difference: panel_logmeanexp() implements the best
#' logLik(ulls,repMargin=1,first="aver",aver="logmeanexp")
#' logLik(ulls,repMargin=1,first="aggr",aver="mean",se=TRUE)
#' @export
setMethod(
  "logLik",
  signature=signature(object="matrix"),
  definition=function (object, repMargin, first = "aver",
                       aver = "logmeanexp", se = FALSE) {
    if (first=="aver" & aver=="logmeanexp") {
      logmeanexps <- apply(object,MARGIN=ifelse(repMargin==1,2,1),
                          FUN=pomp::logmeanexp,se=TRUE)
      summed <- sum(logmeanexps[1,])
      ll <- summed
      ll.se <- sqrt(sum(unname(logmeanexps[2,])^2))
    }
    if (first=="aggr" & aver=="logmeanexp") {
      logmeanexps <- logmeanexp(apply(object,MARGIN=repMargin,FUN=sum),se=TRUE)
      ll <- unname(logmeanexps[1])
      ll.se <- unname(logmeanexps[2])
    }
    if (aver=="mean") {
      means <- apply(object,MARGIN=ifelse(repMargin==1,2,1),FUN=mean)
      vars <- apply(object,MARGIN=ifelse(repMargin==1,2,1),FUN=var)
      ll <- sum(means)
      ll.se <- sqrt(sum(vars))
    }
    if (!se) ll else c(ll,se=ll.se)
  }
)
