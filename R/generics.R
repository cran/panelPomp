## generic functions

#' @title Extract parameters (coefficients) of a panel model
#' @description \code{pparams()} is a generic function that extracts parameter
#' (coefficient) values from objects returned by panel modeling functions. While
#' the named \code{numeric} vector format is useful and possible via S4 methods
#' for \code{coef()}, alternative formats capturing the panel structure can be
#' implemented via \code{pparams()}.
#' @param object an object for which extraction of panel model parameter
#' (coefficient) values is meaningful.
#' @param ... additional arguments.
#' @details This is a generic function: methods can be defined for it.
#' @return Parameter (coefficient) values extracted from the panel model
#' \code{object}.
#'
#' \pparamsReturn
#' @example examples/prw.R
#' @example examples/pparams.R
#' @keywords internal
#' @seealso \link{panelPomp_methods}
#' @author Carles \Breto
#' @export
setGeneric(name = "pparams",
           def = function(object, ...) standardGeneric("pparams"))

#' @title Extract units of a panel model
#' @description \code{unitobjects()} is a generic function that extracts a list
#' of objects corresponding to units of panel objects returned by panel modeling
#' functions.
#' @param object an object for which extraction of panel units is meaningful.
#' @param ... additional arguments.
#' @details This is a generic function: methods can be defined for it.
#' @return Units extracted from the panel model \code{object}.
#'
#' \unitobjectsReturn
#' @example examples/prw.R
#' @example examples/unitobjects.R
#' @keywords internal
#' @seealso \link{panelPomp_methods}
#' @author Carles \Breto
#' @export
setGeneric(name = "unitobjects",
           def = function(object, ...) standardGeneric("unitobjects"))

#' @title Extract log likelihood of units of panel models
#' @description \code{unitlogLik()} is a generic function that extracts the log
#' likelihood for each unit of panel objects returned by panel modeling functions.
#' While the \code{numeric} value with the log likelihood for the entire panel
#' is useful and possible via S4 methods \code{logLik()}, the contributions to it
#' by panel units can be implemented via \code{unitlogLik()}.
#' @param object an object for which log likelihood values for units can be extracted.
#' @param ... additional arguments.
#' @details This is a generic function: methods can be defined for it.
#' @return Log likelihood extracted for each unit of the panel model \code{object}.
#'
#' \unitlogLikReturn
#' @example examples/pfrw.R
#' @example examples/unitlogLik.R
#' @keywords internal
#' @seealso \link{pfilter}
#' @author Carles \Breto
#' @export
setGeneric(name = "unitlogLik",
           def = function(object, ...) standardGeneric("unitlogLik"))
