#' An S4 class to represent the ihw output.
#'
#' @slot df  A data.frame that collects the input data, including the vector of p values and the covariate, the group assignment, as well as outputs (weighted p-values, adjusted p-values)
#' @slot weights  A (nbins X nfolds) matrix of the weight assigned to each stratum
#' @slot alpha   Numeric, the nominal significance level at which the FDR is to be controlled
#' @slot nbins  Integer, number of distinct levels into which the hypotheses were stratified
#' @slot nfolds  Integer, number of folds for pre-validation procedure
#' @slot regularization_term Numeric vector, the final value of the regularization parameter within each fold
#' @slot m_groups Integer vector, number of hypotheses tested in each stratum
#' @slot penalty  Character, "uniform deviation" or "total variation"
#' @slot covariate_type  Character, "ordinal" or "nominal"
#' @slot adjustment_type  Character, "BH" or "bonferroni"
#' @slot reg_path_information  A data.frame, information about the whole regularization path. (Currently not used, thus empty)
#' @slot solver_information  A list, solver specific output, e.g. were all subproblems solved to optimality? (Currently empty list)
#'
#' @return The different methods applied to an ihw_oldResult object can return the following:
#' @return 1) A vector
#'     of length equal to the number of hypotheses tested (e.g. the adjusted p-value or the weight
#'     of each hypothesis).
#' @return 2) A matrix of dimension equal to nfolds x nbins (e.g. the weight of
#'     each stratum, fold combination, set by specifying levels_only=TRUE).
#' @return 3) A vector of length 1 (usually a parameter of the ihw_oldResult object such as nfolds or
#'     the total number of rejections).
#' @return 4) A data.frame (as.data.frame) or just console output (show) for the extended Base generics.
#' @return See section below for the individual methods.
#' @examples
#'
#' save.seed <- .Random.seed; set.seed(1)
#' X   <- runif(n = 20000, min = 0.5, max = 4.5)       # Covariate
#' # Is the null hypothesis (mean=0) true or false ?
#' H   <- rbinom(n = length(X), size = 1, prob = 0.1)  
#' Z   <- rnorm(n = length(X), mean = H * X)           # Z-score
#' .Random.seed <- save.seed
#'
#' pvalue <- 1 - pnorm(Z)                              # pvalue
#' ihw_res <- ihw(pvalue, covariates = X, alpha = 0.1)
#' rejections(ihw_res)
#' colnames(as.data.frame(ihw_res))
#'
#' @seealso ihw
#' @import methods

ihw_oldResult <- setClass("ihw_oldResult",
         	      slots = list(
           		     df = "data.frame",
           		     weights = "matrix",
           		     alpha = "numeric",
                   nbins = "integer",
                   nfolds = "integer",
                   regularization_term = "numeric",
                   m_groups = "integer",
                   penalty = "character",
                   covariate_type = "character",
                   adjustment_type = "character",
                   reg_path_information = "data.frame",
           		     solver_information= "list"))

#-----------------------------adjusted p-values extraction---------------------------------------------------------#
adj_pvalues_ihw_oldResult <- function(object){
  object@df$adj_pvalue
}

#' @rdname ihw_oldResult-class
setGeneric("adj_pvalues", function(object) standardGeneric("adj_pvalues"))

#' @describeIn ihw_oldResult Extract adjusted pvalues
#' @export
setMethod("adj_pvalues", signature(object="ihw_oldResult"),
          adj_pvalues_ihw_oldResult)

#--------------------------- weights extraction --------------------------------------------------------------------#


weights_ihw_oldResult <-function(object, levels_only = FALSE){
  if (levels_only) {
    object@weights
  } else {
    object@df$weight #TODO: Storing redundant information right now
  }
}

#' @param object,x A ihw_oldResult object as returned by a call to ihw(...)
#' @param levels_only Logical, if FALSE, return a vector of weights (thresholds) with one weight
#'    (threshold) for each hypothesis, otherwise return a nfolds x nbins matrix of weights (thresholds)
#' @param row.names,optional See ?base::as.data.frame for a description of these arguments.
#' @param ... Parameters passed in to individual methods
#' @describeIn ihw_oldResult Extract weights
#' @export
setMethod("weights", signature(object="ihw_oldResult"),
          weights_ihw_oldResult)


#--------------------------- threshold extraction -------------------------------------------------------------------#


thresholds_ihw_oldResult <-function(object, levels_only = FALSE){
  t <- get_bh_threshold(na.exclude(weighted_pvalues(object)), alpha(object),
                        mtests = sum(as.numeric(m_groups(object))))
  t*weights(object, levels_only = levels_only)
}

#' @rdname ihw_oldResult-class
setGeneric("thresholds", function(object,...) standardGeneric("thresholds"))

#' @describeIn ihw_oldResult Calculate ihw thresholds
#' @export
setMethod("thresholds", signature(object="ihw_oldResult"),
          thresholds_ihw_oldResult)



#--------------------------- p-value extraction ---------------------------------------------------------------------#

pvalues_ihw_oldResult <- function(object){
  object@df$pvalue
}

#' @rdname ihw_oldResult-class
setGeneric("pvalues", function(object) standardGeneric("pvalues"))

#' @describeIn ihw_oldResult Extract pvalues
#' @export
setMethod("pvalues", signature(object="ihw_oldResult"),
          pvalues_ihw_oldResult)


#--------------------------- weighted p-value extraction -------------------------------------------------------------#

weighted_pvalues_ihw_oldResult <- function(object){
  object@df$weighted_pvalue
}

#' @rdname ihw_oldResult-class
setGeneric("weighted_pvalues", function(object) standardGeneric("weighted_pvalues"))

#' @describeIn ihw_oldResult Extract weighted pvalues
#' @export
setMethod("weighted_pvalues", signature(object="ihw_oldResult"),
          weighted_pvalues_ihw_oldResult)


#---------------------------  covariate extraction ----------------------------------------------------------#

covariates_ihw_oldResult <- function(object){
  object@df$covariate
}

#' @rdname ihw_oldResult-class
setGeneric("covariates", function(object) standardGeneric("covariates"))

#' @describeIn ihw_oldResult Extract covariates
#' @export
setMethod("covariates", signature(object="ihw_oldResult"),
          covariates_ihw_oldResult)

#---------------------------- covariate type -----------------------------------------------------------------#
covariate_type_ihw_oldResult <- function(object){
  object@covariate_type
}

#' @rdname ihw_oldResult-class
setGeneric("covariate_type", function(object) standardGeneric("covariate_type"))

#' @describeIn ihw_oldResult Extract type of covariate ("ordinal" or "nominal")
#' @export
setMethod("covariate_type", signature(object="ihw_oldResult"),
          covariate_type_ihw_oldResult)

#----------------- extract stratification variable------------------------------------------------------------#
groups_factor_ihw_oldResult <- function(object){
	object@df$group
}

#' @rdname ihw_oldResult-class
setGeneric("groups_factor", function(object) standardGeneric("groups_factor"))


#' @describeIn ihw_oldResult Extract factor of stratification (grouping) variable
#' @export
setMethod("groups_factor", signature(object="ihw_oldResult"),
          groups_factor_ihw_oldResult)

#----------------- nfolds extraction ----------------------------------------------------------------------#
nfolds_ihw_oldResult <- function(object) object@nfolds

#' @rdname ihw_oldResult-class
setGeneric("nfolds", function(object) standardGeneric("nfolds"))

#' @describeIn ihw_oldResult Extract number of folds
#' @export
setMethod("nfolds", signature(object="ihw_oldResult"),
          nfolds_ihw_oldResult)


#----------------- nbins extraction ----------------------------------------------------------------------#
nbins_ihw_oldResult <- function(object) object@nbins

#' @rdname ihw_oldResult-class
setGeneric("nbins", function(object) standardGeneric("nbins"))

#' @describeIn ihw_oldResult Extract number of bins
#' @export
setMethod("nbins", signature(object="ihw_oldResult"),
          nbins_ihw_oldResult)


#----------------- nominal alpha extraction ----------------------------------------------------------------------#
alpha_ihw_oldResult <-function(object) object@alpha

#' @rdname ihw_oldResult-class
setGeneric("alpha", function(object) standardGeneric("alpha"))

#' @describeIn ihw_oldResult Extract nominal significance (alpha) level
#' @export
setMethod("alpha", signature(object="ihw_oldResult"),
          alpha_ihw_oldResult)

#----------------- adjustment type extraction ----------------------------------------------------------------------#
adjustment_type_ihw_oldResult <-function(object) object@adjustment_type

setGeneric("adjustment_type", function(object) standardGeneric("adjustment_type"))

setMethod("adjustment_type", signature(object="ihw_oldResult"),
          adjustment_type_ihw_oldResult)


#----------------- rejections ------------------------------------------------------------------------------------#


# TODO: Extend this to groupwise calculation (i.e. rejections per stratum)
rejections_ihw_oldResult <- function(object){
  sum(rejected_hypotheses(object), na.rm=TRUE)
}

#' @rdname ihw_oldResult-class
setGeneric("rejections", function(object,...) standardGeneric("rejections"))

#' @describeIn ihw_oldResult Total number of rejected hypotheses by ihw procedure
#' @export
setMethod("rejections", signature(object="ihw_oldResult"),
          rejections_ihw_oldResult)



rejected_hypotheses_ihw_oldResult <- function(object){
  adj_pvalues(object) <= alpha(object)
}

#' @rdname ihw_oldResult-class
setGeneric("rejected_hypotheses", function(object,...) standardGeneric("rejected_hypotheses"))

#' @describeIn ihw_oldResult Get a boolean vector of the rejected hypotheses
#' @export
setMethod("rejected_hypotheses", signature(object="ihw_oldResult"),
          rejected_hypotheses_ihw_oldResult)


#---------------- regularization terms ---------------------------
regularization_term_ihw_oldResult <-function(object) object@regularization_term

#' @rdname ihw_oldResult-class
setGeneric("regularization_term", function(object) standardGeneric("regularization_term"))

#' @describeIn ihw_oldResult Extract vector of regularization parameters used for each stratum
#' @export
setMethod("regularization_term", signature(object="ihw_oldResult"),
          regularization_term_ihw_oldResult)

#---------------- m_groups --------------------------------------------------
m_groups_ihw_oldResult <-function(object) object@m_groups

#' @rdname ihw_oldResult-class
setGeneric("m_groups", function(object) standardGeneric("m_groups"))

#' @describeIn ihw_oldResult Extract total number of hypotheses within each stratum
#' @export
setMethod("m_groups", signature(object="ihw_oldResult"),
          m_groups_ihw_oldResult)

#--------------- convenience methods ------------------------------------------------------------------------------#

#' @rdname ihw_oldResult-class
as.data.frame_ihw_oldResult <-function(x,row.names=NULL, optional=FALSE, ...){
        as.data.frame(x@df, row.names=row.names, optional=optional)
      }

#' @describeIn ihw_oldResult Coerce ihw_oldResult to data frame
#' @importFrom BiocGenerics as.data.frame
#' @export
setMethod("as.data.frame", "ihw_oldResult",as.data.frame_ihw_oldResult)

#' @describeIn ihw_oldResult Return number of p-values
#' @export
setMethod("nrow", "ihw_oldResult", function(x) nrow(x@df))

#' @describeIn ihw_oldResult Convenience method to show ihw_oldResult object
#' @importFrom methods show
#' @export
setMethod("show", signature(object="ihw_oldResult"), function(object) {
  adj_type <- adjustment_type(object)
  if (adj_type == "BH"){
    typeI_error <- "FDR"
  } else if (adj_type == "bonferroni"){
    typeI_error <- "FWER"
  } else {
    stop("Adjustment method appears to be invalid, corrupted IHWold object.")
  }
  cat("ihw_oldResult object with", nrow(object),"hypothesis tests \n")
  cat("Nominal",  typeI_error, "control level:", alpha(object),"\n")
  cat("Split into", nbins(object),"bins, based on an", covariate_type(object), "covariate\n")
})



#------------------ not exported stuff ----------------------------------------------------------------------------#

##### FDR estimate #############################################################
plugin_fdr.ihw_oldResult <- function(object) {
  ts <- thresholds(object)
  m_groups <- table(groups_factor(object))
  sum(ts*m_groups)/rejections(object, method="thresholds")
}

setGeneric("plugin_fdr", function(object,...) standardGeneric("plugin_fdr"))


setMethod("plugin_fdr", signature(object="ihw_oldResult"),
          plugin_fdr.ihw_oldResult)


##### #############################################################
stratification_breaks.ihw_oldResult <- function(object) {
  ts <- thresholds(object)
  groups <- groups_factor(object)
  filterstat_list <- split(covariates(object), groups)
  filterstats <- sapply(filterstat_list, max)
  filterstats
}

setGeneric("stratification_breaks", function(object,...) standardGeneric("stratification_breaks"))


setMethod("stratification_breaks", signature(object="ihw_oldResult"),
          stratification_breaks.ihw_oldResult)


######## temporary: number of pvals in each stratum #############################
stratum_sizes <- function(object) table(groups_factor(object))

############# validity ##########################################################
setValidity( "ihw_oldResult", function( object ) {
	return(TRUE)
} )


per_bin_fdrs <- function(obj) {
  ts <- thresholds(obj)
  groups <- groups_factor(obj)
  pvals <- pvalues(obj)
  pv_list <- split(pvals, groups)
  print(sapply(pv_list,length))
  fdrs <- mapply(function(t,pvec) length(pvec)*t/max(1, sum(pvec <= t)), ts, pv_list)
  return(fdrs)
}


