#' uniConSigData.
#'
#' @title pre-calculated penalization factors and molecular concepts coefficients
#' @name preCal.data
#' @docType data
#' @author Xu Chi <xuchi_chisquare@hotmail.com>
#' @description This file provides the coefficients for uniConSig calculation. <br> The details of the coefficients can be found in the package "uniConSig"
#' @usage data(preCal.data)
#' @keywords datasets
#' @format A list of characters. Each element contains information for a single <br>
#' gene in a tab delimited format. After splitting, the first <br>
#' element is Entrez Gene ID, the second one is a penalization factor <br>
#' used in uniConSig, the rest are the molecular concepts which are <br>
#' related to this gene and their coefficients.
#' @examples
#' library(uniConSigData)
#' data(preCal.data)
#' length(preCal.data)
#' class(preCal.data)
#' preCal.data[[1]]
NULL
