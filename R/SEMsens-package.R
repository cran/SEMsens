#' A Tool for Sensitivity Analysis in Structural Equation Modeling
#'
#' This package is to help researchers perform and report sensitivity analysis
#'     in structural equation modeling using a phantom variable approach
#'     proposed by (Harring, McNeish, & Hancock, 2017).
#'
#' The package covers sensitivity analysis using ant colony optimization and
#'     other meta-heuristic optimization methods (in development) to automatically
#'     search a phantom variable, if there is any, that meets the optimization function.
#'     The current package includes three main functions and they are
#'     \code{\link{gen.sens.pars}} function that generates sensitivity parameters,
#'     \code{\link{sa.aco}} function that performs sensitivity analysis, and
#'     \code{\link{sens.tables}} function that summarizes sensitivity analysis
#'     results.
#'
#' @author Walter Leite, Zuchao Shen
#'
#' Maintainer: Walter Leite \href{mailto: walter.leite@coe.ufl.edu}{walter.leite@coe.ufl.edu}
#'     (University of Florida)
#'
#' @docType package
"_PACKAGE"
