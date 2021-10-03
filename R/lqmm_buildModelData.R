#' Build model data
#'
#' Helper function for [lqmm] and the predict for allowing custom predictions
#'
#' @param modelCall The [lqmm]-call
#' @param data The `data.frame` that the object has been called with
#' @param extraModelFrameArgs A `list` with the extra arguments to [stats::model.frame]
#' @param existingModelDF If the data is to have the same structure as the previous
#'  we need to provide the original data with corresponding `factor` variables so that
#'  we don't drop factors that should have been retained.
#'
#' @return A list with `dataMix`, `mmf`, `mmf_df`, `group`, `ngroups` and some attributes
#'  used in the [lqmm] fit.
buildModelData <- function(modelCall, data, extraModelFrameArgs = list(), existingModelDF) {
  stopifnot(is.list(extraModelFrameArgs))

  # Handles nested environments
  fixed <- retrieveExpression(modelCall$fixed)
  random <- retrieveExpression(modelCall$random)

  groupFormula <- asOneSidedFormula(modelCall$group)

  # Build basic input data.frame
  mfArgs <- c(list(formula = asOneFormula(random, fixed, groupFormula[[2]]),
                   data = retrieveExpression(data),
                   na.action = modelCall$na.action),
              extraModelFrameArgs)
  dataMix <- do.call("model.frame", mfArgs, envir = parent.frame())

  origOrder <- row.names(dataMix)
  if (!is.null(modelCall$contrasts)) {
    stop("Contrasts are not yet implemented")
  }

  # Reorder into groups as assigned by the group formula
  group <- model.frame(groupFormula, dataMix)
  ord <- order(unlist(group, use.names = FALSE))
  group <- group[ord, , drop = TRUE]
  ngroups <- length(unique(group))

  dataMix <- dataMix[ord, , drop = FALSE]
  revOrder <- match(origOrder, row.names(dataMix))

  # Prep mmf
  mmf_df = model.frame(fixed, dataMix)

  # Fix factors so that they match the original model's factors
  if (!missing(existingModelDF)) {
    stopifnot(is.data.frame(existingModelDF))
    factorVariables <- names(existingModelDF)[sapply(existingModelDF, is.factor)]
    for (varname in factorVariables) {
      newDataLevels <- levels(newdata[[varname]])
      orgDataLevels <- levels(existingModelDF[[varname]])
      existingLevels <- newDataLevels %in% orgDataLevels
      if (any(!existingLevels)) {
        stop("The variable ", varname, " has levels not present in original model data: ",
             paste(newDataLevels[existingLevels], collapse = ", "))
      }

      if (length(orgDataLevels) > length(newDataLevels)) {
        newdata[[varname]] <- factor(newdata[[varname]], levels = orgDataLevels)
      }
    }
  }
  mmf <- model.matrix(as.formula(fixed), data = mmf_df)

  structure(list(dataMix = dataMix,
                 mmf = mmf,
                 mmf_df = mmf_df,
                 group = group,
                 ngroups = ngroups),
            class = "lqmm.processed_input_data",
            # Save variables to append to the fit object
            # These variables serve no other purpose in the fit
            mfArgs = mfArgs,
            revOrder = revOrder,
            origOrder = origOrder)
}

#' Build the mixed effects model
#'
#' @inheritParams buildModelData
#' @param processedInputData Output from [buildModelData]
#'
#' @return Generates the mixed effect data and also returns the `contrasts`
#'  info from [getLqmmContrasts].
buildMixedEffectMatrix <- function(modelCall, processedInputData) {
  stopifnot(inherits(processedInputData, "lqmm.processed_input_data"))

  mmr <- model.frame(modelCall$random, processedInputData$dataMix)
  mmr <- model.matrix(as.formula(modelCall$random), data = mmr)

  structure(mmr,
            # keeping the contrasts for use in predict - not yet implemented
            contrasts = getLqmmContrasts(contr = attr(mmr, "contr"),
                                         mmf_df = processedInputData$mmf_df))
}

#' Get contrast info
#'
#' Retrieves the contrast infor
#'
#' @param contr From the [stats::model.matrix] attribute
#' @param mmf_df The model frame from [buildModelData]
#'
#' @return Contrast info
getLqmmContrasts <- function(contr, mmf_df) {
  auxContr <- lapply(mmf_df, function(el) {
    if (inherits(el, "factor") && length(levels(el)) > 1) contrasts(el)
  })

  contr <- c(contr, auxContr[is.na(match(names(auxContr), names(contr)))])
  contr <- contr[!unlist(lapply(contr, is.null))]

  return(contr)
}
