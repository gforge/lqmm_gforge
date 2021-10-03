buildModelData <- function(modelCall, data, extraModelFrameArgs = list()) {
  stopifnot(is.list(extraModelFrameArgs))

  # Handles nested environments
  fixed <- retrieveExpression(modelCall$fixed)
  random <- retrieveExpression(modelCall$random)

  groupFormula <- asOneSidedFormula(modelCall$group)

  # Build basic input data.frame
  mfArgs <- c(list(formula = asOneFormula(random, fixed, groupFormula[[2]]),
                   data = data,
                   na.action = modelCall$na.action),
              extraModelFrameArgs)
  dataMix <- do.call("model.frame", mfArgs, envir = parent.frame())

  origOrder <- row.names(dataMix)
  for (i in names(modelCall$contrasts)) modelCall$contrasts(dataMix[[i]]) = modelCall$contrasts[[i]]

  # Reorder into groups as assigned by the group formula
  group <- model.frame(groupFormula, dataMix)
  ord <- order(unlist(group, use.names = FALSE))
  group <- group[ord, , drop = TRUE]
  ngroups <- length(unique(group))

  dataMix <- dataMix[ord, , drop = FALSE]
  revOrder <- match(origOrder, row.names(dataMix))

  # Prep mmf
  mmf_df = model.frame(fixed, dataMix)
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

buildMixedEffectMatrix <- function(modelCall, processedInputData) {
  stopifnot(inherits(processedInputData, "lqmm.processed_input_data"))

  mmr <- model.frame(modelCall$random, processedInputData$dataMix)
  mmr <- model.matrix(as.formula(modelCall$random), data = mmr)

  structure(mmr,
            # keeping the contrasts for use in predict - not yet implemented
            contrasts = getLqmmContrasts(contr = attr(mmr, "contr"),
                                         mmf_df = processedInputData$mmf_df))
}

getLqmmContrasts <- function(contr, mmf_df) {
  auxContr <- lapply(mmf_df, function(el) {
    if (inherits(el, "factor") && length(levels(el)) > 1) contrasts(el)
  })

  contr <- c(contr, auxContr[is.na(match(names(auxContr), names(contr)))])
  contr <- contr[!unlist(lapply(contr, is.null))]

  return(contr)
}
