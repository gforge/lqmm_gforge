##########################################################################################
# New generics

boot <- function(object, R = 50, seed = round(runif(1, 1, 10000)), startQR = FALSE) UseMethod("boot")
extractBoot <- function(object, which = "fixed") UseMethod("extractBoot")
predint <- function(object, level = 0, alpha = 0.05, R = 50, seed = round(runif(1, 1, 10000))) UseMethod("predint")
