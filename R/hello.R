# Goodness of Fit
#
# This library presents 15 criteria of goodness of fit
# The criteria can be calculated collectively or seperately.
#
#' @title Goodness of Fit
#'
#' @description Package calculates 15 different goodness of fit statistics
#'
#' # Imports ----
#' @importFrom graphics plot
#' @importFrom stats cor sd var
#'
#' @param Observations, Predicts, nTermInAppr=2, ndigit = 3, RMSE=TRUE, RRMSE = TRUE, SDR = TRUE, CoV = TRUE, PC = TRUE, PI=TRUE, ME=TRUE, RAE=TRUE, MRAE=TRUE, MAPE=TRUE, MAD=TRUE, RSq=TRUE, ARSq=TRUE, AIC=TRUE, CAIC=TRUE
#'
#' @return NULL
#'
#' @examples GoF( c(1, 2, 3, 4, 5), c(0.9, 1.97, 3.01, 4.03, 5.1))
#'
#' @export
# Goodness of Fit ----
GoF <- function(Observations, Predicts,
                nTermInAppr=2,
                ndigit = 3,
                RMSE= TRUE,
                RRMSE = TRUE,
                SDR = TRUE,
                CV = TRUE,
                PC = TRUE,
                PI= TRUE,
                ME= TRUE,
                RAE= TRUE,
                MRAE= TRUE,
                MAPE= TRUE,
                MAD= TRUE,
                RSq= TRUE,
                ARSq= TRUE,
                AIC= TRUE,
                CAIC= TRUE
                )
{
  k <- 0
  criterion <- character()
  value <- numeric()

  if(RMSE) {
    k<-k+1
    criterion[k] <- "Root mean square error (RMSE)"
    value[k] <- gofRMSE(Observations, Predicts, ndigit)
  }

  if(RRMSE) {
    k<-k+1
    criterion[k] <- "Relative root mean square error (RRMSE)"
    value[k] <- gofRRMSE(Observations, Predicts, ndigit)
  }

  if(SDR) {
    k<-k+1
    criterion[k] <- "Standard deviation ratio (SDR)"
    value[k] <- gofSDR(Observations, Predicts, ndigit)
  }

  if(CV) {
    k<-k+1
    criterion[k] <- "Coefficient of variation (CV)"
    value[k] <- gofCV(Observations, Predicts, ndigit)
  }

  if(PC) {
    k<-k+1
    criterion[k] <- "Pearson's correlation coefficients (PC)"
    value[k] <- gofPC(Observations, Predicts, ndigit)
  }

  if(PI) {
    k<-k+1
    criterion[k] <- "Performance index (PI)"
    value[k] <- gofPI(Observations, Predicts, ndigit)
  }

  if(ME) {
    k<-k+1
    criterion[k] <- "Mean error (ME)"
    value[k] <- gofME(Observations, Predicts, ndigit)
  }

  if(RAE) {
    k<-k+1
    criterion[k] <- "Relative approximation error (RAE)"
    value[k] <- gofRAE(Observations, Predicts, ndigit)
  }

  if(MRAE) {
    k<-k+1
    criterion[k] <- "Mean relative approximation error (MRAE)"
    value[k] <- gofMRAE(Observations, Predicts, ndigit)
  }

  if(MAPE) {
    k<-k+1
    criterion[k] <- "Mean absolute percentage error (MAPE)"
    value[k] <- gofMAPE(Observations, Predicts, ndigit)
  }

  if(MAD) {
    k<-k+1
    criterion[k] <- "Mean absolute deviation (MAD)"
    value[k] <- gofMAD(Observations, Predicts, ndigit)
  }

  if(RSq) {
    k<-k+1
    criterion[k] <- "Coefficient of determination (Rsq)"
    value[k] <- gofRSq(Observations, Predicts, ndigit)
  }

  if(ARSq) {
    k<-k+1
    criterion[k] <- "Adjusted coefficient of determination (ARsq)"
    value[k] <- gofARSq(Observations, Predicts, nTermInAppr, ndigit)
  }

  if(AIC) {
    k<-k+1
    criterion[k] <- "Akaike's information cCriterion (AIC)"
    value[k] <-gofAIC(Observations, Predicts, nTermInAppr, ndigit)
  }

  if(CAIC) {
    k<-k+1
    criterion[k] <- "Corrected Akaike's information criterion (CAIC)"
    value[k] <-caic <- gofCAIC(Observations, Predicts, nTermInAppr, ndigit)
  }

  plot(Observations, Predicts)

  return (data.frame(criterion, value))
}


# 1. SDR: Standard Deviation Ratio ----
#' @export
gofSDR <- function (Obs, Prd, dgt=3){
  StandardDeviationRatio <- round(sd((Obs-Prd))/sd(Obs), digits=dgt)
  return(StandardDeviationRatio)
}

# 2. CV: Coefficient of Variation ----
#' @export
gofCV <- function (Obs, Prd, dgt=3){
  CoeficientOfVariation <- round(sd(Obs - Prd)*100/mean(Obs), digits=2)
  return(CoeficientOfVariation)
}

# 3. RRMSE: Relative Root Mean Square Error ----
#' @export
gofRRMSE <- function (Obs, Prd, dgt=3){
  RelativeRootMeanSquareError=round(sqrt(mean((Obs-Prd)^2))*100/mean(Obs), digits=dgt)
  return(RelativeRootMeanSquareError)
}

# 4. PC: Pearson's Correlation Coefficients ----
#' @export
gofPC <- function (Obs, Prd, dgt=3){
  PearsonCorrelation=round(cor(Obs, Prd), digits = 3)
  return(PearsonCorrelation)
}

# 5. RMSE: Root Mean Square Error ----
#' @export
gofRMSE <- function (Obs, Prd, dgt=3){
  RootMeanSquareError=round(mean((Obs - Prd)^2), digits = 3)
  return(RootMeanSquareError)
}

# 6. PI: Performance Index ----
#' @export
gofPI <- function (Obs, Prd, dgt=3){
  RRMSE = round(sqrt(mean((Obs-Prd)^2))*100/mean(Obs), digits=dgt)
  PerformanceIndex=round(RRMSE/(1+cor(Obs, Prd)), digits = 3)
  return(PerformanceIndex)
}

# 7. ME: Mean Error ----
#' @export
gofME <- function (Obs, Prd, dgt=3){
  MeanError=round(mean(Obs - Prd), digits = 3)
  return(MeanError)
}


# 8. RAE: Global Relative Approximation Error ----
#' @export
gofRAE <- function (Obs, Prd, dgt=3){
  RelativeApproximationError=round(sum((Obs - Prd)^2)/sum(Obs^2), digits = dgt)
  return(RelativeApproximationError)
}

# 9. MRAE: Mean Relative Approximation Error ----
#' @export
gofMRAE <- function (Obs, Prd, dgt=3){
  ErrorSq = (Obs - Prd) ^2
  n = length(Obs)
  MeanRelativeApproximationError=round(sqrt(sum(ErrorSq))/sqrt(n*sum(Obs^2)), digits=dgt)
  return(MeanRelativeApproximationError)
}

# 10. MAPE: Mean Absolute Percentage Error ----
#' @export
gofMAPE <- function (Obs, Prd, dgt=3){
  MeanAbsolutePercentageError=round(mean(abs((Obs-Prd)/Obs))*100, digits=dgt)
  return(MeanAbsolutePercentageError)
}

# 11. MAD: Mean Absolute Deviation ----
#' @export
gofMAD <- function (Obs, Prd, dgt=3){
  MeanAbsoluteDeviation=round(sum(abs(Obs-Prd)/(length(Obs))), digits=dgt)
  return(MeanAbsoluteDeviation)
}

# 12. CoD / RSq: Coefficient of Determination (R-Squared) ----
#' @export
gofCoD <- function (Obs, Prd, dgt=3){
  CoefficientofDetermination=round(1-(sum((Obs-Prd)^2)/(var(Obs)*(length(Obs)-1))),
                                   digits=dgt)
  return(CoefficientofDetermination)
}

#' @export
gofRSq <- function (Obs, Prd, dgt=3){
  return(gofCoD(Obs, Prd, dgt))
}

# 13. ACoD / ARSq: Adjusted Coefficient of Determination (Adjusted R-Squared) ----
#' @export
gofACoD <- function (Obs, Prd, nTermInAppr=2, dgt=3){
    RSq = round(1-(sum((Obs-Prd)^2)/(var(Obs)*(length(Obs)-1))), digits=dgt)
    n=length(Obs)
    AdjustedCoefficientofDetermination=round(1-((1- (RSq))*(n-1)/(n-nTermInAppr-1)),
                                             digits=dgt)
    return(AdjustedCoefficientofDetermination)
}

#' @export
gofARSq <- function (Obs, Prd, nTermInAppr=2, dgt=3){
  ARsquared <- gofACoD(Obs, Prd, nTermInAppr, dgt)
  return(ARsquared)
}

# 14. AIC: Akaike's Information Criterion ----
#' @export
gofAIC <- function (Obs, Prd, nTermInAppr=2, dgt=3){
    n=length(Obs)
    AkaikesInformationCriterion=round(n*log(mean((Obs-Prd)^2),
                                      base=exp(1))+2*nTermInAppr, digits=dgt)
    if(n/nTermInAppr <= 40){
      warning("<> Use gofCAIC function (Corrected Akaike's Information Criterion for n/k is not greater than 40!")
    }
  return(AkaikesInformationCriterion)
}


# 15. CAIC: Corrected Akaike's Information Criterion ----
#' @export
gofCAIC <- function (Obs, Prd, nTermInAppr=2, dgt=3){
    n=length(Obs)
    k=nTermInAppr
    CorrectedAkaikesInformationCriterion=round(n*log(mean((Obs-Prd)^2),
                                                     base=exp(1))+(2*k)+(2*k*(k+1)/(n-k-1)), digits=dgt)
    if(n/k > 40) warning ("<> Use gofAIC (Akaike's Information Criterion) for n/k is greater than 40!")
    return(CorrectedAkaikesInformationCriterion)
}

# End
