\name{Goodness of Fit}

\alias{GoF}

\title{Goodness of Fit}

\usage{
GoF(Observations, Predicts,
                nTermInAppr = 2,
                ndigit = 3,
                RMSE = TRUE,
                RRMSE = TRUE,
                SDR = TRUE,
                CV = TRUE,
                PC = TRUE,
                PI = TRUE,
                ME = TRUE,
                RAE = TRUE,
                MRAE = TRUE,
                MAPE = TRUE,
                MAD = TRUE,
                RSq = TRUE,
                ARSq = TRUE,
                AIC = TRUE,
                CAIC = TRUE)
}

\description{
Tests predicted and observed values for the goodness of fit with various criteria. The goodness of fit tests are used to test how well the model fits. Measures of goodness of fit typically summarize the argument between targets or observed values and the values expected or predicted under the model in question.
}

\references{
Comparison of the Predictive Capabilities of Several Data Mining Algorithms and Multiple Linear Regression in the Prediction of Body Weight by Means of Body Measurements in the Indigenous Beetal Goat of Pakistan, Ecevit Eyduran, Daniel Zaborski, Abdul Waheed, Senol Celik, Koksal Karadas and Wilhelm Grzesiak, Pakistan J. Zool., vol. 49(1), pp 257-265, 2017.
%% ~put references to the literature/web site here ~
}

\author{
Prof. Dr. Ecevit Eyduran,
TA. Alper Gulbe

}
\arguments{
  \item{Observations}{
  Observed values or target vector.
}
  \item{Predicts}{
Predicted values. Values produced or fitted by approximation or regression.
}
  \item{nTermInAppr}{
Number of terms used in approximation or regression model. Generally 2 for simple linear model. Default is 2.
}  \item{ndigit}{
Number of digits in decimal places. Default is 3.
}  \item{RMSE}{
Whether to show Root Mean Square Error statistics. Default is TRUE.
}  \item{RRMSE}{
Whether to show Relative Root Mean Square Error statistics. Default is TRUE.
}  \item{SDR}{
Whether to show Standard Deviation Ratio statistics. Default is TRUE.
}  \item{CV}{
Whether to show Coefficient of Variance statistics. Default is TRUE.
}  \item{PC}{
Whether to show Pearson's Correlation Coefficients statistics. Default is TRUE.
}  \item{PI}{
Whether to show Performance Index statistics. Default is TRUE.
}  \item{ME}{
Whether to show Mean Error statistics. Default is TRUE.
}  \item{RAE}{
Whether to show Global Relative Approximation Error statistics. Default is TRUE.
}  \item{MRAE}{
Whether to show Modified Relative Approximation Error statistics. Default is TRUE.
}  \item{MAPE}{
Whether to show Mean Absolute Percentage Error statistics. Default is TRUE.
}  \item{MAD}{
Whether to show Mean Absolute Deviation statistics. Default is TRUE.
}  \item{RSq}{
Whether to show Coefficient of Determination (R-Squared) statistics. Default is TRUE.
}  \item{ARSq}{
Whether to show Adjusted Coefficient of Determination (Adjusted R-Squared) statistics. Default is TRUE.
Warning: nTermInAppr must be supplied.
}  \item{AIC}{
Whether to show Akaike's Information Criterion statistics. Default is TRUE.
Warning: nTermInAppr must be supplied.
}  \item{CAIC}{
Whether to show Corrected Akaike's Information Criterion statistics. Default is TRUE.
Warning: nTermInAppr must be supplied.
}
}

\examples{
# dummy inputs, independent variable
# integers from 0 to 9
inputs <- 0:9

# dummy targets/observed values, dependent variable
# a product of 2*times inputs minus 5 with some normal noise
targets <- -5 + inputs*2 + rnorm(10)

# linear regression model
model<-lm(targets~inputs)

# About the model
summary(model)

# Number of Terms
n = length(model$coefficients)

# model's predicted values against targets
predicted<-model$fitted.values

# using library ehaGoF for goodness of fit.
library(ehaGoF)

# Goodness of fit statistics
GoF(targets, predicted, nTermInAppr=n)
}
