\name{Goodness of Fit - Pearson's Correlation Coefficients}

\alias{gofPC}

\title{
Pearson's Correlation Coefficients
}

\description{
Calculates and returns Pearson's correlation coefficients (PC).
}

\usage{
gofPC(Obs, Prd, dgt = 3)
}

\arguments{
  \item{Obs}{
Observed or measured values or target vector.
}
  \item{Prd}{
  Predicted or fitted values by the model. Values produced by approximation or regression.
}
  \item{dgt}{
  Number of digits in decimal places. Default is 3.
}
}

\value{
\item{PearsonCorrelation}{Pearson's correlation coefficients (PC)}
}

\references{
OBILOR Esezi Isaac, AMADI Eric Chikweru,
Test for Significance of Pearson’s Correlation Coefficient,
International Journal of Innovative Mathematics, Statistics & Energy Policies 6(1):11-23, Jan-Mar, 2018.


Reza Soleimani, Amir Hossein Saeedi Dehaghani, Alireza Bahadori,
A new decision tree based algorithm for prediction of hydrogen sulfide solubility in various ionic liquids,
Journal of Molecular Liquids, Volume 242, 2017, Pages 701-713, ISSN 0167-7322, https://doi.org/10.1016/j.molliq.2017.07.075.
(http://www.sciencedirect.com/science/article/pii/S0167732217305123)

%% ~put references to the literature/web site here ~
}
\author{
Prof. Dr. Ecevit EYDURAN,
TA. Alper GULBE
%%  ~~who you are~~
}
\examples{
# dummy inputs, independent variable
# integers from 0 to 19
inputs <- 0:19

# dummy targets/observed values, dependent variable
# a product of 2*times inputs minus 5 with some normal noise
targets <- -5 + inputs*1.2 + rnorm(20)

# linear regression model
model<-lm(targets~inputs)

# model's predicted values against targets
predicted<-model$fitted.values

# using library ehaGoF for goodness of fit.
library(ehaGoF)

# Goodness of fit - Pearson's correlation coefficient
gofPC(targets, predicted)
}

