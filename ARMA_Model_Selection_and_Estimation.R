################################################################################
#                                                                              #
#           W4451 - Applied project using R, summer semester 2022              #
#                                                                              #
#              Problem 1: ARMA Model Selection and Estimation                  #
#                                                                              #
################################################################################

# Information:

# Use this file to save your code for solving Problem 1 of the applied project.
# Remove '1234567' in line 17 and instead, replace it with the matriculation 
# number of one of your group members.
#1,a)
# Start at line 35 with your own code.

MatrNr <- 6857763

# Set up your current R session for working on Problem 1
# NOTE: online connection required!

if (!require(devtools)) {
  install.packages("devtools", repos = "http://cran.us.r-project.org")
}

devtools::source_gist(
  "https://gist.github.com/dschulz13/3d4fdd7a5cc1067520e38c406207fe37",
  verbose = FALSE, quiet = TRUE,
  local = globalenv(),
  filename = "P1-Setup.R"
)

# <--------------- Begin in the next line with your own code  ---------------> #

#Aufgabe 1 b)-------------------------------------------------------------------

# erforderliche Pakete laden
library(zoo)
library(ggplot2)
library(forecast)

# xt Zeitreihe Diagramm plotten
plot_xt <- autoplot.zoo(xt) +
  xlab("Zeit") +
  ylab("xt") +
  ggtitle("Die in xt gespeicherten simulierten Zeitreihen")
plot_xt

#Aufgabe 1 c)-------------------------------------------------------------------

# autocorrelation von lag = 0,...,8 berechnen
acf(xt, lag.max = 8, plot=FALSE)$acf

# korrelogramm von xt Reihe mit lag = 0,...,8 zeichnen
acf_xt <- ggAcf(coredata(xt), lag.max = 8) + 
  ggtitle("Korrelogramm von xt Reihe mit lag = 0,...,8")
acf_xt

# Anzahl von Beobachtungen
length(xt)

# mittelwert von xt
mean(xt)

# varianz von xt
var(xt)

#Aufgabe 1 d)-------------------------------------------------------------------

# alle 36 ARMA(p,q) Modelle für p = 0,1,...,5 und q = 0,1,...,5 erläutern
n=751
p.max=5
q.max=5
AIC = BIC = matrix(0, nrow = p.max + 1, ncol=q.max + 1)
for(p in 0:5) {
  for(q in 0:5) {
    model_ARMA=arima(xt,order=c(p,0,q))
    AIC[p+1, q+1]=model_ARMA$aic
    BIC[p+1, q+1]=-2*model_ARMA$loglik+log(n)*(p+q)
  }
}

#beste Modell laut AIC:
ARMA_AIC = arima(xt,order=c(2,0,2))
ARMA_AIC

#minimum AIC
min(ARMA_AIC) #-> 6418.907
ARMA_AIC #zeile 2 -> p=2 ; spalte 2 -> q=2 => ARMA(2,2)-Modell

#beste Modell laut BIC:
ARMA_BIC = arima(xt,order=c(2,0,2))
ARMA_BIC

#minimum BIC
min(ARMA_BIC) #-> 6433.392
ARMA_BIC #zeile 2 -> p=2 ; spalte 2 -> q=2 => ARMA(2,2)-Modell

#Aufgabe 1 e)-------------------------------------------------------------------

# ausgewählte ARMA Modell von d)
model_ARMA=arima(xt,order=c(2,0,2))
model_ARMA

#Aufgabe 1 f)-------------------------------------------------------------------

#summe von autocovarianz berechnen
ar_coef = sum(ARMA_BIC$coef[1:2])
ma_coef = sum(ARMA_BIC$coef[2:4])
sig2_eps = ARMA_BIC$sigma2
test2 = ((1 + ar_coef) / (1 - ma_coef))^2 * sig2_eps
test2

#varianz von unserer sample mean var(x??)
n=751
ar_coef = sum(ARMA_BIC$coef[1:2])
ma_coef = sum(ARMA_BIC$coef[2:4])
sig2_eps = ARMA_BIC$sigma2
test3 = ((1 + ar_coef) / (1 - ma_coef))^2 * sig2_eps / n
test3
