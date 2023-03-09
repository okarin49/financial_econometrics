################################################################################
#                                                                              #
#           W4451 - Applied project using R, summer semester 2022              #
#                                                                              #
#                 Problem 3: GARCH Models and Backtesting                      #
#                                                                              #
################################################################################

# Information:

# Use this file to save your code for solving Problem 3 of the applied project.

# Set up your current R session for working on Problem 3
# NOTE: online connection required!

if (!require(devtools)) {
  install.packages("devtools", repos = "http://cran.us.r-project.org")
}

devtools::source_gist(
  "https://gist.github.com/dschulz13/3d4fdd7a5cc1067520e38c406207fe37", 
  local = globalenv(), verbose = FALSE, quiet = TRUE,
  filename = "P2_3_4_Setup.R")

# <--------------- Begin in the next line with your own code  ---------------> #

# erforderliche Pakete laden
library(zoo)
library(ggplot2)
library(ggpubr)
library(forecast)
library(fGarch)
library(ufRisk)

### a.

# daten einlesen
AMZ.DE <- read.csv("AMZ.DE.csv", header = TRUE)

# fünf erste und letzte Beobachtungen anzeigen
head(AMZ.DE)
tail(AMZ.DE)

# plot der Schlusskurse
Date <- as.Date(AMZ.DE$Date)
AMZ.DE_Close <- zoo(AMZ.DE$Close, order.by = Date)

plot_AMZ.DE_Close <- autoplot.zoo(AMZ.DE_Close) +
  xlab("Jahr") +
  ylab("Schlusskurse in Euro") +
  ggtitle("Tagesschlusskurse von Amazon in Deutschland")
plot_AMZ.DE_Close

### b.

# Log-Renditen aus den Schlusskursen berechnen
AMZ.DE_returns <- diff(log(AMZ.DE_Close))

# plot der Renditen
plot_AMZ.DE_returns <- autoplot.zoo(AMZ.DE_returns) +
  xlab("Jahr") +
  ylab("Log-Renditen") +
  geom_hline(yintercept = mean(AMZ.DE_returns), color = "red") +
  ggtitle("Log-Renditen der Schlusskurse von Amazon Deutschland")
plot_AMZ.DE_returns

# plot von ACF der Log-Renditen
plot_acf_1 <- ggAcf(as.numeric(AMZ.DE_returns)) +
  ggtitle("Korrelogramm der Log-Renditen")
plot_acf_1

# plot von ACF der quadrierten Log-Renditen
AMZ.DE_squared_returns <- AMZ.DE_returns^2
plot_acf_2<- ggAcf(as.numeric(AMZ.DE_squared_returns)) +
  ggtitle("Korrelogramm der quadrierten Log-Renditen")
plot_acf_2

# beide Korrelogramme in derselben Grafik anzeigen
plot_acf <- ggarrange(plot_acf_1, plot_acf_2, ncol = 1, nrow = 2)
plot_acf

### c.

# renditen daten in training und test daten aufteilen
n_test <- 250
n_AMZ.DE_returns <- length(AMZ.DE_returns)
n_train <- n_AMZ.DE_returns - n_test

# training data
AMZ.DE_returns_training <- head(AMZ.DE_returns, n_train)
AMZ.DE_returns_training

# test data
AMZ.DE_returns_test <- tail(AMZ.DE_returns, n_test)
AMZ.DE_returns_test

# fit GARCH(1,1)
garch11 <- garchFit(~ garch(1, 1), AMZ.DE_returns_training, trace = FALSE)
garch11

# fit APARCH(1,1)
aparch11 <- garchFit(~ aparch(1, 1), AMZ.DE_returns_training, trace = FALSE)
aparch11

# fit GARCH-t(1,1)
garch11_t <- garchFit(~ garch(1, 1), AMZ.DE_returns_training, trace = FALSE, 
                      cond.dist = "std")
garch11_t

# fit APARCH-t(1,1)
aparch11_t <- garchFit(~ aparch(1, 1), AMZ.DE_returns_training, trace = FALSE, 
                       cond.dist = "std")
aparch11_t

# GARCH(1,1) BIC
garch11@fit$ics[["BIC"]]

# APARCH(1,1) BIC
aparch11@fit$ics[["BIC"]]

# GARCH-t(1,1) BIC
garch11_t@fit$ics[["BIC"]]

# APARCH-t(1,1) BIC
aparch11_t@fit$ics[["BIC"]]

# beste Modell laut BIC angeben
garch_opt <- garchFit(~ garch(1, 1), AMZ.DE_returns_training, trace = FALSE, 
                      cond.dist = "std")
garch_opt

### d.

# one-step rolling forecasts GARCH(1,1)
forecast_garch11 <- varcast(as.numeric(AMZ.DE_Close), model = "sGARCH", 
                            distr = "norm", garchOrder = c(1, 1))

# one-step rolling forecasts APARCH(1,1)
forecast_aparch11 <- varcast(as.numeric(AMZ.DE_Close), model = "apARCH",
                             distr = "norm", garchOrder = c(1, 1))

# one-step rolling forecasts GARCH-t(1,1)
forecast_garch11_t <- varcast(as.numeric(AMZ.DE_Close), model = "sGARCH",
                              distr = "std", garchOrder = c(1, 1))

# one-step rolling forecasts APARCH-t(1,1)
forecast_aparch11_t <- varcast(as.numeric(AMZ.DE_Close), model = "apARCH",
                               distr = "std", garchOrder = c(1, 1))

# plot der geschätzten 97.5% VaR und 97.5% ES Reihen gegen Testverluste (GARCH(1,1))
Date_test <- tail(Date, n_test)
VaR99_garch11 <- zoo(forecast_garch11$VaR.v, order.by = Date_test)
VaR975_garch11 <- zoo(forecast_garch11$VaR.e, order.by = Date_test)
ES975_garch11 <- zoo(forecast_garch11$ES, order.by = Date_test)
Loss_test <- tail(-AMZ.DE_returns, 250)

violations_garch11 <- Loss_test > VaR975_garch11
Date_violations_garch11 <- Date_test[violations_garch11]
VaR_violations_garch11 <- VaR975_garch11[violations_garch11]

n_violations_garch11 <- sum(violations_garch11) # erwartete Anzahl von "violations"
n_violations_garch11

data_garch11 <- data.frame(Date = Date_test, Var975 = VaR975_garch11, 
                           ES975 = ES975_garch11, Loss = Loss_test)
data_garch11_vio <- data.frame(Date = Date_violations_garch11, 
                               VaR = VaR_violations_garch11)

colors <- c("Loss_test" = "grey64", "VaR975" = "red", "ES975" = "green")

plot_loss_garch11 <- ggplot(data_garch11, aes(x = Date)) +
  geom_segment(aes(y = Loss, xend = Date, yend = 0, color = "Loss_test")) +
  geom_line(aes(y = VaR975_garch11, color = "VaR975")) +
  geom_line(aes(y = ES975_garch11, color = "ES975")) +
  geom_point(data = data_garch11_vio, aes(x = Date, y = VaR), color = "blue", size = 3,
             pch = 13) +
  scale_color_manual(values = colors,
                     labels = c("Verluste", "97.5%-VaR", "97.5%-ES"),
                     name = "Reihen") +
  xlab("Monat und Jahr") +
  ylab("Verlust, 97.5%-VaR und 97.5%-ES") +
  ggtitle("Testverluste zusammen mit den 97.5% Risikomaßen - GARCH (1,1)")
plot_loss_garch11

# plot der geschätzten 97.5% VaR und 97.5% ES Reihen gegen Testverluste (APARCH(1,1))
VaR99_aparch11 <- zoo(forecast_aparch11$VaR.v, order.by = Date_test)
VaR975_aparch11 <- zoo(forecast_aparch11$VaR.e, order.by = Date_test)
ES975_aparch11 <- zoo(forecast_aparch11$ES, order.by = Date_test)
Loss_test <- tail(-AMZ.DE_returns, 250)

violations_aparch11 <- Loss_test > VaR975_aparch11
Date_violations_aparch11 <- Date_test[violations_aparch11]
VaR_violations_aparch11 <- VaR975_aparch11[violations_aparch11]

n_violations_aparch11 <- sum(violations_aparch11) # erwartete Anzahl von "violations"
n_violations_aparch11

data_aparch11 <- data.frame(Date = Date_test, Var975 = VaR975_aparch11, 
                           ES975 = ES975_aparch11, Loss = Loss_test)
data_aparch11_vio <- data.frame(Date = Date_violations_aparch11, 
                                VaR = VaR_violations_aparch11)

colors <- c("Loss_test" = "grey64", "VaR975_aparch11" = "red", 
            "ES975_aparch11" = "green")

plot_loss_aparch11 <- ggplot(data_aparch11, aes(x = Date)) +
  geom_segment(aes(y = Loss, xend = Date, yend = 0, color = "Loss_test")) +
  geom_line(aes(y = VaR975_aparch11, color = "VaR975_aparch11")) +
  geom_line(aes(y = ES975_aparch11, color = "ES975_aparch11")) +
  geom_point(data = data_aparch11_vio, aes(x = Date, y = VaR), color = "blue", size = 3,
             pch = 13) +
  scale_color_manual(values = colors,
                     labels = c("Verluste", "97.5%-VaR", "97.5%-ES"),
                     name = "Reihen") +
  xlab("Monat und Jahr") +
  ylab("Verlust, 97.5%-VaR und 97.5%-ES") +
  ggtitle("Testverluste zusammen mit den 97.5% Risikomaßen - APARCH (1,1)")
plot_loss_aparch11

# plot der geschätzten 97.5% VaR und 97.5% ES Reihen gegen Testverluste (GARCH-t(1,1))
VaR99_garch11_t <- zoo(forecast_garch11_t$VaR.v, order.by = Date_test)
VaR975_garch11_t <- zoo(forecast_garch11_t$VaR.e, order.by = Date_test)
ES975_garch11_t <- zoo(forecast_garch11_t$ES, order.by = Date_test)
Loss_test <- tail(-AMZ.DE_returns, 250)

violations_garch11_t <- Loss_test > VaR975_garch11_t
Date_violations_garch11_t <- Date_test[violations_garch11_t]
VaR_violations_garch11_t <- VaR975_garch11_t[violations_garch11_t]

n_violations_garch11_t <- sum(violations_garch11_t) # erwartete Anzahl von "violations"
n_violations_garch11_t

data_garch11_t <- data.frame(Date = Date_test, Var975 = VaR975_garch11_t, 
                              ES975 = ES975_garch11_t, Loss = Loss_test)
data_garch11_t_vio <- data.frame(Date = Date_violations_garch11_t, 
                                  VaR = VaR_violations_garch11_t)

colors <- c("Loss_test" = "grey64", "VaR975_garch11_t" = "red", 
            "ES975_garch11_t" = "green")

plot_loss_garch11_t <- ggplot(data_garch11_t, aes(x = Date)) +
  geom_segment(aes(y = Loss, xend = Date, yend = 0, color = "Loss_test")) +
  geom_line(aes(y = VaR975_garch11_t, color = "VaR975_garch11_t")) +
  geom_line(aes(y = ES975_garch11_t, color = "ES975_garch11_t")) +
  geom_point(data = data_garch11_t_vio, aes(x = Date, y = VaR), color = "blue", size = 3,
             pch = 13) +
  scale_color_manual(values = colors,
                     labels = c("Verluste", "97.5%-VaR", "97.5%-ES"),
                     name = "Reihen") +
  xlab("Monat und Jahr") +
  ylab("Verlust, 97.5%-VaR and 97.5%-ES") +
  ggtitle("Testverluste zusammen mit den 97.5% Risikomaßen - GARCH-t (1,1)")
plot_loss_garch11_t

# plot der geschätzten 97.5% VaR und 97.5% ES Reihen gegen Testverluste (APARCH-t(1,1))
VaR99_aparch11_t <- zoo(forecast_aparch11_t$VaR.v, order.by = Date_test)
VaR975_aparch11_t <- zoo(forecast_aparch11_t$VaR.e, order.by = Date_test)
ES975_aparch11_t <- zoo(forecast_aparch11_t$ES, order.by = Date_test)
Loss_test <- tail(-AMZ.DE_returns, 250)

violations_aparch11_t <- Loss_test > VaR975_aparch11_t
Date_violations_aparch11_t <- Date_test[violations_aparch11_t]
VaR_violations_aparch11_t <- VaR975_aparch11_t[violations_aparch11_t]

n_violations_aparch11_t <- sum(violations_aparch11_t) # erwartete Anzahl von "violations"
n_violations_aparch11_t

data_aparch11_t <- data.frame(Date = Date_test, Var975 = VaR975_aparch11_t, 
                            ES975 = ES975_aparch11_t, Loss = Loss_test)
data_aparch11_t_vio <- data.frame(Date = Date_violations_aparch11_t, 
                                VaR = VaR_violations_aparch11_t)

colors <- c("Loss_test" = "grey64", "VaR975_aparch11_t" = "red", 
            "ES975_aparch11_t" = "green")

plot_loss_aparch11_t <- ggplot(data_aparch11_t, aes(x = Date)) +
  geom_segment(aes(y = Loss, xend = Date, yend = 0, color = "Loss_test")) +
  geom_line(aes(y = VaR975_aparch11_t, color = "VaR975_aparch11_t")) +
  geom_line(aes(y = ES975_aparch11_t, color = "ES975_aparch11_t")) +
  geom_point(data = data_aparch11_t_vio, aes(x = Date, y = VaR), color = "blue", size = 3,
             pch = 13) +
  scale_color_manual(values = colors,
                     labels = c("Verluste", "97.5%-VaR", "97.5%-ES"),
                     name = "Reihen") +
  xlab("Monat und Jahr") +
  ylab("Verlust, 97.5%-VaR and 97.5%-ES") +
  ggtitle("Testverluste zusammen mit den 97.5% Risikomaßen - APARCH-t (1,1)")
plot_loss_aparch11_t

### e.

# unconditional coverage test fürs GARCH(1,1) Modell
covtest(forecast_garch11)

# unconditional coverage test fürs APARCH(1,1) Modell
covtest(forecast_aparch11)

# unconditional coverage test fürs GARCH-t(1,1) Modell
covtest(forecast_garch11_t)

# unconditional coverage test fürs APARCH-t(1,1) Modell
covtest(forecast_aparch11_t)

### f.

# traffic light test fürs GARCH(1,1) Modell
trafftest(forecast_garch11)

# traffic light test fürs APARCH(1,1) Modell
trafftest(forecast_aparch11)

# traffic light test fürs GARCH-t(1,1) Modell
trafftest(forecast_garch11_t)

# traffic light test fürs APARCH-t(1,1) Modell
trafftest(forecast_aparch11_t)
