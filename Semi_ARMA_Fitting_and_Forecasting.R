################################################################################
#                                                                              #
#           W4451 - Applied project using R, summer semester 2022              #
#                                                                              #
#               Problem 2: Semi-ARMA Fitting and Forecasting                   #
#                                                                              #
################################################################################

# Information:

# Use this file to save your code for solving Problem 2 of the applied project.

# Set up your current R session for working on Problem 2 
# NOTE: online connection required!

if (!require(devtools)) {
  install.packages("devtools", repos = "http://cran.us.r-project.org")
}

devtools::source_gist(
  "https://gist.github.com/dschulz13/3d4fdd7a5cc1067520e38c406207fe37", 
  local = globalenv(), verbose = FALSE, quiet = TRUE,
  filename = "P2_3_4_Setup.R")

# <--------------- Begin in the next line with your own code  ---------------> #

library('zoo')
library('ggplot2')
library('smoots')
library('dplyr')
library('forecast')

#Aufgabe 2 a)-------------------------------------------------------------------

data <- read.csv("macro_monthly.csv")                                           #Éinlesen der Daten 

head(data)

pce <- data$pce[-(469:490)]                                                     #Adjustierung bezüglich der Zeitkomponente (nur bis zum 01.12.2019)

pce <- ts(pce, start = c(1981, 1), frequency = 12)        

#Plot des PCE auf die Zeit 
plot_pce <- autoplot.zoo(pce) + xlab("Jahr") + ylab("PCE") + ggtitle("PCE von Anfang 1981 bis Ende 2019")
plot_pce


#Aufgabe 2 b)-------------------------------------------------------------------

#Plot des Log PCE auf die Zeit
log_pce <- log(pce)
plot_log_pce <- autoplot.zoo(log_pce) + xlab("Jahr") + ylab("Log PCE") + ggtitle("Log PCE von Anfang 1981 bis Ende 2019")
plot_log_pce


#Aufgabe 2 c)-------------------------------------------------------------------

#Aufteilen des Datensatzes in training- und test-Daten
pce_training <- log(pce[-(464:468)])
pce_training <- ts(pce_training, start = c(1981, 1), frequency = 12)
plot_log_pce_training <- autoplot.zoo(pce_training) + xlab("Jahr") + ylab("Log PCE") + ggtitle("Log PCE von Anfang 1981 bis Mitte 2019 - Training Set")
plot_log_pce_training

length(pce_training)

pce_testing <- log(pce [464:468])
pce_testing <- ts(pce_testing, start = c(2019, 8), frequency = 12)              
#plot_pce_testing <- autoplot.zoo(as.zoo(pce_testing)) + xlab('a') + ylab('b')
#plot_pce_testing 

#Smoots Schätzung
est <- msmooth(pce_training)
bwidth <- est$b0
bwidth

trend <- fitted(est)
df_pce_training <- data.frame(t = time(pce_training), trend = trend)

plot_pce_training <- plot_log_pce_training +                                             
  geom_line(data = df_pce_training, aes(x = t, y = trend), color = "red", size = 0.7) +
  ggtitle(paste0("Log-transformierte PCE Zeitreihe (schwarz) mit geschätzem lokal linearen Trend (rot)"))
plot_pce_training

#Risiduenbestimmung und -abbildung
res <- resid(est)

plot_res_pce_training <- autoplot.zoo(res) +
  xlab("Jahr") +
  ylab("Residuum Wert") +
  ggtitle("Residuumzeitreihe") +
  geom_hline(yintercept = 0, color = "blue")
plot_res_pce_training


#Aufgabe 2 d)-------------------------------------------------------------------

test_linear <- confBounds(est, plot = TRUE, x = c(time(pce_training)), xlab = "Jahr", 
                          ylab = "Log PCE")
test_linear
test_linear$b.ub


#Aufgabe 2 e)-------------------------------------------------------------------

acf_res_pce_training <- ggAcf(as.numeric(res)) +
  ggtitle("Korrelogramm der trendbereinigten Zeitreihe")
acf_res_pce_training


#Aufgabe 2 f)-------------------------------------------------------------------

model_selection <- modelCast(est)                                               
?modelCast


#Aufgabe 2 g)-------------------------------------------------------------------

p_max <- q_max <- 5
p <- 0:p_max
q <- 0:q_max

mase <- matrix(NA, ncol = q_max + 1, nrow = p_max + 1)
rownames(mase) <- paste0("p=", p)
colnames(mase) <- paste0("q=", q)
rmsse <- mase

for (p0 in p) {
  for (q0 in q) {
    backt <- rollCast(log_pce, p = p0, q = q0, plot = FALSE)
    mase[p0 + 1, q0 + 1] <- backt$MASE
    rmsse[p0 + 1, q0 + 1] <- backt$RMSSE
  }
}

pq_opt_mase <- unname(which(mase == min(mase), arr.ind = TRUE) - 1)
p_opt_mase <- pq_opt_mase[1]
q_opt_mase <- pq_opt_mase[2]

pq_opt_rmsse <- unname(which(rmsse == min(rmsse), arr.ind = TRUE) - 1)
p_opt_rmsse <- pq_opt_rmsse[1]
q_opt_rmsse <- pq_opt_rmsse[2]

pq_opt_mase
pq_opt_rmsse

mase
rmsse


backtest <- rollCast(log_pce, p = 1, q = 1, plot = FALSE)                   
accuracy(backtest$fcast.roll[1, ], x = backtest$y.out)


#Aufgabe 2h)--------------------------------------------------------------------

p_opt <- 1
q_opt <- 3

est2 <- msmooth(log_pce)                                                        
est2$b0

res2 <- resid(est2)
arma_opt <- arima(res2, order = c(p_opt, 0, q_opt), include.mean = FALSE)
arma_opt


et <- arma_opt$residuals

df <- data.frame(et = c(et))
plot_qq <- ggplot(df, aes(sample = et)) +
  geom_qq() +
  geom_qq_line(color = "red") +
  xlab("Theoretische Quantile") +
  ylab("Stichproben Quantile") +
  ggtitle("Normal Q-Q Plot")
plot_qq


#Aufgabe 2i)--------------------------------------------------------------------

set.seed(123)
fc_pce <- modelCast(est2, p = p_opt, q = q_opt, h = 5, method = "boot")
fc_pce

#Retransformation der log Daten
fc_normal_pce <- exp(fc_pce)
fc_normal_pce

#Plot der zukünftigen Schätzungen und deren Intervalle 
fcast <- ts(fc_normal_pce[1, ], start = c(2020, 1), frequency = 12)
df <- data.frame(
  t = c(time(fcast)),
  fcast = fc_normal_pce[1, ],
  LB = fc_normal_pce[2, ],
  UB = fc_normal_pce[3, ]
)

plot_est<- autoplot.zoo(as.zoo(window(pce, start = c(2018, 1)))) +            
  geom_ribbon(data = df, aes(x = t, ymin = LB, ymax = UB), 
              fill = alpha("red", 0.2), color = NA) +
  geom_line(data = df, aes(x = t, y = fcast), col = "red") +
  xlab("Jahr") +
  ylab("PCE") +
  ggtitle("Ende der PCE-Zeitreihe von Anfang 2018 bis Ende 2019 mit Prognose für 5 weitere Monate")
plot_est