################################################################################
#                                                                              #
#           W4451 - Applied project using R, summer semester 2022              #
#                                                                              #
#           Problem 4: Detecting and Analyzing the COVID-19 Crisis             #
#                                                                              #
################################################################################

# Information:

# Use this file to save your code for solving Problem 4 of the applied project.

# Set up your current R session for working on Problem 4
# NOTE: online connection required!

if (!require(devtools)) {
  install.packages("devtools", repos = "http://cran.us.r-project.org")
}

devtools::source_gist(
  "https://gist.github.com/dschulz13/3d4fdd7a5cc1067520e38c406207fe37", 
  local = globalenv(), verbose = FALSE, quiet = TRUE,
  filename = "P2_3_4_Setup.R")

# <--------------- Begin in the next line with your own code  ---------------> #

#Aufgabe 4 a)-------------------------------------------------------------------

library('fGarch')
library('zoo')
library('ggplot2')
library('ggpubr')

data <- read.csv('nasdaq.csv', header = TRUE)
Price <- data$Close

Date <- as.Date(data$Date, format = "%Y-%m-%d")

Price <- zoo(Price, order.by = Date)

n_cp <- length(Price)    
n_ret <- n_cp - 1  

Returns <- diff(log(Price))  
head(Returns)
Returns_v <- coredata(Returns)

plot_price <- autoplot.zoo(Price) +
  xlab("Jahr") +
  ylab("NASDAQ Composite") +
  ggtitle("Tägliche NASDAQ Composite Schlusskurse von Jan 2015 bis Dez 2021")
plot_price

plot_returns <- autoplot.zoo(Returns) +                                         
  xlab("Jahr") +
  ylab("NASDAQ Composite log-Renditen") +
  ggtitle("Tägliche NASDAQ Composite Renditen von Jan 2015 bis Dez 2021")
plot_returns

plot_comb <- ggarrange(plot_price, plot_returns , ncol = 1, nrow = 2, align = "v")
plot_comb

#Genaue graphische Analyse des Covid Shocks 
plot_crisis <- autoplot.zoo(window(Price, start = '2020-02-01', end = '2020-04-15')) +
  xlab("2020") +
  ylab("NASDAQ Composite") +
  ggtitle("Tägliche NASDAQ Composite Schlusskurse vom 2020-02-01 bis zum 2020-04-15")
plot_crisis

Close_subset <- window(Price, start = "2020-01-01", 
                       end = "2020-03-31")
min_ind <- which(Close_subset == min(Close_subset))
max_ind <- which(Close_subset == max(Close_subset))

time(Close_subset)[min_ind]
time(Close_subset)[max_ind]

plot_comb_2 <- ggarrange(plot_price, plot_returns ,plot_crisis, ncol = 1, nrow = 3, align = "v")
plot_comb_2


#Aufgabe 4 b)-------------------------------------------------------------------

bw <- 25

RetDiff <- matrix(0, 3, n_cp - 2 * bw+1)
DV <- c((1:bw) * 0, (1:bw)^0)                                                   #DV teilt die Daten in links- und rechtsseitige
uh <- (-(bw - 1):0) / (bw - 0.5)                                                #Kernel Variable
mu <- 2                                                                         #mu = 2 wird vorgeschlagen
KWh <- (1 + uh)^mu * (-uh)^(mu - 1)                                             #linke Hälfte der Gewichte nach Mueller/Wang 1994
KW <- c(KWh, rev(KWh))                                                          #ganzheitliche Geichte, continuous, wenn mu >= 2


for(i in (bw + 1):(n_cp - bw+1)) {
  Yi <- Returns_v[(i - bw):(i + bw-1)]
  Mi <- lm(Yi ~ DV, weights = KW)                                               #lm Funktion mit einem DV
  RetDiff[1, i - bw] <- Mi$coefficients[1]                                      #Intercept
  RetDiff[2, i - bw] <- Mi$coefficients[2]                                      #Der Koeffizient des DV, für die Differenz der Fälle DV=0 und DV=1
}

RetDiff[3, ] <- RetDiff[1, ] + RetDiff[2, ]                                     #Die schätzung bezüglich DV = 1

#Das Minimum der Differenzen entspricht dem Anfang der Covid-19 Krise
ns <- ((bw + 1):(n_cp - bw))[RetDiff[2, ] == min(RetDiff[2, ])] 
ns  

data[ns+1, ]                                                                    #Anfang der Krise

#Informationen zum Ende der Krise
ne <- ((bw + 1):(n_cp - bw))[RetDiff[2, ] == max(RetDiff[2, ])] 
ne    

data[ne, ]                                                                      #Ende der Krise 

data[ne + 1, ]                                                                  #Anfang der Erholungsphase


#Aufgabe 4 c)-------------------------------------------------------------------

pre_return <- Returns[1:(ns-1)]                                                 #Nutzung der ersten (ns-1) Renditen um ein GARCH-t(1,1) Model zu fitten
pre_return.GARCH <- garchFit(~ garch(1, 1), pre_return, cond.dist = "std", trace = FALSE) #GARCH-t-Model 
pre_return.GARCH

#Degrees of Freedom
df <- pre_return.GARCH@fit$par[["shape"]]                                       #Geschätzte df
df                                                                              #df = 4.874459 

#Angepasste konditionelle Volatilitäten
Vol.est <- pre_return.GARCH@sigma.t

#Geschätzte konditionelle Volatilitäten für die übrigen Beobachtungen
K <- n_cp - ns
GARCH.Pred <- predict(pre_return.GARCH, n.ahead = K)
Vol.Pred <- GARCH.Pred$standardDeviation

#Kombination der Volatilitäten zu einem Vektor 
Vol.all <- c(Vol.est, Vol.Pred)


#t-Quantilswerte
alpha <- 0.0005                                                                 #Festlegung des Signifikanzniveaus
qt.L <- qt(alpha, df)
qt.L

#Extrem negative Renditen
Date.Extreme.N <- Date[-1][Returns_v < Vol.all * qt.L]                          #Daten mit extrem negativen Renditen
Date.Extreme.N                                                                  #Daten: "2020-03-09" "2020-03-12" "2020-03-16"

round(Returns[Date.Extreme.N], digits = 6)                                      #Gerundete Renditen zu den entsprechenden Daten
#Extrem negative Renditen: -0.075666  -0.099099  -0.131492.

#Extrem positive Renditen
qt.U <- -qt.L
Date.Extreme.P <- Date[-1][Returns_v > Vol.all * qt.U]                          #Daten mit extrem positiven Renditen
Date.Extreme.P                                                                  #Daten: "2020-03-12", "2020-03-24"
round(Returns[Date.Extreme.P], digits = 6) 
#Extrem positive Renditen: 0.089347,   0.078085.


#Aufgabe 4 d)-------------------------------------------------------------------

date_min <- Date[1]                                                             #Das erste Datum
date_max <- tail(Date, 1)                                                       #Das letzte Datum

#1.Plot: Index Schlusskurse
p1 <- autoplot.zoo(Price) +
  ylab("") +
  xlab("Jahr") +
  ggtitle("Tägliche NASDAQ Composite Schlusskurse von Jan 2015 bis Dez 2021") +
  geom_vline(xintercept = Date[ns + 1], color = 2) +
  geom_vline(xintercept = Date[ne + 1], color = 7) +
  xlim(date_min, date_max) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
p1

#2. Plot: Renditezeitreihe 
p2 <- autoplot.zoo(Returns) +
  ylab("") +
  xlab("Jahr") +
  ggtitle("Tägliche NASDAQ Composite Renditen mit 99.9% Konfidenzintervall von Jan 2015 bis Dez 2021") +
  geom_vline(xintercept = Date[ns + 1], color = 2) +
  geom_vline(xintercept = Date[ne + 1], color = 7) +
  geom_line(data = data.frame(x = Date[-1], y = Vol.all * qt.L), aes(x = x, y = y), color = 4) +
  geom_line(data = data.frame(x = Date[-1], y = Vol.all * qt.U), aes(x = x, y = y), color = 4) +
  xlim(date_min, date_max) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
p2

#3. Plot: Differenzen zwischen links und rechtsseitigen Renditedurchschnitten
p3 <- ggplot(data.frame(x = Date[(bw + 1):(n_cp - bw+1)], y = RetDiff[2, ]), aes(x = x, y = y)) +
  geom_line() +
  ggtitle("Differenzen der rechts- und linksseitigen lokal gewichteter Renditedurchschnitte") +
  geom_vline(xintercept = Date[ns + 1], color = 2) +
  geom_vline(xintercept = Date[ne + 1], color = 7) +
  xlim(date_min, date_max) +
  xlab("Jahr") +
  ylab("") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
p3

#Gemeinsamer Plot
plot_crisis <- ggarrange(p1, p2, p3, ncol = 1, nrow = 3, align = "v")
plot_crisis


#Aufgabe 4 e)-------------------------------------------------------------------

#GARCH-t(1, 1) für die Subreihe nach der COVID-19 Krise
post_return <- Returns[(ne + 1):(n_cp - 1)]                                     # Die Nachkrisenrenditen
post_return.GARCH <- garchFit(~ garch(1, 1), post_return, cond.dist = "std", trace = FALSE)  # Anpassung des GARCH-t Modells

#Vergleich der GARCH Modelle 
pre_return.GARCH@fit$matcoef
post_return.GARCH@fit$matcoef

#Aufgabe 4 f)-------------------------------------------------------------------

retdif_plot1 <- ggplot(data.frame(x = Date[(bw + 1):(n_cp - bw+1)], y = RetDiff[1, ]), aes(x = x, y = y)) +
  geom_line() +
  ggtitle("Linksseitige lokal-gewichtete Renditedurchschnitte") +
  geom_vline(xintercept = Date[ns + 1], color = 2) +
  geom_vline(xintercept = Date[ne + 1], color = 7) +
  xlim(date_min, date_max) +
  xlab("Jahr") +
  ylab("") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
retdif_plot1

retdif_plot2 <- ggplot(data.frame(x = Date[(bw + 1):(n_cp - bw+1)], y = RetDiff[2, ]), aes(x = x, y = y)) +
  geom_line() +
  ggtitle("Differenzen der rechts- und linksseitigen Renditedurchschnitte") +
  geom_vline(xintercept = Date[ns + 1], color = 2) +
  geom_vline(xintercept = Date[ne + 1], color = 7) +
  xlim(date_min, date_max) +
  xlab("Jahr") +
  ylab("") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
retdif_plot2

retdif_plot3 <- ggplot(data.frame(x = Date[(bw + 1):(n_cp - bw+1)], y = RetDiff[3, ]), aes(x = x, y = y)) +
  geom_line() +
  ggtitle("Rechsseitige lokal-gewichtete Renditedurchschnitte") +
  geom_vline(xintercept = Date[ns + 1], color = 2) +
  geom_vline(xintercept = Date[ne + 1], color = 7) +
  xlim(date_min, date_max) +
  xlab("Jahr") +
  ylab("") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
retdif_plot3

retdif_comp <- ggarrange(retdif_plot2, retdif_plot3, retdif_plot1, ncol = 1, nrow = 3, align = "v")
retdif_comp
