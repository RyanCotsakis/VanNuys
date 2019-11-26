#Language setting, just in case, TODO delete
#library(readr)
#Sys.setlocale(locale="English")

library(evd)
library(evdbayes)

load("VanNuys.RData")

date = seq(as.Date("1973-03-31"),as.Date("2015-05-31"),by="days") # Creates a date vector
wind = as.numeric(VanNuys[,1])
ffwi = as.numeric(VanNuys[,2])

# day = seq(1,15402)
month=months(date, abbreviate=FALSE)
month_as_numbers = match(month,month.name)

plot(date,wind)
plot(date,ffwi)

plot(month_as_numbers,wind)
plot(month_as_numbers,ffwi)

wind_measured = !is.na(wind) # Boolean vector: if wind was measured
ffwi_measured = !is.na(ffwi) # Boolean vector: if ffwi was measured

# Auto Correlation Functions
acf(wind[wind_measured], lag.max = 400)
acf(ffwi[ffwi_measured], lag.max = 400)


#Fit GEV-wind
wind_maxima = apply(matrix(wind[wind_measured], ncol=365), 1, max)
fit_GEV_wind = fgev(wind_maxima,method="Nelder-Mead")
qqplot(qgumbel(c(1:42)/43),wind_maxima)
plot(fit_GEV_wind)
fit_GEV_wind
plot(profile(fit_GEV_wind))

#Fit GEV-ffwi
ffwi_maxima = apply(matrix(ffwi[ffwi_measured], ncol=365), 1, max)
fit_GEV_ffwi = fgev(ffwi_maxima,method="Nelder-Mead")
qqplot(qgumbel(c(1:42)/43),ffwi_maxima)
plot(fit_GEV_ffwi)
fit_GEV_ffwi
plot(profile(fit_GEV_ffwi))

#Stationarity
plot(apply(matrix(wind[wind_measured], ncol = 20), 1, mean))
plot(apply(matrix(ffwi[ffwi_measured], ncol = 120), 1, mean))

#fit GPD-ffwi
qu.min = quantile(ffwi[ffwi_measured], 0.5)
qu.max = quantile(ffwi[ffwi_measured],(length(ffwi[ffwi_measured])-30)/length(ffwi[ffwi_measured]))
mrlplot(ffwi[ffwi_measured], tlim=c(qu.min, qu.max))
tcplot(ffwi[ffwi_measured],tlim=c(qu.min, qu.max))

th = 80
fpot(dat,threshold=th)