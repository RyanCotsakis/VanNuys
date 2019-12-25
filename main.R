#Language setting, just in case, TODO delete
#library(readr)
#Sys.setlocale(locale="English")

# TODO: Check dates of maxima to see if there is independence between years (eg Dec 31 to Jan 1)

library(evd)
library(evdbayes)
library(ismev)
library(tseries)

set.seed(42) # 42: The answer to everything

load("VanNuys.RData")

date = seq(as.Date("1973-03-31"),as.Date("2015-05-31"),by="days") # Creates a date vector
wind = as.numeric(VanNuys[,1])
ffwi = as.numeric(VanNuys[,2])

month=months(date, abbreviate=FALSE)
month_as_numbers = match(month,month.name)
year=as.numeric(format(date,"%Y"))
year = year + (month_as_numbers > 7)
day_of_year=as.numeric(format(date,"%j"))

plot(date,wind, main = "Wind Speeds - Time Series",
     xlab = "Year",
     ylab = "Wind Speed [mph]")
plot(date,ffwi, main = "FFWI - Time Series",
     xlab = "Year",
     ylab = "FFWI")

# plot(month_as_numbers,wind)
# plot(month_as_numbers,ffwi)

wind_measured = !is.na(wind) # Boolean vector: if wind was measured
ffwi_measured = !is.na(ffwi) # Boolean vector: if ffwi was measured

# Jitter
wind = jitter(wind,amount = 1.5)
ffwi = jitter(ffwi,amount = 1.5)
wind[!wind_measured] = 0
ffwi[!ffwi_measured] = 0
plot(year[wind_measured], wind[wind_measured])
plot(year[ffwi_measured], ffwi[ffwi_measured])

# Auto Correlation Functions
# acf(wind[wind_measured], lag.max = 400)
# acf(ffwi[ffwi_measured], lag.max = 400)


#
# --- UNIVARIATE WIND ---
#
wind_maxima = array(0,43) # Initialize array
for(j in 1:43){
  wind_maxima[j]=max(wind[year==1972+j])
}
wind_maxima # Print values
wind_maxima = wind_maxima[-15] #remove year 1987 due to many NA's
wind_maxima_years = c(1973:2015)
wind_maxima_years = wind_maxima_years[-15]
plot(wind_maxima_years, wind_maxima, main = "Wind Speeds - Anual Maxima",
     xlab = "Year",
     ylab = "Wind Speed [mph]")

fit_GEV_wind = fgev(wind_maxima,method="Nelder-Mead")
qqplot(qgumbel(c(1:42)/43),wind_maxima)
plot(fit_GEV_wind)
fit_GEV_wind
plot(profile(fit_GEV_wind))

# Get new indices for POT
th = 23
wind = as.numeric(VanNuys[,1])
wind = jitter(wind, amount = 0.6)
wind[!wind_measured] = 0
wind_pot_indices = wind_measured & year >= 1998 & (day_of_year < 120 | day_of_year > 320)
plot(day_of_year[wind > th],wind[wind>th] - th,
     main = "Wind Speeds Exceeding 23 mph - Seasonal Variation",
     xlab = "Day of the Year",
     ylab = "Threshold Exceedance [mph]")
adf.test(wind[wind_pot_indices])

qu.min = quantile(wind[wind_pot_indices], 0.5)
qu.max = quantile(wind[wind_pot_indices],(length(wind[wind_pot_indices])-30)/length(wind[wind_pot_indices]))
mrlplot(wind[wind_pot_indices], tlim=c(qu.min, qu.max))
tcplot(wind[wind_pot_indices],tlim=c(qu.min, qu.max))
exiplot(wind[wind_pot_indices], tlim=c(20,30))

fit_PP_wind = fpot(wind[wind_pot_indices],threshold=th,model="pp")
plot(fit_PP_wind)
fit_PP_wind

fit_PP_wind_gumble = fpot(wind[wind_pot_indices],threshold=th,model="pp", shape=0)
plot(fit_PP_wind_gumble)
fit_PP_wind_gumble


#
# --- UNIVARIATE FFWI ---
#
# TODO: check which quantile we dont exceed until ~1992
ffwi_maxima = array(0,43) # Initialize array
for(j in 1:43){
  ffwi_maxima[j]= max(ffwi[year == 1972 +j])
}
ffwi_maxima # Print values
ffwi_maxima = ffwi_maxima[-15]
ffwi_maxima = ffwi_maxima[-(1:8)]
ffwi_maxima_years = c(1973:2015)
ffwi_maxima_years = ffwi_maxima_years[-15]
ffwi_maxima_years = ffwi_maxima_years[-(1:8)]
plot(ffwi_maxima_years,ffwi_maxima, main = "FFWI - Anual Maxima",
     xlab = "Year",
     ylab = "FFWI")
qqplot(qgumbel(c(1:42)/43),ffwi_maxima)
fit_GEV_ffwi = fgev(ffwi_maxima,method="Nelder-Mead")
plot(fit_GEV_ffwi)
fit_GEV_ffwi
plot(profile(fit_GEV_ffwi))

# Get new indices for POT
th = 70
plot(day_of_year[ffwi > th],ffwi[ffwi>th]-th, main = "FFWIs Exceeding 70 - Seasonal Variation",
     xlab = "Day of the Year",
     ylab = "Threshold Exceedance")
ffwi_pot_indices = ffwi_measured & year >= 1998 & (day_of_year < 120 | day_of_year > 320)
adf.test(ffwi[ffwi_pot_indices])

qu.min = quantile(ffwi[ffwi_pot_indices], 0.5)
qu.max = quantile(ffwi[ffwi_pot_indices],(length(ffwi[ffwi_pot_indices])-30)/length(ffwi[ffwi_pot_indices]))
mrlplot(ffwi[ffwi_pot_indices], tlim=c(qu.min, qu.max))
tcplot(ffwi[ffwi_pot_indices],tlim=c(qu.min, qu.max))
exiplot(ffwi[ffwi_pot_indices], tlim=c(40,80))

fit_PP_ffwi = fpot(ffwi[ffwi_pot_indices],threshold=th,model="pp")
plot(fit_PP_ffwi)
fit_PP_ffwi

fit_PP_ffwi_gumble = fpot(ffwi[ffwi_pot_indices],threshold=th,model="pp", shape=0)
plot(fit_PP_ffwi_gumble)
fit_PP_ffwi_gumble

## R largest:
# r = 9
# ffwi_biggest = matrix(0,nrow=43,ncol=10)
# for(j in 1:43){
#   ffwi_biggest[j,1:10] = sort(ffwi[year==1972+j], decreasing = TRUE)[1:10]
# }
# ffwi_biggest = ffwi_biggest[9:42 , 1:r]
# ffwi_biggest.vec = as.vector(t(ffwi_biggest))
# year.vec = rep(1:34, each=r)
# fit1 = rlarg.fit(xdat = ffwi_biggest)
# par(mfrow = c(3,3))
# rlarg.diag(fit1)


# The following was replicated in the code for some reason...
# {
# #gumble profile 
# fit_GEV_ffwi = fgev(ffwi_maxima,method="Nelder-Mead")
# qqplot(qgumbel(c(1:42)/43),ffwi_maxima)
# plot(fit_GEV_ffwi)
# fit_GEV_ffwi
# plot(profile(fit_GEV_ffwi))
# }


#Stationarity
# plot(apply(matrix(wind[wind_measured], ncol = 20), 1, mean))
# plot(apply(matrix(ffwi[ffwi_measured], ncol = 120), 1, mean))


#exceedence rate
# exceedence_per_day = array(0,366)
# for(j in 1:366){
#   numerator = sum(ffwi[day_of_year == 3 && ffwi_measured])
#   exceedence_per_day[j] = length(ffwi[day_of_year==j & ffwi>th])/length(ffwi[day_of_year==j])
# }
# plot(c(1:366),exceedence_per_day)

# abline(v=0,col=2,lty=2) # 0 not inside of the confidence intervals

# deal with varying variance
# plot(date[ffwi>th],ffwi[ffwi>th])
# fit_GPD_ffwi_linkfit = gpd.fit(ffwi[ffwi_measured],threshold=th,ydat=as.matrix(scale(date)),sigl=1,siglink=exp)
# gpd.diag(fit_GPD_ffwi_linkfit) #TODO: ??

# 
# ## eliminating bad data & discreteness
# # create elimination vectors
# good_dates = (year>1992) & (ffwi_measured) & ((day_of_year<150) | (day_of_year>280))
# plot(good_dates)
# 
# #eliminate discreteness
# # ffwi_jittered = jitter(ffwi,amount = 1.5)
# # ffwi=ffwi_jittered
# 
# # choose threshold, fit GPD for good dates and check fit
# qu.min = quantile(ffwi[good_dates], 0.5)
# qu.max = quantile(ffwi[good_dates],(length(ffwi[good_dates])-30)/length(ffwi[good_dates]))
# mrlplot(ffwi[good_dates], tlim=c(qu.min, qu.max))
# tcplot(ffwi[good_dates],tlim=c(qu.min, qu.max))
# 
# th = 55 #TODO: check threshold assumptions
# fit_GPD_ffwi = fpot(ffwi[good_dates],threshold=th)
# 
# # check fit
# plot(fit_GPD_ffwi) #prob. plot wiggling out, QQ wiggling, density ok, return wiggling
# 

##
# fit PP-ffwi
#choose threshold as before, compare with...
# # u = quantile(ffwi[ffwi_measured],0.95)
# fit_PP_ffwi = fpot(ffwi[good_dates],threshold=th,model="pp")
# # plot(fit_PP_ffwi) #same as before

##
## --- CHAPTER 3: BIVARIATE ANALYSIS ---
##

# 3.1 - Maxima
bivariate_years = ffwi_maxima_years
wind_maxima = wind_maxima[-(1:8)]
M = cbind(bivariate_years, ffwi_maxima, wind_maxima)

# One step approach
?fbvevd
fit_bi_gev = fbvevd(M[,-1], model = "log") # TODO: try different values of model
# See which models give asymmetry ==> BAD!
plot(fit_bi_gev) # "which" argument tells which plots to generate
fit_bi_gev

# Two step approach
fit_bi_frechet.mar1 = fgev(ffwi_maxima)
fit_bi_frechet.mar2 = fgev(wind_maxima)
rescaled_ffwi_maxima = qgev(pgev(ffwi_maxima, loc=fit_bi_frechet.mar1$param["loc"], scale=fit_bi_frechet.mar1$param["scale"],
        shape=fit_bi_frechet.mar1$param["shape"]), loc=1, scale=1, shape=1)
rescaled_wind_maxima = qgev(pgev(wind_maxima, loc=fit_bi_frechet.mar2$param["loc"], scale=fit_bi_frechet.mar2$param["scale"],
        shape=fit_bi_frechet.mar2$param["shape"]), loc=1, scale=1, shape=1)
fit_bi_frechet = fbvevd(cbind(rescaled_ffwi_maxima, rescaled_wind_maxima), cscale = TRUE, cshape = TRUE, cloc = TRUE,
        loc1 = 1, scale1 = 1, shape1 = 1)
plot(fit_bi_frechet)
fit_bi_frechet

# 3.2 - POT
q_ffwi = quantile(ffwi[ffwi_pot_indices], 0.95) # Seems to be okay to use 0.95. (stefano did that)
q_wind = quantile(wind[wind_pot_indices], 0.95)
bi_pot_indices = wind_pot_indices & ffwi_pot_indices
M_all = cbind(year, ffwi, wind)
fit_bi_pot = fbvpot(M_all[bi_pot_indices,-1], c(q_ffwi, q_wind), model = "log")
plot(fit_bi_pot)
fit_bi_pot

# EVIR
library(evir)
fit_bi_pot_evir = gpdbiv(ffwi[bi_pot_indices], wind[bi_pot_indices], q_ffwi, q_wind)
interpret.gpdbiv(fit_bi_pot_evir, q_ffwi, q_wind)
