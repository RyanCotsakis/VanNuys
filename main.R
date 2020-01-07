library(evd)
library(evdbayes)
library(ismev)
library(tseries)
library(evir)

set.seed(42) # 42: The answer to everything

load("VanNuys.RData")


#
# --- DATA PROCESSING ---
#

date = seq(as.Date("1973-03-31"),as.Date("2015-05-31"),by="days") # Creates a date vector same length as wind and ffwi
wind = as.numeric(VanNuys[,1])
ffwi = as.numeric(VanNuys[,2])

month = months(date, abbreviate=FALSE)
month_as_numbers = match(month,month.name)
year = as.numeric(format(date,"%Y"))
year = year + (month_as_numbers > 7) # makes the year go from August to July to ensure that maxima occur near the middle of the annual period
day_of_year=as.numeric(format(date,"%j"))

plot(date,wind, main = "Wind Speeds - Time Series",
     xlab = "Year",
     ylab = "Wind Speed [mph]")
plot(date,ffwi, main = "FFWI - Time Series",
     xlab = "Year",
     ylab = "FFWI")

wind_measured = !is.na(wind) # Boolean vector: if wind was measured
ffwi_measured = !is.na(ffwi) # Boolean vector: if ffwi was measured

# Jitter
wind = jitter(wind,amount = 1.5)
ffwi = jitter(ffwi,amount = 1.5)
wind[!wind_measured] = 0
ffwi[!ffwi_measured] = 0


#
# --- UNIVARIATE WIND ---
#

# Annual Maxima Approach
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
plot(fit_GEV_wind)
fit_GEV_wind
plot(profile(fit_GEV_wind))

# POT Approach
th = 23
wind = as.numeric(VanNuys[,1])
wind = jitter(wind, amount = 0.6)
wind[!wind_measured] = 0
wind_pot_indices = wind_measured & year >= 1998 & (day_of_year < 120 | day_of_year > 320)
plot(day_of_year[wind > th],wind[wind>th] - th,
     main = "Wind Speeds Exceeding 23 mph - Seasonal Variation",
     xlab = "Day of the Year",
     ylab = "Threshold Exceedance [mph]")
plot(day_of_year[wind_measured & year >= 1998],wind[wind_measured & year >= 1998],
     main = "Wind Speeds - Seasonal Variation",
     xlab = "Day of the Year",
     ylab = "Wind Speed [mph]")
adf.test(wind[wind_pot_indices])

qu.min = quantile(wind[wind_pot_indices], 0.5)
qu.max = quantile(wind[wind_pot_indices],(length(wind[wind_pot_indices])-30)/length(wind[wind_pot_indices]))
mrlplot(wind[wind_pot_indices], tlim=c(qu.min, qu.max))
tcplot(wind[wind_pot_indices],tlim=c(qu.min, qu.max))
exiplot(wind[wind_pot_indices], tlim=c(20,30))

# Generate the approximate Retern Level Profiles
gpd.prof(gpd.fit(wind[wind_pot_indices],threshold =  23), m=10, xlow=34.5, xup=41, npy=164, conf=0.95)
title("Wind Speeds 10 Year Return Level Profile")
gpd.prof(gpd.fit(wind[wind_pot_indices],threshold =  23), m=100, xlow=38, xup=51.5, npy=164, conf=0.95)
title("Wind Speeds 100 Year Return Level Profile")

# Decluster to get more accurate parameters
exceedances_wind = pot(wind[wind_pot_indices],threshold=th)
declustered_wind = decluster(exceedances_wind$data,3,picture=FALSE)
fit_PP_declustered_wind = fpot(declustered_wind,threshold=th,model="pp")
plot(fit_PP_declustered_wind)
fit_PP_declustered_wind

#
# --- UNIVARIATE FFWI ---
#

# Annual Maxima Approach
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
fit_GEV_ffwi = fgev(ffwi_maxima,method="Nelder-Mead")
plot(fit_GEV_ffwi)
fit_GEV_ffwi
plot(profile(fit_GEV_ffwi))

# POT Approach
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

# Generate the approximate Retern Level Profiles
gpd.prof(gpd.fit(ffwi[ffwi_pot_indices],threshold =  70), m=10, xlow=95, xup=120, npy=164, conf=0.95)
title("FFWI 10 Year Return Level Profile")
gpd.prof(gpd.fit(ffwi[ffwi_pot_indices],threshold =  70), m=100, xlow=105, xup=160, npy=164, conf=0.95)
title("FFWI 100 Year Return Level Profile")

# Decluster
exceedances_ffwi = pot(ffwi[ffwi_pot_indices],threshold=th)
declustered_ffwi = decluster(exceedances_ffwi$data,3,picture=FALSE)
fit_PP_declustered_ffwi = fpot(declustered_ffwi ,threshold=th,model="pp")
plot(fit_PP_declustered_ffwi)
fit_PP_declustered_ffwi


##
## --- BIVARIATE ANALYSIS ---
##

# 3.1 - Maxima
bivariate_years = ffwi_maxima_years # There are fewer FFWI observations, because the time series starts later
wind_maxima = wind_maxima[-(1:8)] # Vectors need to be the same length
M = cbind(bivariate_years, ffwi_maxima, wind_maxima)

# One step approach
?fbvevd # Lists possible choices for model parameter
fit_bi_gev = fbvevd(M[,-1], model = "log")
plot(fit_bi_gev) # "which" argument tells which plots to generate
fit_bi_gev

# POT
q_ffwi = quantile(ffwi[ffwi_pot_indices], 0.95)
q_wind = quantile(wind[wind_pot_indices], 0.95)
bi_pot_indices = wind_pot_indices & ffwi_pot_indices
M_all = cbind(year, ffwi, wind)
fit_bi_pot = fbvpot(M_all[bi_pot_indices,-1], c(q_ffwi, q_wind), model = "log")
plot(fit_bi_pot)
fit_bi_pot

# EVIR
fit_bi_pot_evir = gpdbiv(ffwi[bi_pot_indices], wind[bi_pot_indices], q_ffwi, q_wind)
interpret.gpdbiv(fit_bi_pot_evir, q_ffwi, q_wind)
