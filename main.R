#Language setting, just in case, TODO delete
#library(readr)
#Sys.setlocale(locale="English")

# TODO: Check dates of maxima to see if there is independence between years (eg Dec 31 to Jan 1)

library(evd)
library(evdbayes)
library(ismev)

load("VanNuys.RData")

date = seq(as.Date("1973-03-31"),as.Date("2015-05-31"),by="days") # Creates a date vector
wind = as.numeric(VanNuys[,1])
ffwi = as.numeric(VanNuys[,2])

# Analyze NAs
plot(date, is.na(ffwi))
plot(date, is.na(wind))

# day = seq(1,15402)
month=months(date, abbreviate=FALSE)
year=as.numeric(format(date,"%Y"))
day_of_year=as.numeric(format(date,"%j"))
month_as_numbers = match(month,month.name)
# 
# plot(date,wind)
# plot(date,ffwi)
# 
# plot(month_as_numbers,wind)
# plot(month_as_numbers,ffwi)

wind_measured = !is.na(wind) # Boolean vector: if wind was measured
ffwi_measured = !is.na(ffwi) # Boolean vector: if ffwi was measured
wind[!wind_measured] = 0
ffwi[!ffwi_measured] = 0

# Auto Correlation Functions
acf(wind[wind_measured], lag.max = 400)
acf(ffwi[ffwi_measured], lag.max = 400)
exiplot(ffwi[ffwi_measured], tlim=c(50,90))
# 
# # Get rid of ALL NA's
# wind = wind[wind_measured]
# ffwi = ffwi[ffwi_measured]
# 

#Fit GEV-wind
wind_maxima = array(0,43)
for(j in 1:43){
  wind_maxima[j]=max(wind[year==1972+j])
}
fit_GEV_wind = fgev(wind_maxima,method="Nelder-Mead")
qqplot(qgumbel(c(1:42)/43),wind_maxima)
plot(fit_GEV_wind)
fit_GEV_wind
plot(profile(fit_GEV_wind))

#Fit GEV-ffwi
ffwi_maxima = array(0,43)
for(j in 1:43){
  ffwi_maxima[j]= max(ffwi[year == 1972 +j])
}
fit_GEV_ffwi = fgev(ffwi_maxima,method="Nelder-Mead")
qqplot(qgumbel(c(1:42)/43),ffwi_maxima)
plot(fit_GEV_ffwi)
fit_GEV_wind
plot(profile(fit_GEV_wind))

## R largest:
r = 9
ffwi_biggest = matrix(0,nrow=43,ncol=10)
for(j in 1:43){
  ffwi_biggest[j,1:10] = sort(ffwi[year==1972+j], decreasing = TRUE)[1:10]
}
ffwi_biggest = ffwi_biggest[9:42 , 1:r]
ffwi_biggest.vec = as.vector(t(ffwi_biggest))
year.vec = rep(1:34, each=r)
fit1 = rlarg.fit(xdat = ffwi_biggest)
par(mfrow = c(3,3))
rlarg.diag(fit1)

#gumble profile
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

th = 62 #TODO: check threshold assumptions 
fit_GPD_ffwi = fpot(ffwi[ffwi_measured],threshold=th)

# check fit
plot(fit_GPD_ffwi) #prob. plot wiggling out, QQ wiggling, density ok, return wiggling
#wiggling because of discreteness, step distance 3 dealt with in next chapter

#exceedence rate
plot(day_of_year[ffwi>th],ffwi[ffwi>th])
exceedence_per_day = array(0,366)
for(j in 1:366){
  numerator = sum(ffwi[day_of_year == 3 && ffwi_measured])
  exceedence_per_day[j] = sum(ffwi[day_of_year==j & ffwi>th])/sum(ffwi[day_of_year==j])
}


# profile log-likelihood for both parameters
par(mfrow=c(1,2))
plot(profile(fit_GPD_ffwi))
abline(v=0,col=2,lty=2) # 0 not inside of the confidence intervals
par(mfrow=c(1,1))

# fit shape=0
fit_GPD_ffwi_gum = fpot(ffwi[ffwi_measured], threshold=th, shape=0)
plot(fit_GPD_ffwi_gum) # almost same, maybe a little worse

# deal with varying variance
plot(date[ffwi>th],ffwi[ffwi>th])
fit_GPD_ffwi_linkfit = gpd.fit(ffwi[ffwi_measured],threshold=th,ydat=as.matrix(scale(date)),sigl=1,siglink=exp)
gpd.diag(fit_GPD_ffwi_linkfit) #TODO: ??

##
## eliminating bad data & discreteness
# create elimination vectors
good_dates = (year>1992) & (ffwi_measured) & ((day_of_year<150) | (day_of_year>280))
plot(good_dates)

#eliminate discreteness
ffwi_jittered = jitter(ffwi,amount = 1.5)
ffwi=ffwi_jittered

# choose threshold, fit GPD for good dates and check fit
qu.min = quantile(ffwi[good_dates], 0.5)
qu.max = quantile(ffwi[good_dates],(length(ffwi[good_dates])-30)/length(ffwi[good_dates]))
mrlplot(ffwi[good_dates], tlim=c(qu.min, qu.max))
tcplot(ffwi[good_dates],tlim=c(qu.min, qu.max))

th = 55 #TODO: check threshold assumptions 
fit_GPD_ffwi = fpot(ffwi[good_dates],threshold=th)

# check fit
plot(fit_GPD_ffwi) #prob. plot wiggling out, QQ wiggling, density ok, return wiggling


##
# fit PP-ffwi
#choose threshold as before, compare with...
u = quantile(ffwi[ffwi_measured],0.95)
fit_PP_ffwi = fpot(ffwi[good_dates],threshold=th,model="pp")
plot(fit_PP_ffwi) #same as before

##
## --- CHAPTER 3: BIVARIATE ANALYSIS ---
##

# 3.1 - Maxima
y = c(1973:2015)
M = cbind(y, ffwi_maxima, wind_maxima)
M = M[9:42,]
M_all = cbind(year, ffwi, wind)

fit_bi_gev = fbvevd(M[,-1], model = "log") # TODO: try different values of model
plot(fit_bi_gev, mar = 2) # "which" argument tells which plots to generate
plot(fit_bi_gev)

q_ffwi = quantile(ffwi, 0.95)
q_wind = quantile(wind, 0.95)
fit_bi_pot = fbvpot(M_all[,-1], c(q_ffwi, q_wind), model = "log")
