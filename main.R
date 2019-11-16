load("VanNuys.RData")

date = seq(as.Date("1973-03-31"),as.Date("2015-05-31"),by="days") # Creates a date vector
wind = as.numeric(VanNuys[,1])
ffwi = as.numeric(VanNuys[,2])
# day = seq(1,15402)
# month=months(date)

# plot(date,wind)
# plot(date,ffwi)

wind_measured = !is.na(wind) # Boolean vector: if wind was measured
ffwi_measured = !is.na(ffwi) # Boolean vector: if ffwi was measured

# Auto Correlation Functions
acf(wind[wind_measured])
acf(ffwi[ffwi_measured])
