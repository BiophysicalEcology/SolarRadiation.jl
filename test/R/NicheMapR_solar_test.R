library(NicheMapR) # load the NicheMapR package

longlat <- c(-89.40123, 43.07305)
EC <- 0.0167238 # Eccentricity of the earth's orbit (current value 0.0167238, ranges between 0.0034 to 0.058)
# Aerosol extinction coefficient profile
# the values extracted from GADS for Madison
TAI <- c(0.269904738,0.266147825, 0.262442906, 0.258789404, 0.255186744, 0.251634356, 0.248131676, 0.2412732, 0.234606887, 0.228128378, 0.221833385, 0.215717692, 0.20977715, 0.204007681, 0.198405272, 0.187685927, 0.177588357, 0.168082846, 0.159140695, 0.150734206, 0.142836655, 0.135422274, 0.128466227, 0.12194459, 0.115834329, 0.110113284, 0.104760141, 0.099754417, 0.09507644, 0.090707328, 0.086628967, 0.082823998, 0.07927579, 0.075968428, 0.072886691, 0.070016034, 0.067342571, 0.064853053, 0.062534858, 0.060375964, 0.058364941, 0.056490925, 0.054743609, 0.053113222, 0.051590514, 0.050166738, 0.046408775, 0.045302803, 0.044259051, 0.043271471, 0.042334415, 0.041442618, 0.040591184, 0.039775572, 0.038991583, 0.038235345, 0.037503301, 0.036792197, 0.036099067, 0.034101935, 0.033456388, 0.032817888, 0.032184949, 0.031556287, 0.030930816, 0.030307633, 0.029065372, 0.027825562, 0.027205981, 0.026586556, 0.025967391, 0.025348692, 0.024114005, 0.023498886, 0.021669152, 0.021066668, 0.019292088, 0.018144698, 0.016762709, 0.015451481, 0.014949794, 0.014224263, 0.013093462, 0.012670686, 0.012070223, 0.011164062, 0.010241734, 0.009731103, 0.009507687, 0.009212683, 0.008965785, 0.008827751, 0.008710756, 0.008574128, 0.008462605, 0.008446967, 0.008539475, 0.009015237, 0.009748444, 0.010586023, 0.011359647, 0.011901268, 0.012062153, 0.011735443, 0.010882215, 0.009561062, 0.007961182, 0.006438984, 0.005558204, 0.006133532, 0.009277754)
elev <- 226 # altitude (m)
slope <- 0 # slope (degrees, range 0-90)
aspect <- 180 # aspect (degrees, 0 = North, range 0-360)
hori <- rep(0, 24) # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
lamb <- 1 # Return wavelength-specific solar radiation output?
IUV <- 1 # Use gamma function for scattered solar radiation? (computationally intensive)
REFL <- 0.15 # soil reflectance
P_atmos <- get_pressure(elev)
input <- as.matrix(c(longlat, EC, elev, slope, aspect, lamb, IUV, REFL, P_atmos))
micro <- micro_global(loc = longlat, EC = EC, TAI = TAI, elev = elev, 
                      slope = slope, aspect = aspect, hori = hori, lamb = lamb, 
                      IUV = IUV, REFL = REFL, solonly = 1, runshade = 0,
                      clearsky = 1)

metout <- as.data.frame(micro$metout) # retrieve above ground microclimatic conditions, min shade
drlam <- as.data.frame(micro$drlam)
drrlam <- as.data.frame(micro$drrlam)
srlam <- as.data.frame(micro$srlam)

write.csv(input, file = 'c:/git/SolarRadiation.jl/test/data/input.csv')
write.csv(hori, file = 'c:/git/SolarRadiation.jl/test/data/horizon.csv')
write.csv(TAI, file = 'c:/git/SolarRadiation.jl/test/data/TAI.csv')
write.csv(metout$SOLR, file = 'c:/git/SolarRadiation.jl/test/data/global.csv')
write.csv(metout$ZEN, file = 'c:/git/SolarRadiation.jl/test/data/zenith.csv')
write.csv(drlam, file = 'c:/git/SolarRadiation.jl/test/data/drlam.csv')
write.csv(drrlam, file = 'c:/git/SolarRadiation.jl/test/data/drrlam.csv')
write.csv(srlam, file = 'c:/git/SolarRadiation.jl/test/data/srlam.csv')

plot(metout$SOLR, type = 'l')
