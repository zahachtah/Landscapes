# Landscapes

Julia code to calculate response capacity in a metacommunity. Use by first calling:

``` julia

using Landscapes

```

then make a landscape:

```julia

Landscapes.makeXYLandscape(LandscapeNo,NoiseNo,NoSites,save)

```

LandscapeNo creates a ID and folder so julia knows which to use. 
NoiseNo uses a premade anual variation in temperature to be able to replicate this series
NoSites is the number of sites int he landscape, if you use negative integer for this they will be alligned in a 1D gradient
save=1 saves the data, otherwise the function will return an array of XY values


First one needs to create a landscape
