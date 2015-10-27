# Landscapes

Julia code to calculate response capacity in a metacommunity. Use by first calling:

``` julia

using Landscapes

```

then make a landscape:

```julia

Landscapes.makeXYLandscape(LandscapeNo,NoiseNo,NoSites,save)

```

* LandscapeNo creates a ID and folder so julia knows which to use. 
* NoiseNo uses a premade anual variation in temperature to be able to replicate this series
* NoSites is the number of sites int he landscape, if you use negative integer for this they will be alligned in a 1D gradient
* save=1 saves the data, otherwise the function will return an array of XY values

There are two ways to run the actual model:

``` julia
Landscapes.Go(NoSpecies,LandscapeNo,repl,NoiseNo,alpha,Tend,Poisson,save)

```

will make a run using the supplied values (repl is for landscape scenarios that have several replicates). Poisson=0 runs with no dispersal, Poisson=1 using dispersal and Poisson=2 uses dispersal but disallows traversal. The other way is to first get the input parameters adn then run the model:

``` julia

p=Landscapes.getP()
p.Tend=1200
Landscapes.Go(p)

```

