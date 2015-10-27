# Landscapes

Julia code to calculate response capacity in a metacommunity. Use by first calling:

``` julia

using Landscapes

```

then make a landscape:

```julia

Landscapes.makeXYLandscape(LandscapeNo,NoiseNo,NoSites,Extent,save)

```

* LandscapeNo creates a ID and folder so julia knows which to use. 
* NoiseNo uses a premade anual variation in temperature to be able to replicate this series
* NoSites is the number of sites in the landscape, if you use negative integer for this they will be alligned in a 1D gradient, e.g. -100
* Extent is the maximum distances in aither 2D or 1D
* save=1 saves the data, otherwise the function will return an array of XY values

You can create landscapes outside as long as they are stored in InData/Landscapes/L(number)/XY.h5 as a HDF5 file format with values for XY and Ext (Extent)

There are two ways to run the actual model:

``` julia

Landscapes.Go(NoSpecies,LandscapeNo,repl,NoiseNo,alpha,Tend,Poisson,save)

```

will make a run using the supplied values (repl is for landscape scenarios that have several replicates). Poisson=0 runs with no dispersal, Poisson=1 using dispersal and Poisson=2 uses dispersal but disallows traversal. The other way is to first get the input parameters, tweak them and then run the model:

``` julia

p=Landscapes.getP()
p.Tend=1200
Landscapes.Go(p)

```

For both ways, if you have save=1 they are saved as outData/out(landscapeNo)(repl).h5
