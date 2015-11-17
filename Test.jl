cd("/Users/raukhur/Documents/Github/Landscapes")

push!(LOAD_PATH, "/Users/raukhur/Documents/Github/Landscapes")
using Landscapes
using Gadfly
Pkg.clone()
include("Landscapes.jl")
reload("Landscapes")
p=Landscapes.getP()
p.inData="/Users/raukhur/Google Drive/Landscapes/inData"
p.outData="/Users/raukhur/Google Drive/Landscapes/outData"
p.NoLandscape=113
Landscapes.GetData(p)
N=300
plot(x=1:10,y=rand(10))
plot(x=p.XY[1,:],y=p.XY[2,:])
out=Landscapes.Go(p)
pwd()
readdir()
p.tempGrad
Pkg.add("ProgressMeter")
getP()
pwd()
readdir()
p
LOAD_PATH
