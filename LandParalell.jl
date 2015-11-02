#module Land
cd("/Users/Raukhur/Documents/Github/Landscapes")
using Landscapes

addprocs(4)
nprocs()
@everywhere cd("/Users/Raukhur/Documents/Github/Landscapes")
@everywhere using Landscapes, ODE, HDF5
@time @sync @parallel for i=1:6
  p=Landscapes.getP()
  p.inData="/Users/raukhur/Google Drive/Landscapes/inData"
  p.outData="/Users/raukhur/Google Drive/Landscapes/outData"
  p.NoLandscape=100#i
  Landscapes.GetData(p)
  p.Tend=3000
  p.Poisson=2
  out=Landscapes.Go(p)
  p.Poisson=0
  out=Landscapes.Go(p)
end

X=Landscapes.getResult(20,0,"/Users/raukhur/Google Drive/Landscapes/outData")
using Gadfly
p=Landscapes.getP()
plot(layer(x=p.z,y=X[end,:,10],Geom.point))


