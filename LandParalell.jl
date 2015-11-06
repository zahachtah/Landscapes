#module Land
cd("/Users/Raukhur/Documents/Github/Landscapes")
using Landscapes

addprocs(4)
nprocs()
#pp=Landscapes.getP()
#pp.inData="/Users/raukhur/Google Drive/Landscapes/inData"
#pp.outData="/Users/raukhur/Google Drive/Landscapes/outData"
@everywhere cd("/Users/Raukhur/Documents/Github/Landscapes")
@everywhere using Landscapes, ODE, HDF5
@time @sync @parallel for i=1:2
  out=Landscapes.Go(i,1200)
end


out=GetData(1, 1,1)

out[3]

using Gadfly

#Questions:
# do species traverse?
# do species spread?
# Do we get the right south immigration?

reload("Landscapes")
p=Landscapes.getP()
N=2
x0,xcc0,Ta0=Landscapes.getResult(N,0,"/Users/Raukhur/Documents/Github/Landscapes/outData")
x1,xcc1,Ta1=Landscapes.getResult(N,1,"/Users/Raukhur/Documents/Github/Landscapes/outData")
x2,xcc2,Ta2=Landscapes.getResult(N,2,"/Users/Raukhur/Documents/Github/Landscapes/outData")

Ta0
xc=xcc1[1,:,:]

M0=Landscapes.moments(x0,p)
S=1
plot(x=1:1200,y=Ta0[:,S]-M0[:,S,2],Geom.point)

plot(x=1:100,y=x2[end,:,10]-x1[end,:,10])

