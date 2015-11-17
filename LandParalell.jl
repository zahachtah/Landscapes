#module Land
cd("/Users/Raukhur/Documents/Github/Landscapes")
push!(LOAD_PATH, "/Users/raukhur/Documents/Github/Landscapes")
using Landscapes
i=1

# one seed wheight 1 mg, and one plant 1 g dry mass
# 100 plants per m2 and 10 seeds per plant
# wetland 100*100*100 plants = 1 (K)
# one seed density 1/1000/(100*100)

i=1

addprocs(4)
nprocs()
#pp=Landscapes.getP()
#pp.inData="/Users/raukhur/Google Drive/Landscapes/inData"
#pp.outData="/Users/raukhur/Google Drive/Landscapes/outData"
@everywhere cd("/Users/Raukhur/Documents/Github/Landscapes")
@everywhere using Landscapes, ODE, HDF5
@time @sync @parallel for i=1:1
  out=Landscapes.Go(i,3000)
end


out=GetData(1, 1,1)

using Gadfly

#Questions:
# do species traverse?
# do species spread?
# Do we get the right south immigration?

reload("Landscapes")
p=Landscapes.getP()
p.NoLandscape=3
p.inData="/Users/Raukhur/Documents/Github/Landscapes/inData"
p.outData="/Users/Raukhur/Documents/Github/Landscapes/outData"
p=Landscapes.GetData(p)
N=4
f="/Users/Raukhur/Dropbox (Personal)/Wetlandmodelling/outData2"

x0,xcc0,Ta0=Landscapes.getResult(N,0,"/Users/Raukhur/Documents/Github/Landscapes/outData")
x1,xcc1,Ta1,SDX1=Landscapes.getResult(N,1,f)
x2,xcc2,Ta2=Landscapes.getResult(N,2,"/Users/Raukhur/Documents/Github/Landscapes/outData")
id=sortperm(p.XY[2,:][:])

x1

T=3000
  M=Landscapes.moments(x1,p)
  MS=Landscapes.moments(SDX1,p)
  SM,TT,Dist,DS=Landscapes.addSouth(x1,T,p)

plot(layer(x=p.z,y=DS*TT,Geom.line,Theme(default_color=color("orange"))),layer(x=p.z,y=x1[T,:,1],Geom.line),layer(x=p.z,y=x1[T,:,100],Geom.line))
plot(layer(x=p.XY[2,:],y=M[T,:,2]',Geom.point),layer(x=[0.0 p.XY[2,end]],y=[SM SM],Geom.line,Theme(default_color=color("green"))),layer(x=[0.0 p.XY[2,end]],y=[MS[T,1,2] MS[T,1,2]],Geom.line,Theme(default_color=color("orange"))))

1/Landscapes.dispersal(150000.0,5000.0,0.5)




Pkg.update()
