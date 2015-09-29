#module Land
addprocs(7)
nprocs()
@everywhere cd("/home/x_jonno/Landscapes")
@everywhere using Landscapes, ODE, HDF5
@time @sync @parallel for i=1:14
    Landscapes.Go(100,i,3,50,1)
end

cd("/Users/Raukhur/Google Drive/juliabox/Landscapes")
cd("/Users/jon/Google Drive/juliabox/Landscapes")
using Landscapes
reload("Landscapes")
X,p=Landscapes.Go(100,1,2,400,0)
@time Landscapes.Go(100,1,2,200,0)
120/200*2*3000/60
# try to test calling sim 200 times

using Gadfly
plot(x=linspace(-20,20,101),y=DS.*T)

plot(x=linspace(-20,20,101),y=Dist.*T)

plot(x=linspace(0,400,1200),y=X[:,:,30],Theme(default_point_size=0.3mm,highlight_width=0mm))

randn(100,1)*0.5

D=GetData()

XY=Landscapes.GetData(1)
SM,T,Dist,DS=addSouth(X,p)
XY=Landscapes.makeXYLandscape(1,100,0)
Landscapes.makeXYLandscape(2,50,500.0,1)
Landscapes.makeNoise(2,10000,1)

southD=0.0
Landscapes.dispersal(rand(p.NoSpecies,1)[1],0.1)
62/200*3000/60

X[end,:,10]

@time(for i=1:200 tout,yout=Landscapes.sim(0+p.XY[2,10]*10.0-5.0,p,rand(p.NoSpecies,1)*0.01,[0.0;180.0]) end)

typeof(yout)

yout[1]

