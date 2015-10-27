cd("/Users/raukhur/Documents/Github/Landscapes")
using Landscapes
using Gadfly
using Distances
reload("Landscapes")



NoSpecies=100 # Set number of species
NoLandscape=1 #Select Landscape type
repl=1 # Select replicate of Landscape type
NoNoise=1 #Select noise series for replication
Tend=600 # Length of simulation
alpha=0.9 # Competition coefficient
Poisson=1 # 0=No colonization, 1=Poisson colonization, 2=Poisson colonization with no traversal of species
_save=0 #0: returns variables and does not save, 1: saves and gives no return

Landscapes.makeXYLandscape(1,1,-10,10000.,1)
Landscapes.makeXYLandscape(2,1,-10,100000.,1)
Landscapes.makeXYLandscape(3,1,100,10000.,1)
Landscapes.makeXYLandscape(4,1,100,100000.,1)


p1=Landscapes.getP(1,1,1,1,0.8)
X1=Landscapes.Go(p1,Tend,Poisson,_save)
p2=Landscapes.getP(2,1,1,1,0.8)
X2=Landscapes.Go(p2,Tend,Poisson,_save)
p3=Landscapes.getP(3,1,1,1,0.8)
X3=Landscapes.Go(p3,Tend,Poisson,_save)
p4=Landscapes.getP(4,1,1,1,0.8)
X4=Landscapes.Go(p4,Tend,Poisson,_save)

X=(X3[1])
p=p3
site=100
id=sortperm(p.XY[2,:][:])
plot(layer(x=p.z,y=X[end,:,id[site]],Geom.point))

p

p=p3
  p.dispersalAlpha=500
  p.seedPerBiomass=300*300*1000
  p.D=Landscapes.connectivity(p.XY,p.dispersalAlpha,p.dispersalC)
  D=pairwise(Euclidean(),p.XY,p.XY)
  P=zeros(size(p.D))
  N=100
  for i=1:N
    P+=Landscapes.PoissonRnd(p)
  end
  plot(layer(x=D[:]./1000, y=log10(1+p.D[:].*p.seedPerBiomass),Geom.line,Theme(default_color=color("orange"))),layer(x=(D[:]./1000),y=log10(1+P[:]/N),Geom.point))

plot(layer(x=D[:], y=p.D[:]*p.seedPerBiomass,Geom.point))

p

ISD[find(isnan(ISD))]=0.0
ntISD[find(isnan(ntISD))]=0.0
spy(squeeze(sum(ISD[:,:,id]-ntISD[:,:,id],1),1))


ntIE
spy(squeeze(X[:,:,id[2]],3))


spy(squeeze(X[end,:,id],1))

I=Landscapes.immigration(X[250,:,:],p,id[end])*p.reprod
Landscapes.PoissonImm(I[:],p,1,0)
methods(Landscapes.PoissonImm)
plot(x=p.z,y=Landscapes.immigration(X[250,:,:],p,id[end])*p.reprod)


# check default parameters (LandscapeNo::Int64,NoiseSeries::Int64,NoSpecies::Int64,alpha::Float64)
p=LandscapesP(LandscapeNo,replicate,NoiseSeries,alpha)
show(p)
p
# Check loading of landscapes
XY,N=Landscapes.GetData(p.NoLandscape,p.repl,p.NoNoise)
plot(x=XY[1,:],y=XY[2,:],Geom.point)
plot(x=1:length(N),y=p.noise,Geom.point,Theme(default_point_size=0.3mm,highlight_width=0mm))

# Create connectivity matrix
C=Landscapes.connectivity(p.XY,p.dispersalAlpha,p.dispersalC)
spy(C)

(typeof(1.0)==Int64 || typeof(1.0)==Float64)

#(NoSpecies::Int64,Landscape::Int64,NoiseSeries::Int64,Tend::Int64,extent::Float64,dispKernel::Float64,_Save::Int64)

Pkg.update()

maximum(IE)


Z=squeeze(M[:,:,2],3)
M=Landscapes.moments(X,p)
plot(layer(x=1:300,y=Z,Geom.point,Theme(default_point_size=0.3mm,highlight_width=0mm)))
size(X)
C=p.D
threshold=1e-8
C[find(C.<threshold)]=0.0
C
I=Landscapes.immigration(X,C,p,1)*p.reprod
plot(layer(x=p.z,y=I/maximum(I),Geom.line),layer(x=p.z,y=X[end,:,n]/maximum(X[end,:,n]),Geom.line, Theme(default_color=color("orange"))))


X=ones(Float64,2,NoSpecies,p.NoSites)
for i=1:p.NoSites
  for j=1:NoSpecies
    X[1,j,i]=exp(-(p.tempGrad[i]-p.z[j])^2/p.TWidth)
    X[2,j,i]=0.0
  end
end
X
plot(x=p.z,y=X[1,:,1])
SM,T,Dist,DS=addSouth(X,2,p) # CHECK DISTANCES AND ACTUAL IMMIGRATION FROM SOUTH!!
plot(x=p.z,y=DS,Geom.point)

T
DS
SM
Dist
p.z
x=ones(Float64,p.NoSpecies)*0.1/p.NoSpecies
dy=Landscapes.dy(0.0,x[:],20.0,p)
plot(x=p.z,y=dy,Geom.point)

p.extent
x=ones(p.NoSpecies)*0.5*p.overWinter/p.NoSpecies
tout,yout=Landscapes.sim(20.0,p,x,[0.0;180.0])
plot(x=p.z,y=yout[80])

p.XY
#test for speed of sim
x=ones(Float64,p.NoSpecies)*0.01
x=XCS[1,:,1]

X,XCS,p,timed=Landscapes.Go(100,1,2,100,0,0.01,0)
simtime=sum(timed)*p.NoSites
#(NoSpecies::Int64,Landscape::Int64,NoiseSeries::Int64,Tend::Int64,extent::Float64,dispKernel::Float64,_Save::Int64)
size(X)
b=@elapsed Landscapes.Go(100,1,2,100,0,0.01,0)
b*2/100*3000/60
# Overhead of Dispersal in percent
100*(b-simtime)/b

b*2/100*3000/60
LandscapeNo=82
rX,rXCS,rIE,rp,rtimed=Landscapes.Go(NoSpecies,LandscapeNo,repl,NoiseSeries,100,1,0.01,0)
X,XCS,IE,p,timed=Landscapes.Go(NoSpecies,LandscapeNo,repl,NoiseSeries,100,0,0.01,0)
XX=zeros(100,100)
for i=1:100
  for j=1:100
    f=sum(rIE[:,i,j],1)-sum(IE[:,i,j],1)
    XX[i,j]=f[1]#X[end,i,j]*abs(XCS[end,i,j]-1)#
  end
end
p.XY[2,:][:]
id=sortperm(p.XY[2,:][:])
XX

XX=XX[:,id]

spy(XX)


plot(x=p.z,y=-XX[:,7])

Pkg.add("Interact")
using Interact
Pkg.update()

@manipulate for i=1:100
  Gadfly.plot(x=p.z,y=XX[:,i])
end

