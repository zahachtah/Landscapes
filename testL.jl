cd("/Users/raukhur/Google Drive/juliabox/Landscapes")
using Landscapes
using Gadfly
reload("Landscapes")

# NEED TO CHECK TEMP GRADIENT PER 500 KM!!!!!!!!

NoSpecies=100 # Set number of species
LandscapeNo=1 #Select Landscape type
repl=1 # Select replicate of Landscape type
NoiseSeries=1 #Select noise series for replication
Tend=600 # Length of simulation
alpha=0.9 # Competition coefficient
Poisson=1 # 0=No colonization, 1=Poisson colonization, 2=Poisson colonization with no traversal of species
_save=0 #0: returns variables and does not save, 1: saves and gives no return


X,XCS,IE,ISD,p,tim=Landscapes.Go(NoSpecies,LandscapeNo,repl,NoiseSeries,Tend,Poisson,_save)
nX,nXCS,nIE,nISD,np,ntim=Landscapes.Go(NoSpecies,LandscapeNo,repl,NoiseSeries,Tend,0,_save)
ntX,ntXCS,ntIE,ntISD,ntp,nttim=Landscapes.Go(NoSpecies,LandscapeNo,repl,NoiseSeries,Tend,2,_save)

site=45
id=sortperm(p.XY[2,:][:])
plot(layer(x=p.z,y=X[end,:,id[site]],Geom.point),layer(x=p.z,y=nX[end,:,id[site]],Geom.point,Theme(default_color=color("orange"))),layer(x=p.z,y=ntX[end,:,id[site]],Geom.point,Theme(default_color=color("green"))))

Landscapes.PoissonRnd(p,10000.0*p.reprod*Landscapes.dispersal(50.0,p.dispersalAlpha,p.dispersalC))

Landscapes.PoissonRand(p,0.1)
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

