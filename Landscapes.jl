module Landscapes
export Go,GetData,addSouth,dispersal,sim, dy, immigration

using ODE
using HDF5
using Distances
using Distributions
using ProgressMeter

type par
  NoSpecies::Int64
  NoSites::Int64
  NoLandscape::Int64
  Tend::Int64
  repl::Int64
  NoNoise::Int64
  Poisson::Int64
  r::Float64
  m::Float64
  ext::Float64
  alpha::Float64
  reprod::Float64
  dispersalAlpha::Float64
  dispersalC::Float64
  localTvar::Float64
  overWinter::Float64
  seedPerBiomass::Float64
  extent::Float64
  TWidth::Float64
  z::Array{Float64}
  XY::Array{Float64,2}      # Coordinates
  sDist::Float64            # added distance to southern pop
  D::Array{Float64,2}       # Connectivity probability to get between patches (calculated with dispersal function)
  CCstart::Int64
  CCamp::Float64
  CCk::Float64              # time of half CC after CCstart
  tempSlope::Float64        # degrees per km
  tempGrad::Array{Float64}  # Array of site starting temperatures
  noise::Array{Float64}     # Noise series uploaded from a file
  inData::ASCIIString       # Assign input folder path
  outData::ASCIIString
end

include("LandscapeHelper.jl")

function dy(t::Float64,x::Array{Float64,1},T::Float64,p::par)
  dx=zeros(Float64,length(x))
  S::Float64
  S=sum(x)*p.alpha
  # @inbounds Add after debugging
  @inbounds for i=1:p.NoSpecies #adjust for intraspecific alpha
     dx[i]=((1.0-(S+(1.0-p.alpha)*x[i]))*exp(-(T-p.z[i])*(T-p.z[i])/p.TWidth)-p.m)*x[i]
  end
  return dx
end

function simE(T::Float64,p::par,x::Array{Float64,1},tr::Array{Float64,1})
  for i in 0:180
    x=x.+dy(0.0,x,T,p)
  end
  return 0.0,x
end

function Go(L::Int64,T::Int64)
  p=Landscapes.getP()
  p.NoLandscape=L
  p.inData="/Users/Raukhur/Documents/Github/Landscapes/inData"
  p.outData="/Users/Raukhur/Documents/Github/Landscapes/outData"
  p=Landscapes.GetData(p)
  if L>4
    p.XY=p.XY*1000
    p.extent=p.extent*1000
  end
  p.seedPerBiomass=1/Landscapes.dispersal(10000.0,5000.0,0.5)
  p.dispersalAlpha=5000.0
  p.Tend=T
  out=Landscapes.Go(p)
  p.Poisson=2
  out=Landscapes.Go(p)
  p.Poisson=0
  out=Landscapes.Go(p)
end

function Go(p::par)
  srand(1234+p.NoSpecies+p.NoLandscape*10+p.repl*100+p.NoNoise*1000+p.Tend*10000) # sets a random sequence that is different for all
  r=Float64
  x=zeros(Float64,p.NoSpecies,1)+0.5/p.NoSpecies
  X=zeros(Float64,p.Tend,p.NoSpecies,p.NoSites)
  IE=zeros(Float64,p.Tend,p.NoSpecies,p.NoSites)
  ISD=zeros(Float64,p.Tend,p.NoSpecies,p.NoSites)
  XCS=ones(Float64,1,p.NoSpecies,p.NoSites)
  SDX=zeros(Float64,p.Tend,p.NoSpecies)
  I=Array(Float64,p.NoSites)
  SD=Array(Float64,p.NoSites)
  timed=zeros(Float64,p.Tend)
  progress=Progress(p.Tend,1)
  Tactual=zeros(Float64,p.Tend,p.NoSites)
  X[1,:,:]=0.5/p.NoSpecies
  T=0.0
  for t=1:p.Tend
    next!(progress)
    TC=CC(t,p) # Get climate change
    SM,TotD,Dist,DS=addSouth(X,t-1,p)
    SDX[t,:]=DS'*TotD
    if t==p.CCstart
      XCS[1,:,:]=X[max(1,t-1),:,:];
      XCS[find(XCS.>0.0)]=1.0;
      #("adjust XCS")
    end
    for j=1:p.NoSites
      if t>1
        SD=DS'*TotD*p.reprod*dispersal(abs(p.XY[2,j]-p.sDist),p.dispersalAlpha,p.dispersalC)
        if p.Poisson==2
          I=immigration(X[t-1,:,:].*XCS[end,:,:],p,j)*p.reprod+SD
        else
          I=immigration(X[t-1,:,:],p,j)*p.reprod+SD
        end
        for k=1:p.NoSpecies
          if p.Poisson>=1 || t<=p.CCstart
            if I[k]>0.0 #self immigration is zero
              # Maybe separate dispersal scaling parameter from seedBiomass
              P=PoissonRnd(p::par,I[k])
              r=X[t-1,k,j]*(1-p.reprod)+P*p.ext #p.ext=min propagule biomass
            else
              r=0.0
              P=0.0
            end
            IE[t-1,k,j]=P #save for return value to check immigration events
            ISD[t-1,k,j]=SD[k]
          elseif p.Poisson==0 && t>p.CCstart
            r=0.0
          else
            r=X[t-1,k,j]*(1-p.reprod)*p.overWinter+I[j]
          end

          x[k]=X[t-1,k,j]*(1-p.reprod)*p.overWinter+r
          # remove call to dispersal out of loop
        end
      x[find(x.<p.ext)]=0.0
      end
      Tactual[t,j]=TC+p.noise[t]+p.tempGrad[j]
      tout,yout=sim(Tactual[t,j],p,x[:],[0.0;180.0])
      X[t,:,j]=yout[end][:]
    end
  end
    file=p.outData*"/out"*string(p.NoLandscape)*"_"*string(p.Poisson)*".h5"
    A=h5open(file,"w") do file
       write(file,"X", X)
       write(file,"Tactual", Tactual)
       write(file,"XCS", XCS)
       write(file,"ISD", ISD)
       write(file,"IE", IE)
       write(file,"SDX", SDX)
      write(file,"sDist",p.sDist)
  end
  return X,XCS,IE,ISD
end

function Go(NoSpecies::Int64,Landscape::Int64,repl::Int64,NoiseSeries::Int64,alpha::Float64,Tend::Int64,PoissonRand::Int64,_Save::Int64)
  p=getP(Landscape,repl,NoiseSeries,Tend,alpha)
  p.Poisson=PoissonRand
  X,XCS,IE,ISD,p,timed=Go(p,_save)
  return X,XCS,IE,ISD,p,timed
end

function sim(T::Float64,p::par,x::Array{Float64,1},tr::Array{Float64,1})
  tout,yout=ode23((t,x)->dy(t,x,T,p),x,tr)
  return tout,yout
end

function CC(t,p)
  TC=p.CCamp*max(0.0,t-p.CCstart)^2/(p.CCk^2+max(0.0,t-p.CCstart)^2)
  return TC
end

function connectivity(XY,a,c)
  N=size(XY,2)
  C=zeros(Float64,N,N)
  D=pairwise(Euclidean(),XY,XY)
  for i=1:N
    for j=1:N
      C[i,j]=dispersal(D[i,j],a,c)
      if i==j
        C[i,j]=0
      end
    end
  end
  return C
end

function connectivity(p::par)
  p.D=connectivity(p.XY,p.dispersalAlpha,p.dispersalC)
  return p.D
end

function PoissonRnd(p::par,I::Float64)
  P=convert(Float64,rand(Poisson(I.*p.seedPerBiomass)))
  return P
end

function PoissonRnd(p::par)
  P=similar(p.D)
  for i=1:size(p.D,1)
    for j=1:size(p.D,2)
      P[i,j]=convert(Float64,rand(Poisson(p.D[i,j].*p.seedPerBiomass)))
    end
  end
  return P
end

function dispersal(D::Float64,a::Float64,c::Float64)
  d=(c/(2*pi*a^2*gamma(2/c))).*exp(-(D./a).^c)
  return d
end

function immigration(X,p,site)
  I=zeros(Float64,1,p.NoSpecies)
  for i=1:p.NoSpecies
    for j=1:p.NoSites
      I[i]=I[i]+p.D[j,site]*X[end,i,j]
    end
  end
  return I
end

function addSouth(X,tt,p)
# Make sure the out variable is the distribution put on a 1:NoSPecies array at the right place
S=zeros(Float64,p.NoSites,1)
Dist=zeros(Float64,p.NoSpecies+1,1)
N=zeros(Float64,p.NoSpecies+1,1)
DS=zeros(Float64,p.NoSpecies,1)
T=0.0
  t=max(tt,1)
  for i=1:p.NoSites
    M=0.0
    #M=sum(X[t,:,i].*p.z')/sum(X[t,:,i])
    M=p.z[indmax(X[t,:,i])]
    T=T+sum(X[t,:,i])
    S[i,:]=M
    id=indmax(X[t,:,i])
    if id<p.NoSpecies/2
      q=X[t,1:p.NoSpecies-((p.NoSpecies/2-id)*2+1),i]
    elseif id>p.NoSpecies/2
      q=X[t,2*(id-p.NoSpecies/2):p.NoSpecies,i]
    else
      q=zeros(Float64,1,p.NoSpecies+1)
      q[2:p.NoSpecies+1]=X[t,:,i]
    end
    dd=(p.NoSpecies-(length(q)-1))/2
    Dist[1+dd:p.NoSpecies+1-dd]=Dist[1+dd:p.NoSpecies+1-dd]+q'./sum(q)
    N[1+dd:p.NoSpecies+1-dd]=N[1+dd:p.NoSpecies+1-dd]+1
  end
  T=T/p.NoSites
  reg=[p.XY[2,:],ones(Float64,p.NoSites,1)']'\S
  SM=reg[2]+p.sDist*reg[1] #CHECK THAT THIS IS RIGHT
  Dist=Dist./N
  idm=indmin((p.z.-SM).^2)
  if p.z[idm]>SM
    idm=idm-1
  end
  for i=1:p.NoSpecies
    if idm-50+i>0 && idm-50+i<=p.NoSpecies
      DS[idm-50+i]=Dist[i]*(SM-p.z[idm])+Dist[i+1]*(1-(SM-p.z[idm]))
    end
  end
  #for i=1:idm
  #  if p.NoSpecies/2+1-(idm-i)>=1
  #    DS[i]=Dist[p.NoSpecies/2+1-(idm-i)]
  #  end
  #end
  #for i=idm+1:p.NoSpecies
  #  if i-idm<p.NoSpecies/2
  #    DS[i]=Dist[p.NoSpecies/2+i-idm+1]
  #  end
  #end
  return SM,T,Dist,DS
end

end
