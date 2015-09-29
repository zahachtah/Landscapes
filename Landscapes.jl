module Landscapes
export Go,GetData,addSouth,dispersal,sim, dy, immigration
#=
!conn-matrix in julia
Check read in of landscapes
Scaling of landscapes
applying microtopology


=#
using ODE
using HDF5
using Distances
using Distributions
using ProgressMeter


immutable par
  NoSpecies::Int64
  NoSites::Int64
  NoLandscape::Int64
  repl::Int64
  NoNoise::Int64
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
  XY::Array{Float64,2}
  sDist::Float64
  D::Array{Float64,2}
  CCstart::Int64
  CCamp::Float64
  CCk::Float64
  tempSlope::Float64
  tempGrad::Array{Float64}
  noise::Array{Float64}
end

  function Base.show(io::IO, E::par)
    N=names(E)
    println()
    println("Simulation parameters")
    println("--------------------------------")
    for i in 1:length(N)
      if (typeof(E.(N[i]))==Int64 || typeof(E.(N[i]))==Float64)
        println(io,N[i]," (",E.(N[i]),")")
      else
        println(io,N[i]," (",typeof(E.(N[i])),")")
      end
    end
  end

include("LandscapeHelper.jl")

function dy(t::Float64,x::Array{Float64,1},T::Float64,p::par)
  dx=zeros(Float64,length(x))
  S::Float64
  ialpha::Float64
  S=sum(x)*p.alpha
  ialpha=1.0-p.alpha
  # @inbounds Add after debugging
  @inbounds for i=1:p.NoSpecies #adjust for intraspecific alpha
     dx[i]=((1.0-(S+(1.0-p.alpha)*x[i]))*exp(-(T-p.z[i])*(T-p.z[i])/p.TWidth)-p.m)*x[i]
  end
  return dx
end

function Go(NoSpecies::Int64,Landscape::Int64,repl::Int64,NoiseSeries::Int64,Tend::Int64,PoissonRand::Int64,_Save::Int64)
  #getP(LandscapeNo::Int64,NoiseSeries::Int64,NoSpecies::Int64,alpha::Float64)
  srand(1234+NoSpecies+Landscape*10+repl*100+NoiseSeries*1000+Tend*10000) # sets a random sequence that is different for all parameters except PoissonRand
  r=Float64
  p=getP(Landscape,repl,NoiseSeries,0.9)
  x=zeros(Float64,p.NoSpecies,1)+0.5/p.NoSpecies
  X=zeros(Float64,Tend,p.NoSpecies,p.NoSites)
  IE=zeros(Float64,Tend,p.NoSpecies,p.NoSites)
  ISD=zeros(Float64,Tend,p.NoSpecies,p.NoSites)
  XCS=ones(Float64,1,p.NoSpecies,p.NoSites)
  I=Array(Float64,p.NoSites)
  SD=Array(Float64,p.NoSites)
  timed=zeros(Float64,Tend)
  progress=Progress(Tend,1)
  X[1,:,:]=0.5/p.NoSpecies
  T=0.0
  for t=1:Tend

    if _Save==0 next!(progress) end
    TC=CC(t,p) # Get climate change
    SM,TotD,Dist,DS=addSouth(X,t-1,p)
    if t==p.CCstart && PoissonRand==2
      XCS[1,:,:]=X[max(1,t-1),:,:];
      XCS[find(XCS.>0.0)]=1.0;
      println("adjust XCS")
    end
    for j=1:p.NoSites
      if t>1
        SD=DS'*TotD*p.reprod*dispersal(p.XY[2,j]+p.sDist,p.dispersalAlpha,p.dispersalC)
        I=immigration(X[t-1,:,:].*XCS[end,:,:],p,j)*p.reprod+SD
        for k=1:p.NoSpecies
          if PoissonRand>=1
            if I[k]>0.0 #self immigration is zero
              P=PoissonRnd(p::par,10000.0*I[k])
              r=X[t-1,k,j]*(1-p.reprod)+P/p.seedPerBiomass
            else
              r=0.0
            end
            IE[t-1,k,j]=P #save for return value to check immigration events
            ISD[t-1,k,j]=SD[k]
          elseif PoissonRand==0
            r=0.0
          else
            r=X[t-1,k,j]*(1-p.reprod)*p.overWinter+I[j]
          end

          x[k]=X[t-1,k,j]*(1-p.reprod)*p.overWinter+r
          # remove call to dispersal out of loop
        end
      x[find(x.<p.ext)]=0.0
      end
      timed[t]=@elapsed tout,yout=sim(TC+p.noise[t]+p.tempGrad[j],p,x[:],[0.0;180.0])
      #println(yout[end][:])
      X[t,:,j]=yout[end][:]
    end
  end
  if _Save==1
    file="outData/out"*string(Landscape)*".h5"
    A=h5open(file,"w") do file
       write(file,"X", X)
       end
  else
    return X,XCS,IE,ISD,p,timed
  end

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

function PoissonRnd(p::par,I::Float64)
  P=0.0
  P=convert(Float64,rand(Poisson(I*p.seedPerBiomass)))
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
    M=sum(X[t,:,i].*p.z')/sum(X[t,:,i])
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
  SM=reg[2]+(1+p.sDist)*reg[1] #CHECK THAT THIS IS RIGHT
  Dist=Dist./N
  idm=indmin((p.z.-SM).^2)
  for i=1:idm
    if p.NoSpecies/2+1-(idm-i)>=1
      DS[i]=Dist[p.NoSpecies/2+1-(idm-i)]
    end
  end
  for i=idm+1:p.NoSpecies
    if i-idm<p.NoSpecies/2
      DS[i]=Dist[p.NoSpecies/2+i-idm+1]
    end
  end
  return SM,T,Dist,DS
end



end
