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
  LDDED::Float64
  LDDVD::Float64
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

"
Starts a full simulation based on parameter type definition. Generate a parameter type and redefine variables as needed like:

   ```
    p=getP()
    p.Tend=2000,
    ```

then run simulation as

   `Go(p)`

"
function Go(p::par)
  srand(1234+p.NoSpecies+p.NoLandscape*10+p.repl*100+p.NoNoise*1000+p.Tend*10000) # sets a random sequence that is different for all
  r=Float64
  LDDf=Float64
  LIE=Int64[]
  x=zeros(Float64,p.NoSpecies,1)+0.5/p.NoSpecies
  X=zeros(Float64,p.Tend,p.NoSpecies,p.NoSites)
  IE=zeros(Float64,p.Tend,p.NoSpecies,p.NoSites)
  ISD=zeros(Float64,p.Tend,p.NoSpecies,p.NoSites)
  XCS=ones(Float64,1,p.NoSpecies,p.NoSites)
  SDX=zeros(Float64,p.Tend,p.NoSpecies)
  I=Array(Float64,p.NoSites)
  SD=Array(Float64,p.NoSites)
  SM=Array(Float64,p.Tend)
  timed=zeros(Float64,p.Tend)
  Tactual=zeros(Float64,p.Tend,p.NoSites)
  X[1,:,:]=0.5/p.NoSpecies
  T=0.0
  for t=1:p.Tend
    # Get climate change
    TC=CC(t,p)
                
    # Get Distribution coming from south "mainland"
    SM[t],TotD,Dist,DS=southDist(X,t-1,p)
    SDX[t,:]=DS'*TotD

    # Remove traversing species if Poisson==2
    if t==p.CCstart && p.Poisson==2
      XCS[1,:,:]=mean(X[p.CCstart-26:p.CCstart-1,:,:],1); ## Better average over last stable 25 years !!
      XCS[find(XCS.>0.0)]=1.0;
    end

    for j=1:p.NoSites
      if t>1

        SD=addSouth(DS'*TotD,j,p)

        I=immigration(X[t-1,:,:].*XCS[end,:,:],p,j)+SD

        for k=1:p.NoSpecies
          if p.Poisson>=1
            if I[k]>0.0
              P=PoissonRnd(p::par,I[k]*p.seedPerBiomass)
              r=P*p.ext #p.ext=min propagule biomass
            else
              r=0.0
              P=0.0
            end
            IE[t-1,k,j]=P #save for return value to check immigration events
            ISD[t-1,k,j]=SD[k]

          elseif p.Poisson==0
            r=0.0
          else
            r=0.0
          end
                            
          ## Calculate between season change
          x[k]=X[t-1,k,j]*p.overWinter+r
          
          ## saves immigration events for later study
          if t>p.CCstart && X[t-1,k,j]==0.0 && x[k]>0.0
            push!(LIE, [t,j,k]...)
          end
                            
        end
                        
        ## Kill biomass less than minimum propagule size                 
        x[find(x.<p.ext)]=0.0
      end
      Tactual[t,j]=TC+p.noise[t]+p.tempGrad[j]
                    
      ## Calculate within season              
      tout,yout=simE(Tactual[t,j],p,x[:],[0.0;180.0])
      X[t,:,j]=yout
    end
  end
  M=moments(X,p)
  #file=p.outData*"/out"*string(p.NoLandscape)*"_"*string(p.Poisson)*"_"*string(p.repl)".h5"
  A=h5open(p.outData,"w")
     A["M","compress",3]=convert(Array{Float32,3},M)
     A["Tactual","compress",3]=convert(Array{Float32,2},Tactual)
     A["LIE","compress",3]=reshape(LIE,(3,div(length(LIE),3)))
     A["END","compress",3]=convert(Array{Float32,2},squeeze(X[end,:,:],1))
     #A["ISD","compress",3]=ISD
     #A["IE","compress",3]=IE
     #A["SDX","compress",3]=SDX
     #A["sDist","compress",3]=p.sDist
  close(A)
  return X,XCS,IE,ISD
end

function RunAll(i,inData,outData)
#Read M
A=h5open(inData*"/M.h5","r")
  M=read(A,"M")
close(A)


  p=Landscapes.getP()
  p.NoLandscape=M[i,4]
  println(M[i,:])
  p.Tend=3000
  p.Poisson=M[i,6]
  p.inData=inData
  p.outData=outData*"/"*string(i)*".h5"
  if M[i,1]==2
     p.r=0.01
     p.m=0.0005
  end
  if M[i,2]==1 #Set alpha
    p.alpha=0.8
  else
    p.alpha=0.5
  end

  p.repl=M[i,3] #set replicate
  p.NoNoise=M[i,3]
  if M[i,5]==2 #Set alpha
    p.LDDED=10000
    p.LDDVD=3000
  end
  p=Landscapes.GetData(p)
  
  return p, Go(p)
end

function sim(T::Float64,p::par,x::Array{Float64,1},tr::Array{Float64,1})
  tout,tout=ode23((t,x)->dy(t,x,T,p),x,tr)
  yout=tout[end][:]
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
  if p.LDDED>0.0
    LDDf=Float64
    LDDf=1/dispersal(p.LDDED,p.LDDVD,p.dispersalC)/p.seedPerBiomass
    p.D=p.D+LDDf*connectivity(p.XY,p.LDDVD,p.dispersalC)
  end
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

function dispersal(D::Array{Float64},a::Float64,c::Float64)
  d=zeros(Float64,length(D))
  for i=1:length(D)
  d[i]=dispersal(D[i],a,c)
  end
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

function immigrationLDD(X,p,site)
  LDDf=Float64
  if p.LDDED>0.0
    LDDf=1/dispersal(p.LDDED,p.LDDVD,p.dispersalC)/p.seedPerBiomass
  end
  I=zeros(Float64,1,p.NoSpecies)
  for i=1:p.NoSpecies
    for j=1:p.NoSites
      I[i]=I[i]+LDDf*dispersal(abs(p.XY[2,j]-p.XY[2,site]),p.LDDVD,p.dispersalC)*X[end,i,j]
    end
  end
  return I
end

function southDist(X,tt,p)
# Make sure the out variable is the distribution put on a 1:NoSPecies array at the right place
S=zeros(Float64,p.NoSites,1)
Dist=zeros(Float64,p.NoSpecies+1,1)
N=zeros(Float64,p.NoSpecies+1,1)
DS=zeros(Float64,p.NoSpecies,1)
SM=0.0
T=0.0
t=max(tt,1)
  for i=1:p.NoSites
    M=0.0
    M=p.z[indmax(X[t,:,i])]
    T=T+sum(X[t,:,i])
    S[i,:]=M
    id=indmax(X[t,:,i])
    if id<p.NoSpecies÷2
      q=X[t,1:p.NoSpecies-((p.NoSpecies÷2-id)*2+1),i]
    elseif id>p.NoSpecies÷2
      q=X[t,2*(id-p.NoSpecies÷2):p.NoSpecies,i]
    else
      q=zeros(Float64,1,p.NoSpecies+1)
      q[2:p.NoSpecies+1]=X[t,:,i]
    end
    dd=(p.NoSpecies-(length(q)-1))÷2
    Dist[1+dd:p.NoSpecies+1-dd]=Dist[1+dd:p.NoSpecies+1-dd]+q'./sum(q)
    N[1+dd:p.NoSpecies+1-dd]=N[1+dd:p.NoSpecies+1-dd]+1
  end
  Dist=Dist./N
  Dist[isnan(Dist)]=0.0 #some values that where zero get Nan, check
  T=T/p.NoSites
  reg=[p.XY[2,:];ones(Float64,p.NoSites,1)']'\S
  Xdist=collect(linspace(abs(p.sDist),100000.0,convert(Int64,(round(100000/abs(p.sDist))))))
  for j=1:length(Xdist)
    dp=dispersal(abs(Xdist[j]),p.dispersalAlpha,p.dispersalC)
    if p.LDDED>0.0
      LDDf=1/dispersal(p.LDDED,p.LDDVD,p.dispersalC)/p.seedPerBiomass
      dp=dp+LDDf*dispersal(abs(Xdist[j]),p.LDDVD,p.dispersalC)
    end
    if p.Poisson==0
        dp=0.0
    end
    SM=reg[2]+(p.sDist-Xdist[j])*reg[1]
    idm=indmin((p.z.-SM).^2)
    if p.z[max(1,idm)]>SM
      idm=max(1,idm-1)
    end
    for i=1:p.NoSpecies
      if idm-50+i>0 && idm-50+i<=p.NoSpecies
        DS[idm-50+i]=DS[idm-50+i]+dp*(Dist[i]*(SM-p.z[idm])+Dist[i+1]*(1-(SM-p.z[idm])))
      end
    end
  end
  if sum(DS)>0.0
      DS=DS/sum(DS)
  end
  return SM,T,Dist,DS
end

function addSouth(D,j,p)
  return D.*dispersal(abs(p.XY[2,j]-p.sDist),p.dispersalAlpha,p.dispersalC)
end

function addSouthLDD(D,j,p)
  LDDf=1/dispersal(p.LDDED,p.LDDVD,p.dispersalC)/p.seedPerBiomass
  return D.*LDDf*dispersal(abs(p.XY[2,j]-p.sDist),p.LDDVD,p.dispersalC)
end

end

