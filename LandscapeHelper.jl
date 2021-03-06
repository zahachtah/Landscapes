function GetData(LandscapeNo, repl,NoiseSeries)
  file="inData/Landscapes/L"*string(LandscapeNo)*"/XY"*string(repl)*".h5"
  XY=h5open(file,"r") do file
       read(file,"XY")
       end
  Ext=h5open(file,"r") do file
       read(file,"Extent")
  end
    file="inData/Noise/N"*string(NoiseSeries)*".h5"
    N=h5open(file,"r") do file
       read(file,"N")
       end
  return XY,N,Ext
end

function GetData(p::par)
  file=p.inData*"/Landscapes/L"*string(p.NoLandscape)*"/XY"*string(p.repl)*".h5"
  p.XY=h5open(file,"r") do file
       read(file,"XY")
       end
  p.NoSites=size(p.XY,2)
  Ext=h5open(file,"r") do file
       read(file,"Extent")
  end
   p.extent=Ext[1]
    file=p.inData*"/Noise/N"*string(p.NoNoise)*".h5"
    p.noise=h5open(file,"r") do file
       read(file,"N")
       end
  p.tempGrad=20-((p.XY[2,:]/p.extent)*2-1.0)*p.tempSlope*p.extent
  connectivity(p)
  D=pairwise(Euclidean(),p.XY,p.XY)
  D=D+eye(D)*maximum(D)
  p.sDist=minimum(p.XY[2,:])-mean(minimum(D,1))
  return p
end

function makeXYLandscape(LandscapeNo::Int64,repl::Int64,NoSites::Int64,Extent::Float64,_save::Int64)
    if NoSites>0
        XY=rand(Float64,2,NoSites).*Extent #XY coordinates in a 0:1 range
    else
        NoSites=abs(NoSites)
        XY=[zeros(NoSites)'; collect(linspace(0,Extent,NoSites))'] # for making 1D landscapes

    end
  if _save==1
        if isdir("inData/Landscapes/L"*string(LandscapeNo))
    else
      mkdir("inData/Landscapes/L"*string(LandscapeNo))
    end
    file="inData/Landscapes/L"*string(LandscapeNo)*"/XY"*string(repl)*".h5"
    A=h5open(file,"w") do file
      write(file,"XY", XY)
      write(file,"Extent", Extent)
       end
  else
    return XY
  end
end

function makeNoise(NoiseNo::Int64,Tend::Int64,_save::Int64)
    #  No evidence of year to year autocorrelation
    #  http://journals.ametsoc.org/doi/pdf/10.1175/1520-0493%281977%29105%3C0009%3AEOTAAS%3E2.0.CO%3B2
    #  Generally normal distributed year to year anomalies
    #  http://www.giss.nasa.gov/research/briefs/hansen_17/
    #  average of 0.5 st.d.
    #  http://www.nature.com/nature/journal/v500/n7462/pdf/nature12310.pdf
  N=randn(1,Tend)*0.5 #
  if _save==1
    if isdir("inData/Noise")
    else
      mkdir("inData/Noise")
    end
    file="inData/Noise/N"*string(NoiseNo)*".h5"
    A=h5open(file,"w") do file
      write(file,"N", N)
       end
  else
    return N
  end
end

function getP()
  getP(1,1,1,600,1,0.8)
end

function getP(NoLandscape::Int64,repl::Int64,NoNoise::Int64,Tend::Int64,Poisson::Int64,alpha::Float64)
    NoSpecies=100
    XY,N,extent=GetData(NoLandscape,repl,NoNoise)
    NoSites=size(XY,2)
    r=0.1
    m=0.05*r
    ext=1.0e-10 #one seed per total plant biomass of wetland
    reprod=0.05
    dispersalAlpha=20.0
    dispersalC=0.5
    LDDED=0.0
    LDDVD=0.0
    overWinter=0.3
    seedPerBiomass=1000*100*100
    TWidth=90.0
    z=collect(linspace(0,40,NoSpecies))
    sDist=0.
    C=connectivity(XY,dispersalAlpha,dispersalC)
    CCstart=200
    CCamp=5.0
    CCk=75.0
    tempSlope=1.0/100000.0
    tempGrad=20-((XY[2,:]/extent)*2-1.0)*tempSlope*extent
    inData="inData"
    outData="outData"
    p=par(NoSpecies,NoSites,NoLandscape,Tend,repl,NoNoise,Poisson,r,m,ext,alpha,reprod,dispersalAlpha,dispersalC,LDDED,LDDVD,overWinter,seedPerBiomass,extent,TWidth,z,XY,sDist,C,CCstart,CCamp,CCk,tempSlope,tempGrad,N,inData,outData)
end

function moments(X,p)
T=size(X,1)
S=size(X,2)
N=size(X,3)
M=zeros(Float64,T,N,5)
X[isnan(X)]=0.0
  for t=1:T
    for n=1:N
      M[t,n,1]=sum(X[t,:,n])
      M[t,n,2]=sum(p.z.*X[t,:,n]')./M[t,n,1]
      M[t,n,3]=sum(p.z.*(M[t,n,2]'-p.z).^2)./M[t,n,1]
      M[t,n,4]=sum(p.z.*(M[t,n,2]'-p.z).^3)./M[t,n,1]
      M[t,n,5]=sum(p.z.*(M[t,n,2]'-p.z).^4)./M[t,n,1]
    end
  end
  return M
end

function extractX(X,a::Int64)
  for i=1:size(X,1)
    for j=1:size(X,2)
      for k=1:size(X,3)
      end
    end
  end
end

function Base.show(io::IO, E::par) # funciton to display par type
  N=fieldnames(E)
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

function getResult(L,P,f)
      file=f*"/out"*string(L)*"_"*string(P)*".h5"
    X=h5open(file,"r") do file
       read(file,"X")
       end
      xcc=h5open(file,"r") do file
       read(file,"XCS")
       end
        Tactual=h5open(file,"r") do file
       read(file,"Tactual")
       end
          SDX=h5open(file,"r") do file
       read(file,"SDX")
       end
       IE=h5open(file,"r") do file
        read(file,"IE")
      end
      ISD=h5open(file,"r") do file
       read(file,"IE")
     end
  return X,xcc,Tactual,SDX,IE,ISD
end
