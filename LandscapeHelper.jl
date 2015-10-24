function GetData(LandscapeNo, repl,NoiseSeries)
  file="inData/Landscapes/L"*string(LandscapeNo)*"/XY"*string(repl)*".h5"
  XY=h5open(file,"r") do file
       read(file,"XY")
       end
  Ext=h5open(file,"r") do file
       read(file,"Ext")
  end
    file="inData/Noise/N"*string(NoiseSeries)*".h5"
    N=h5open(file,"r") do file
       read(file,"N")
       end
  return XY,N,Ext
end

function makeXYLandscape(LandscapeNo::Int64,repl::Int64,NoSites::Int64,Ext::Float64,_save::Int64)
    if NoSites>0
        XY=rand(Float64,2,NoSites).*Ext #XY coordinates in a 0:1 range
    else
        NoSites=abs(NoSites)
        XY=[zeros(NoSites)'; collect(linspace(0,Ext,NoSites))']
        
    end
  if _save==1
        if isdir("inData/Landscapes/L"*string(LandscapeNo))
    else
      mkdir("inData/Landscapes/L"*string(LandscapeNo))
    end
    file="inData/Landscapes/L"*string(LandscapeNo)*"/XY"*string(repl)*".h5"
    A=h5open(file,"w") do file
      write(file,"XY", XY)
      write(file,"Ext", Ext)
       end

    for i=1:10
      file="inData/Landscapes/L"*string(LandscapeNo)*"/SiteVar"*string(i)*".h5"
      SV=rand(Float64,1,NoSites)*2.0-1.0
      A=h5open(file,"w") do file
        write(file,"SV", SV)
      end
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

function getP(NoSpecies,repl,NoSites,XY,D,extent,dispKernel,N)

    p=par(NoSpecies,NoSites,0.05,0.0001,0.9,0.01,dispKernel,0,0.0,0.7,100000.0,extent.0,90.0,linspace(0,40,NoSpecies),XY,0.0,D,200,4.0,75.0,1.0,20+((XY[2,:]/extent)*2-1.0)*5,N)

end


function getP(NoLandscape::Int64,repl::Int64,NoNoise::Int64,Poisson::Int64,alpha::Float64)

    NoSpecies=100
    XY,N,Ext=GetData(NoLandscape,repl,NoNoise)
    NoSites=size(XY,2)
    m=0.05
    ext=1.0e-06
    reprod=0.05
    dispersalAlpha=20.0
    dispersalC=0.5
    localTvar=0.0
    overWinter=0.3
    seedPerBiomass=1000*100*100
    extent=Ext
    TWidth=90.0
    z=linspace(0,40,NoSpecies)
    sDist=500/NoSites/2
    C=connectivity(XY,dispersalAlpha,dispersalC)
    CCstart=200
    CCamp=5.0
    CCk=75.0
    tempSlope=1.0/100000.0
    tempGrad=20-((XY[2,:]/Ext)*2-1.0)*tempSlope*Ext

    p=par(NoSpecies,NoSites,NoLandscape,repl,NoNoise,Poisson,m,ext,alpha,reprod,dispersalAlpha,dispersalC,localTvar,overWinter,seedPerBiomass,extent,TWidth,z,XY,sDist,C,CCstart,CCamp,CCk,tempSlope,tempGrad,N)

end

function moments(X,p)

T=size(X,1)
S=size(X,2)
N=size(X,3)
M=zeros(Float64,T,N,5)
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

