cd("/Users/raukhur/Documents/Github/Landscapes")
using Landscapes
using Gadfly
reload("Landscapes")
p=Landscapes.getP()
p.inData="/Users/raukhur/Google Drive/Landscapes/inData"
p.outData="/Users/raukhur/Google Drive/Landscapes/outData"
p.NoLandscape=113
Landscapes.GetData(p)

plot(x=p.XY[1,:],y=p.XY[2,:])
out=Landscapes.Go(p)
pwd()
readdir()
p.tempGrad
