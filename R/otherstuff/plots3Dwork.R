## Started 10 Nov 2014 ##
## Lizzie, just flapping about trying to learn 3D plots ##

## I did this for the storage effect project ##
## it requires first running WithinYearDynamics.R ##

require(rgl)
require(plot3D)

# http://www.quantmod.com/examples/chartSeries3d/
# http://stackoverflow.com/questions/3979240/r-plotting-a-3d-surface-from-x-y-z

# Graph is the biomass in for species 1 (x), biomass in for species 2 (y) and then the out biomass for species 1 or 2 eventually (z)

plot3d(Bin[,1],Bin[,2],B1out[,,r])
surface3d(Bin[,1],Bin[,2],B1out[,,r])

material3d(col="white") # this does nothing... for me, so far
persp3d(Bin[,1],Bin[,2],B1out[,,r], col = "yellow", axes=FALSE, box=FALSE) #  aspect=c(1, 1, 0.5)

persp(Bin[,1],Bin[,2],B1out[,,r], theta = 40, phi = 40, col = "gold",
    border = NA, shade = 0.5)

# best ones here!
persp3D(z=B1out[,,r], x=Bin[,1], y=Bin[,2], theta=30, phi=30, expand=0.75,
    clab="biomass out")

require(plot3D)
quartz(width=9, height=4)
par(mfrow=c(1,2))
persp3D(z=B1out[,,r], x=Bin[,1], y=Bin[,2], theta=30, phi=30, expand=0.75,
    clab="biomass out Sp1", shade = 0.5, ticktype = "detailed", bty = "b2",
    xlab="inputs Sp1", ylab="input Sp2", zlab="biomass out")

persp3D(z=B2out[,,r], x=Bin[,1], y=Bin[,2], theta=120, phi=30, expand=0.75,
    clab="biomass out Sp2", shade = 0.5, ticktype = "detailed", bty = "b2",
    xlab="inputs Sp1", ylab="input Sp2", zlab="biomass out")

